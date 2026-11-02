#!/usr/bin/env python3
"""Reproduce dilute-gas source correlations and local Sutherland fits.
Requires numpy and scipy. Does not depend on CoolProp or any proprietary data.
All outputs are generated numeric data; As/Ts are OUR fits, not literature constants.
"""
from pathlib import Path
import csv,json
import numpy as np
from scipy.optimize import least_squares

ROOT=Path(__file__).resolve().parent
B=np.array([0.431,-0.4623,0.08406,0.005341,-0.00331])
LJ={"N2":(28.01348,98.94,0.3656),"O2":(31.9988,118.5,0.3428),"air":(28.9586,103.3,0.360)}
SOURCES={
 "Lemmon2004":{"doi":"10.1023/B:IJOT.0000022327.04529.f3","url":"https://trc.nist.gov/refprop/FAQ/NAO.PDF","location":"Eq.2, Tables I and II, printed pp23-24"},
 "Huber2009":{"doi":"10.1063/1.3088050","url":"https://srd.nist.gov/jpcrdreprint/1.3088050.pdf","location":"Eq.11 and Table2, printed p109"},
 "MarreroMason1972":{"doi":"10.1063/1.3253094","url":"https://srd.nist.gov/jpcrdreprint/1.3253094.pdf","location":"Table13 printed p35; H2O-N2 comments p82; Figure76 p83"},
 "Burgess2024":{"doi":"10.6028/NIST.TN.2279","url":"https://nvlpubs.nist.gov/nistpubs/TechnicalNotes/NIST.TN.2279.pdf","location":"Table1a printed p8; equation definition p7; reference25 Hellmann2019"}}

def mu0(species,T):
    T=np.asarray(T,dtype=float)
    if species=='H2O':
        theta=T/647.096
        H=[1.67752,2.20462,0.6366564,-0.241605]
        return 1e-6*100*np.sqrt(theta)/sum(h/theta**i for i,h in enumerate(H))
    M,eps,sigma=LJ[species]
    x=np.log(T/eps)
    omega=np.exp(sum(b*x**i for i,b in enumerate(B)))
    # Published Eq.2 gives microPa.s using M in g/mol and sigma in nm.
    return 1e-6*0.0266958*np.sqrt(M*T)/(sigma*sigma*omega)

def sutherland(T,As,Ts):
    return As*np.asarray(T)**1.5/(np.asarray(T)+Ts)

def diffusion_marrero(T,p=101325.0):
    # Table13 / Figure76: p[atm] D[cm2/s] = 1.87e-6 T[K]^2.072.
    return 1.87e-10*np.asarray(T)**2.072*(101325.0/p)

def diffusion_burgess(T,p=101325.0):
    T=np.asarray(T)
    return 1e-4*np.exp(-10.324-84.181/T+1.620*np.log(T))*(101325.0/p)

if __name__=='__main__':
    train=np.arange(300.0,601.0,1.0)
    validate=np.linspace(300.0,600.0,30001)
    results={'sources':SOURCES,'fit_description':{'model':'mu[Pa.s]=As*T[K]^1.5/(T[K]+Ts)',
      'As_unit':'Pa.s.K^-0.5','Ts_unit':'K','target':'source zero-density viscosity',
      'training':'301 uniformly spaced points at 1 K spacing from300 through600K',
      'objective':'unweighted sum of squared relative residuals (mu_fit/mu_source-1)',
      'validation':'30001 uniformly spaced points from300 through600K',
      'caveat':'fit error is relative to source correlation, not an experimental uncertainty; no validation outside fit interval'},'species':{}}
    rows=[]
    for species in ['N2','O2','H2O','air']:
        y=mu0(species,train)
        opt=least_squares(lambda q:sutherland(train,q[0]*1e-6,q[1])/y-1,[1.5,150.0],
            bounds=([0.01,-299.0],[100,100000]),xtol=1e-14,ftol=1e-14,gtol=1e-14,max_nfev=5000)
        As,Ts=float(opt.x[0]*1e-6),float(opt.x[1])
        err=sutherland(validate,As,Ts)/mu0(species,validate)-1
        old=sutherland(validate,1.67212e-6,170.672)/mu0(species,validate)-1
        out={'As':As,'Ts':Ts,'source':'Huber2009' if species=='H2O' else 'Lemmon2004',
             'max_abs_relative_error_pct':float(np.max(np.abs(err))*100),
             'rms_relative_error_pct':float(np.sqrt(np.mean(err**2))*100),
             'relative_error_pct_at_320K':float((sutherland(320,As,Ts)/mu0(species,320)-1)*100),
             'mu_source_320_Pa_s':float(mu0(species,320)),
             'old_shared_constants_error_320_pct':float((sutherland(320,1.67212e-6,170.672)/mu0(species,320)-1)*100),
             'old_shared_constants_max_error_pct':float(np.max(np.abs(old))*100)}
        results['species'][species]=out
        for T in train:
            rows.append({'species':species,'T_K':float(T),'mu_source_Pa_s':float(mu0(species,T)),
                'mu_sutherland_Pa_s':float(sutherland(T,As,Ts)),
                'relative_error_pct':float((sutherland(T,As,Ts)/mu0(species,T)-1)*100)})
    T=320.0; p=101325.0
    rho=p*0.02801348/(8.31446261815324*T)
    nu=float(sutherland(T,results['species']['N2']['As'],results['species']['N2']['Ts'])/rho)
    results['diffusion']={'recommended_tutorial_reference':'MarreroMason1972',
       'formula_SI':'D[m2/s]=1.87e-10*T[K]^2.072*(101325/p[Pa])',
       'valid_T_K':[282,373],'reference_pressure_Pa':101325,
       'Marrero_uncertainty_pct':4,'Marrero_uncertainty_location':'Table11, printedp34',
       'Marrero_D_298_15':float(diffusion_marrero(298.15)),
       'Marrero_D_300':float(diffusion_marrero(300)),
       'Marrero_D_320':float(diffusion_marrero(320)),
       'Burgess_D_298_15_from_extended_fit':float(diffusion_burgess(298.15)),
       'Burgess_D_320_from_extended_fit':float(diffusion_burgess(320)),
       'Burgess_tabulated_Dref_298_15':2.54e-5,
       'Burgess_caveat':'tabulated Dref is separately fitted near room temperature; extended fit need not exactly reproduce it',
       'rho_N2_ideal_320_kg_m3':rho,'nu_N2_fitted_320_m2_s':nu,
       'equivalent_Sc_320':nu/float(diffusion_marrero(320)),
       'D_reference_divided_by_D_Sc1_320':float(diffusion_marrero(320))/nu,
       'dilute_Y_assumption':'Y_H2O=1e-4 +/-2e-5 in N2; use binary dilute-gas Fick diffusion; ideal-gas pressure scaling'}
    (ROOT/'gas_transport_reference.json').write_text(json.dumps(results,indent=2)+'\n')
    with (ROOT/'gas_viscosity_source_and_fit.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    print(json.dumps(results,indent=2))
