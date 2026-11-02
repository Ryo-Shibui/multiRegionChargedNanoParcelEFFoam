# Gas transport references for the 320 K dilute H2O/N2 tutorial

These coefficients are project fits to published zero-density correlations, not Sutherland constants printed in the cited papers. Reproduce them with `python3 fit_gas_transport.py` (NumPy/SciPy). The CSV contains generated source-curve values and fitted values; the JSON records coefficients and error metrics. These small generated files and this note can be copied into a repository. Do not publish the downloaded source PDFs or extracted full texts with the tutorial.

## Viscosity source and fit

Lemmon and Jacobsen (2004), DOI 10.1023/B:IJOT.0000022327.04529.f3, supplies the dilute-gas viscosity expression in Eq.2 and Tables I-II (printed pp23-24). Its collision-integral parameters are reproduced in the script. Source uncertainties quoted in Section5 are 0.5% for dilute N2 and 1% for dilute O2 above200K and air. Open primary PDF: https://trc.nist.gov/refprop/FAQ/NAO.PDF

Huber et al. (2009), DOI 10.1063/1.3088050, supplies H2O zero-density viscosity in Eq.11 and Table2 (printedp109). We use that expression directly. It describes gaseous H2O at zero density, including temperatures below the normal boiling point. Evaluating a pure-water equilibrium state at1atm below373K would instead select liquid water and would be inappropriate for a dilute vapor species. Open primary PDF: https://srd.nist.gov/jpcrdreprint/1.3088050.pdf

Fit model: mu[Pa.s] = As*T[K]^1.5/(T[K]+Ts). As has units Pa.s.K^-0.5; Ts is in K. The objective is the sum of squared relative residuals on301 points at1K spacing,300-600K inclusive. A separate30001-point uniform grid determines maximum and RMS fit errors. Fit error is additional to source uncertainty and is not an experimental uncertainty.

|Species|As|Ts [K]|Maximum absolute error [%]|Error at320K [%]|
|---|---:|---:|---:|---:|
|N2|1.4525113378005044e-6|123.02751260895496|0.1986075|-0.1018094|
|O2|1.7476228876785538e-6|141.03939217022997|0.1982735|-0.1006461|
|H2O|2.021933233849032e-6|801.806493269838|2.3843763|-1.0943465|
|air (optional)|1.5193549366605398e-6|127.05943716472858|0.1974833|-0.1009610|

The H2O Sutherland approximation is less accurate than the N2/O2 fits; retain its2.385% fit-error bound prominently. These fits are validated only over300-600K and for continuum dilute-gas use, not dense fluids or liquids. They are suitable species inputs for a low-pressure gas-mixture transport model; the mixture rule remains a separate model choice.

At320K the old shared As=1.67212e-6,Ts=170.672 curve differs from each source by +3.8353% (N2), -10.1891% (O2), +87.0029% (gaseousH2O), and +0.1713% (air). The shared curve therefore primarily resembles air and substantially overstates dilute water-vapor viscosity.

## Experimental-compilation diffusion reference

Marrero and Mason (1972), DOI10.1063/1.3253094, Table13 (printedp35) and Figure76 (p83), gives specifically N2-H2O:

p[atm]*D[cm2/s] = 1.87e-6*T[K]^2.072, 282 <= T <=373K.

Equivalently D[m2/s] = 1.87e-10*T[K]^2.072*(101325/p[Pa]). Table11(p34) assigns an estimated uncertainty of +/-4%. Comments onp82 explain that only evaporation-tube measurements were available; O'Connell et al.(1969), with water diffusing downward through N2, were judged most reliable. The compilation did not apply an equimolar-composition correction for this group because composition effects were smaller than measurement uncertainty. We use its binary gas coefficient for the dilute tutorial, with this limitation stated. Open primary PDF: https://srd.nist.gov/jpcrdreprint/1.3253094.pdf

At101325Pa, generated values are D298.15=2.5053623144e-5,D300=2.5376799380e-5,D320=2.9007637639e-5 m2/s. For ideal N2 at320K, the fitted viscosity gives nu=1.7591962449e-5 m2/s. Thus Sc=nu/D=0.60645967; using Sc=1 would reduce D by39.35%, while the reference D is1.6489 times the Sc=1 value. This difference should be clear in a transient cosine-mode decay comparison. Trace Y_H2O=1e-4 +/-2e-5 supports the dilute approximation.

## Independent modern check

Burgess,NISTTN2279(2024), DOI10.6028/NIST.TN.2279, Table1a,p8, gives H2O/N2 from Hellmann(2019), DOI10.1021/acs.jced.9b00822: ln(D/cm2/s)=-10.324-84.181/T+1.620ln(T),150-2000K, at101325Pa. This is a NIST fit to first-principles calculations, not a new experiment. At320K it gives2.8869468234e-5m2/s, only0.48% below the experimental-compilation law. Its separately fitted room-temperature tabulatedDref=2.54e-5m2/s need not equal the broad-range expression exactly. NIST explains its table's2-3% uncertainty as experimental-comparison uncertainty, rather than inherent numerical or fit uncertainty. Open official PDF: https://nvlpubs.nist.gov/nistpubs/TechnicalNotes/NIST.TN.2279.pdf

No university authentication, institutional subscription session, or nonpublic source was used.
