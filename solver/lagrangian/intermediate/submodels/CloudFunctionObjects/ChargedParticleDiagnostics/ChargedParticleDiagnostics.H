/* GPL-3.0-or-later. Native cloud callbacks, separate CSV files per rank. */
#ifndef ChargedParticleDiagnostics_H
#define ChargedParticleDiagnostics_H
#include "CloudFunctionObject.H"
#include "wallPolyPatch.H"
#include "OSspecific.H"
#include "Pstream.H"
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <vector>
namespace Foam
{
template<class CloudType>
class ChargedParticleDiagnostics : public CloudFunctionObject<CloudType>
{
    typedef typename CloudType::parcelType parcelType;
    struct Event
    {
        scalar time, d, charge;
        point position;
        label origProc, origId;
        word patch;
    };
    std::ofstream trajectory_, deposition_;
    scalar restartTime_;
    std::map<std::string, Event> pending_;
    std::set<std::string> deposited_;
    std::string key(const parcelType& p) const
    {
        return std::to_string(p.origProc()) + ":" + std::to_string(p.origId());
    }
    scalar eventTime(const parcelType& p) const
    {
        const Time& t = this->owner().mesh().time();
        return t.value() - t.deltaTValue() + p.stepFraction()*t.deltaTValue();
    }
    void reconcile(const fileName& path, const bool deposition)
    {
        if (!isFile(path)) return;
        std::ifstream input(path.c_str());
        std::vector<std::string> retained;
        std::string line;
        if (std::getline(input, line)) retained.push_back(line);
        while (std::getline(input, line))
        {
            std::istringstream row(line);
            std::string ts, proc, id;
            if (!std::getline(row,ts,',') || !std::getline(row,proc,',')
                || !std::getline(row,id,',')) continue;
            std::istringstream number(ts);
            double t;
            if (!(number >> t)) continue;
            // At t=0 the run starts afresh. Checkpoint state owns its history.
            if (restartTime_ > SMALL && t <= restartTime_ + SMALL)
            {
                retained.push_back(line);
                if (deposition) deposited_.insert(proc + ":" + id);
            }
        }
        input.close();
        std::ofstream output(path.c_str(), std::ios::trunc);
        if (!output) FatalErrorInFunction << "Cannot reconcile " << path << exit(FatalError);
        for (const auto& entry : retained) output << entry << '\n';
    }
    void openFiles()
    {
        if (trajectory_.is_open()) return;
        // A rank never writes another rank's stream.
        fileName dir(this->outputDir()/this->owner().mesh().name()/("processor" + Foam::name(Pstream::myProcNo())));
        mkDir(dir);
        const fileName tf(dir/"trajectory.csv"), df(dir/"deposition.csv");
        reconcile(tf, false);
        reconcile(df, true);
        const bool th = !isFile(tf) || fileSize(tf) == 0;
        const bool dh = !isFile(df) || fileSize(df) == 0;
        trajectory_.open(tf.c_str(), std::ios::app);
        deposition_.open(df.c_str(), std::ios::app);
        if (!trajectory_ || !deposition_)
        {
            FatalErrorInFunction << "Cannot open diagnostics in " << dir
                << exit(FatalError);
        }
        trajectory_ << std::setprecision(17);
        deposition_ << std::setprecision(17);
        const char* header = "time,origProc,origId,x,y,z,d,chargeNumber,patch\n";
        if (th) trajectory_ << header;
        if (dh) deposition_ << header;
    }
    Event event(const parcelType& p, const word& patch) const
    {
        return {eventTime(p), p.d(), p.chargeNumber(), p.position(),
            p.origProc(), p.origId(), patch};
    }
    void row(std::ofstream& out, const Event& e)
    {
        out << e.time << ',' << e.origProc << ',' << e.origId << ','
            << e.position.x() << ',' << e.position.y() << ',' << e.position.z()
            << ',' << e.d << ',' << e.charge << ',' << e.patch << '\n';
    }
    void confirm(const parcelType& p)
    {
        const auto id = key(p);
        auto iter = pending_.find(id);
        // postPatch runs BEFORE localInteraction. Confirm only the retained,
        // inactive parcel AFTER correct() has implemented stick.
        if (iter != pending_.end())
        {
            if (!p.active() && deposited_.insert(id).second)
            {
                openFiles();
                row(deposition_, iter->second);
            }
            pending_.erase(iter);
        }
    }
public:
    TypeName("chargedParticleDiagnostics");
    ChargedParticleDiagnostics
    (
        const dictionary& dict, CloudType& owner, const word& modelName
    )
    : CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
      restartTime_(owner.mesh().time().value()) {}
    ChargedParticleDiagnostics(const ChargedParticleDiagnostics& f)
    : CloudFunctionObject<CloudType>(f), restartTime_(f.restartTime_), deposited_(f.deposited_) {}
    virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
    {
        return autoPtr<CloudFunctionObject<CloudType>>
            (new ChargedParticleDiagnostics(*this));
    }
    virtual bool postMove
    (
        parcelType& p, const scalar, const point&,
        const typename parcelType::trackingData& td
    )
    {
        confirm(p);
        openFiles();
        if (td.depositedParticleContact && deposited_.insert(key(p)).second)
        {
            row(deposition_, event(p, "depositedParticle"));
        }
        row(trajectory_, event(p, word::null));
        return true;
    }
    virtual bool postPatch
    (
        const parcelType& p, const polyPatch& pp,
        const typename parcelType::trackingData&
    )
    {
        // Outlet/processor patches and rebounds cannot count as deposition.
        if (isA<wallPolyPatch>(pp) && p.active() && !deposited_.count(key(p)))
        {
            pending_[key(p)] = event(p, pp.name());
        }
        return true;
    }
    virtual void postEvolve(const typename parcelType::trackingData&)
    {
        // A hit at stepFraction==1 has no following postMove callback.
        for (const auto& p : this->owner()) confirm(p);
        pending_.clear();
        openFiles();
        trajectory_.flush();
        deposition_.flush();
    }
private:
    virtual void write() {}
};
}
#endif
