#ifndef ALIALGDET_H
#define ALIALGDET_H

#include <TNamed.h>
#include <TObjArray.h>
class AliESDtrack;
class AliAlgTrack;
class AliAlgPoint;
class AliAlgSens;
class AliAlgVol;
class AliTrackPointArray;

class AliAlgDet : public TNamed
{
 public:
  AliAlgDet();
  AliAlgDet(const char* name, const char* title="");
  virtual ~AliAlgDet();
  //
  Int_t   GetVolIDMin()              const {return fVolIDMin;}
  Int_t   GetVolIDMax()              const {return fVolIDMax;}

  void    SetVolIDMin(Int_t v)             {fVolIDMin = v;}
  void    SetVolIDMax(Int_t v)             {fVolIDMax = v;}
  //
  void    AcknowledgeNewRun(Int_t run);
  //
  Int_t   VID2SID(Int_t vid)         const {return 0;} //todo
  Int_t   SID2VID(Int_t sid)         const {return 0;} //todo
  AliAlgSens* GetSensor(Int_t id)    const {return (AliAlgSens*)fSensors.UncheckedAt(id);}
  AliAlgVol*  GetVolume(Int_t id)    const {return (AliAlgVol*)fVolumes.UncheckedAt(id);}
  Int_t   GetNSensors()              const {return fSensors.GetEntriesFast();}
  Int_t   GetNVolumes()              const {return fVolumes.GetEntriesFast();}
  //
  Bool_t  VIDofDetector(Int_t id)    const {return id>=fVolIDMin && id<=fVolIDMax;}
  //
  virtual void AddVolume(AliAlgVol* vol);
  virtual void DefineVolumes();
  virtual void DefineMatrices();
  virtual void PrintHierarchy();
  virtual void ExtractSensorMatrices();
  virtual Int_t ProcessPoints(const AliESDtrack* esdTr, AliAlgTrack* algTrack);
  virtual AliAlgPoint* TrackPoint2AlgPoint(int pntId, const AliTrackPointArray* trp);
  //
  virtual AliAlgPoint* GetPointFromPool();
  void    ResetPool();
  UShort_t GetVolumeIDFromSymname(const Char_t *symname);
  //
 protected:
  
  Int_t     fVolIDMin;                   // min volID for this detector
  Int_t     fVolIDMax;                   // max volID for this detector
  Int_t     fNSensors;                   // number of sensors (i.e. volID's)
  Int_t*    fVID2SID;                    //[fNSensors] table of conversion from VID to sid
  //
  Int_t     fPoolNPoints;            // number of points in the pool
  Int_t     fPoolFreePointID;        // id of the last free point in the pool
  TObjArray fPointsPool;             // pool of aligment points
  TObjArray fSensors;                // all sensors of the detector
  TObjArray fVolumes;                // all volumes of the detector  
  ClassDef(AliAlgDet,1);             // base class for detector global alignment
};


#endif
