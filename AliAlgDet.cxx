#include "AliAlgDet.h"
#include "AliAlgSens.h"
#include "AliAlgDet.h"
#include "AliESDtrack.h"
#include "AliAlgTrack.h"
#include "AliLog.h"
#include "AliGeomManager.h"

ClassImp(AliAlgDet)


//____________________________________________
AliAlgDet::AliAlgDet()
  :fVolIDMin(0)
  ,fVolIDMax(0)
  //
  ,fPoolNPoints(0)
  ,fPoolFreePointID(0)
  ,fPointsPool()
  ,fSensors()
  ,fVolumes()
{
  // def c-tor
}

//____________________________________________
AliAlgDet::AliAlgDet(const char* name, const char* title) 
  : TNamed(name,title)
  ,fVolIDMin(0)
  ,fVolIDMax(0)
  //
  ,fPoolNPoints(0)
  ,fPoolFreePointID(0)
  ,fPointsPool()
  ,fSensors()
  ,fVolumes()
{
  // def c-tor
  
}

//____________________________________________
AliAlgDet::~AliAlgDet()
{
  // d-tor
  fSensors.Clear(); // sensors are also attached as volumes, don't delete them here
  fVolumes.Delete(); // here all is deleted
  fPointsPool.Delete();
}


//____________________________________________
Int_t AliAlgDet::ProcessPoints(const AliESDtrack* esdTr, AliAlgTrack* algTrack)
{
  // extract the points corresponding to this detector, recalibrate/realign them to the
  // level of the "starting point" for the alignment/calibration session
  const AliESDfriendTrack* trF=0;
  const AliTrackPointArray* trP=0;
  //
  ResetPool();
  //
  if (!(trF=esdTr->GetFriendTrack()) || !(trP=trF->GetTrackPointArray())) return 0;
  //
  const UShort_t* vidArr = trP->GetVolumeID();
  int np = trP->GetNPoints();  
  int npSel = 0;
  for (int ip=0;ip<np;ip++) {
    if ( !VIDofDetector(vidArr[ip]) ) continue;
    AliAlgPoint* apnt = TrackPoint2AlgPoint(ip, trP);
    algTrack->AddPoint(apnt);
    npSel++;
  }
  //
  return npSel;
}


//____________________________________________
AliAlgPoint* AliAlgDet::TrackPoint2AlgPoint(int pntId, const AliTrackPointArray* trpArr)
{
  // convert the pntId-th point to AliAlgPoint, detectors may override this method
  AliAlgPoint* pnt = GetPointFromPool();
  //
  // convert to detector tracking frame
  UShort_t vid = trpArr->GetVolumeID()[pntId];
  //
  Int_t sid = VID2SID(vid); // sensor index within the detector
  if (!sid) return 0;
  //
  double tra[3],loc[3],glo[3] = {trpArr->GetX()[pntId], trpArr->GetY()[pntId], trpArr->GetY()[pntId]};
  AliAlgSens* sens = GetSensor(sid);
  const TGeoHMatrix& matG2L = sens->GetMatrixG2L(); // local to global matrix
  matG2L.MasterToLocal(glo,loc);
  const TGeoHMatrix& matT2L = sens->GetMatrixT2L();  // matrix for tracking to local frame translation
  matT2L.MasterToLocal(loc,tra);
  //
  pnt->SetXYZTracking(tra[0],tra[1],tra[2]);
  //
  // todo
  return pnt;
  //
}

//_________________________________________________________
void AliAlgDet::AcknowledgeNewRun(Int_t run)
{
  // update parameters needed to process this run
}

//_________________________________________________________
void AliAlgDet::ExtractSensorMatrices()
{
  // extract matrices for sensors
  for (int sid=0;sid<GetNSensors();sid++) {
    int volID = SID2VID(sid);
    TGeoHMatrix *mcurr = new TGeoHMatrix();
    // to do
  }  
}

//_________________________________________________________
AliAlgPoint* AliAlgDet::GetPointFromPool()
{
  // fetch or create new free point from the pool.
  // detector may override this method to create its own points derived from AliAlgPoint
  //
  if (fPoolFreePointID>=fPoolNPoints) { // expand pool
    fPointsPool.AddAtAndExpand(new AliAlgPoint(), fPoolNPoints++);
  }
  //
  AliAlgPoint* pnt = (AliAlgPoint*) fPointsPool.UncheckedAt(fPoolFreePointID++);
  pnt->Clear();
  return pnt;
  //
}

//_________________________________________________________
void AliAlgDet::ResetPool()
{
  // declare pool free
  fPoolFreePointID = 0;
}
 
//_________________________________________________________
void AliAlgDet::DefineVolumes()
{
  // dummy method
  AliError("This method must be implemented by specific detector");
}

//_________________________________________________________
void AliAlgDet::PrintHierarchy()
{
  // dummy method
  AliError("This method must be implemented by specific detector");
}

//_________________________________________________________
void AliAlgDet::AddVolume(AliAlgVol* vol)
{
  // add volume
  if (fVolumes.FindObject(vol->GetName())) {
    AliFatalF("Volume %s was already added to %s",vol->GetName(),GetName());
  }
  fVolumes.AddLast(vol);
  if (vol->IsSensor()) fSensors.AddLast(vol);
  //
}

//_________________________________________________________
UShort_t AliAlgDet::GetVolumeIDFromSymname(const Char_t *symname) 
{
  // volume ID from symname
  if (!symname) return 0;
  //
  for (UShort_t vid=fVolIDMin;vid<fVolIDMax;vid++) {
    Int_t modId;
    AliGeomManager::ELayerID layerId = AliGeomManager::VolUIDToLayer(vid,modId);
    if (layerId>0 && layerId<=AliGeomManager::kLastLayer && 
	modId>=0 && modId<AliGeomManager::LayerSize(layerId)) {
      if (!strcmp(symname,AliGeomManager::SymName(layerId,modId))) return vid;
    }
  }
  return 0;
}

