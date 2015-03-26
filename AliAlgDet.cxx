#include "AliAlgDet.h"
#include "AliAlgSens.h"
#include "AliAlgDet.h"
#include "AliESDtrack.h"
#include "AliAlgTrack.h"
#include "AliLog.h"
#include "AliGeomManager.h"

#include <TString.h>

ClassImp(AliAlgDet)


//____________________________________________
AliAlgDet::AliAlgDet()
  :fVolIDMin(-1)
  ,fVolIDMax(-1)
  ,fNSensors(0)
  ,fSID2VolID(0)
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
  ,fVolIDMin(-1)
  ,fVolIDMax(-1)
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
  const AliESDfriendTrack* trF(0);
  const AliTrackPointArray* trP(0);
  //
  ResetPool();
  //
  if (!(trF=esdTr->GetFriendTrack()) || !(trP=trF->GetTrackPointArray())) return 0;
  //
  int np(trP->GetNPoints());
  int npSel(0);
  AliAlgPoint* apnt(0);
  for (int ip=0;ip<np;ip++) {
    if (!(apnt=TrackPoint2AlgPoint(ip, trP))) continue;
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
  //
  // convert to detector tracking frame
  UShort_t vid = trpArr->GetVolumeID()[pntId];
  Int_t sid = VolID2SID(vid); // sensor index within the detector
  if (!sid) return 0;
  AliAlgPoint* pnt = GetPointFromPool();
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
void AliAlgDet::AddVolume(AliAlgVol* vol)
{
  // add volume
  if (GetVolume(vol->GetSymName())) {
    AliFatalF("Volume %s was already added to %s",vol->GetName(),GetName());
  }
  fVolumes.AddLast(vol);
  if (vol->IsSensor()) {
    fSensors.AddLast(vol);
    Int_t vid = ((AliAlgSens*)vol)->GetVolID();
    if (fVolIDMin<0 || vid<fVolIDMin) fVolIDMin = vid;
    if (fVolIDMax<0 || vid>fVolIDMax) fVolIDMax = vid;
  }
  //
}

//_________________________________________________________
void AliAlgDet::DefineMatrices()
{
  // define transformation matrices. Detectors may override this method
  //
  TGeoHMatrix mtmp;
  //
  TIter next(&fVolumes);
  AliAlgVol* vol(0);
  while ( (vol=(AliAlgVol*)next()) ) {
    // modified global-local matrix
    const TGeoHMatrix* g2l = AliGeomManager::GetMatrix(vol->GetSymName());
    if (!g2l) AliFatalF("Failed to find G2L matrix for %s",vol->GetSymName());
    vol->SetMatrixG2L(*g2l);
    //
    // ideal global-local matrix
    if (!AliGeomManager::GetOrigGlobalMatrix(vol->GetSymName(),mtmp)) 
      AliFatalF("Failed to find ideal G2L matrix for %s",vol->GetSymName());
    vol->SetMatrixG2LIdeal(mtmp);
    //
    if (vol->IsSensor()) { // tracking-local matrix
      AliAlgSens* sens = (AliAlgSens*)vol;
      const TGeoHMatrix* t2l = AliGeomManager::GetTracking2LocalMatrix(sens->GetVolID());
      sens->SetMatrixT2L(*t2l);
    }
  }
  //
}

//_________________________________________________________
void AliAlgDet::SortSensors()
{
  // build local tables for internal numbering
  fNSensors = fSensors.GetEntriesFast();
  if (!fNSensors) {
    AliWarning("No sensors defined");
    return;
  }
  fSensors.Sort();
  fSID2VolID = new Int_t[fNSensors]; // cash id's for fast binary search
  for (int i=0;i<fNSensors;i++) fSID2VolID[i] = GetSensor(i)->GetVolID();
  //
}

//_________________________________________________________
void AliAlgDet::Init()
{
  // define hiearchy, initialize matrices
  DefineVolumes();
  SortSensors();    // VolID's must be in increasing order
  DefineMatrices();
}

//_________________________________________________________
Int_t AliAlgDet::VolID2SID(Int_t vid) const 
{
  // find SID corresponding to VolID
  int mn(0),mx(fNSensors-1);
  while (mx>=mn) {
    int md( (mx+mn)>>1 ), vids(GetSensor(md)->GetVolID());
    if (vid<vids)      mx = md-1;
    else if (vid>vids) mn = md+1;
    else return md;
  }
  return -1;
}

//____________________________________________
void AliAlgDet::Print(const Option_t *opt) const
{
  // print info
  TString opts = opt;
  opts.ToLower();
  printf("Detector:%s %d volumes %d sensors {VolID: %d-%d}\n",
	 GetName(),GetNVolumes(),GetNSensors(),GetVolIDMin(),GetVolIDMax());
  if (opts.Contains("long")) for (int iv=0;iv<GetNVolumes();iv++) GetVolume(iv)->Print(opt);
  //
}
