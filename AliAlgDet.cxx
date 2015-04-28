#include "AliAlgDet.h"
#include "AliAlgSens.h"
#include "AliAlgDet.h"
#include "AliAlgSteer.h"
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
  ,fObligatory(kTRUE)
  ,fTrackFlagSel(0)
  ,fNPointsSel(0)
  //
  ,fSensors()
  ,fVolumes()
  //
  ,fNPoints(0)
  ,fPoolNPoints(0)
  ,fPoolFreePointID(0)
  ,fPointsPool()
  ,fAlgSteer(0)
{
  // def c-tor
  SetUniqueID(AliAlgSteer::kUndefined); // derived detectors must override this
  fAddError[0] = fAddError[1] = 0;
}

//____________________________________________
AliAlgDet::AliAlgDet(const char* name, const char* title) 
  : TNamed(name,title)
  ,fVolIDMin(-1)
  ,fVolIDMax(-1)
  ,fNSensors(0)
  ,fSID2VolID(0)
  //
  ,fSensors()
  ,fVolumes()
  //
  ,fNPoints(0)
  ,fPoolNPoints(0)
  ,fPoolFreePointID(0)
  ,fPointsPool()
  ,fAlgSteer(0)
{
  // def c-tor
  SetUniqueID(AliAlgSteer::kUndefined); // derived detectors must override this  
  fAddError[0] = fAddError[1] = 0;
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
Int_t AliAlgDet::ProcessPoints(const AliESDtrack* esdTr, AliAlgTrack* algTrack, Bool_t inv)
{
  // Extract the points corresponding to this detector, recalibrate/realign them to the
  // level of the "starting point" for the alignment/calibration session.
  // If inv==true, the track propagates in direction of decreasing tracking X 
  // (i.e. upper leg of cosmic track)
  //
  const AliESDfriendTrack* trF(esdTr->GetFriendTrack());
  const AliTrackPointArray* trP(trF->GetTrackPointArray());
  //
  int np(trP->GetNPoints());
  int npSel(0);
  AliAlgPoint* apnt(0);
  for (int ip=0;ip<np;ip++) {
    if (!SensorOfDetector(trP->GetVolumeID()[ip])) continue;
    if (!(apnt=TrackPoint2AlgPoint(ip, trP))) continue;
    algTrack->AddPoint(apnt);
    if (inv) apnt->SetInvDir();
    npSel++;
    fNPoints++;
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
  if (!sid<0) return 0;
  AliAlgSens* sens = GetSensor(sid);
  if (sens->GetSkip()) return 0;
  AliAlgPoint* pnt = GetPointFromPool();
  pnt->SetSensor(sens);
  //
  double tra[3],loc[3],glo[3] = {trpArr->GetX()[pntId], trpArr->GetY()[pntId], trpArr->GetZ()[pntId]};
  //  const TGeoHMatrix& matL2G = sens->GetMatrixL2G(); // local to global matrix
  const TGeoHMatrix& matL2Gor = sens->GetMatrixL2GOrig(); // local to global orig matrix
  matL2Gor.MasterToLocal(glo,loc);
  const TGeoHMatrix& matT2L = sens->GetMatrixT2L();  // matrix for tracking to local frame translation
  matT2L.MasterToLocal(loc,tra);
  //
  // convert error
  TGeoHMatrix hcov;
  Double_t hcovel[9];
  const Float_t *pntcov = trpArr->GetCov()+pntId*6; // 6 elements per error matrix
  hcovel[0] = double(pntcov[0]);
  hcovel[1] = double(pntcov[1]);
  hcovel[2] = double(pntcov[2]);
  hcovel[3] = double(pntcov[1]);
  hcovel[4] = double(pntcov[3]);
  hcovel[5] = double(pntcov[4]);
  hcovel[6] = double(pntcov[2]);
  hcovel[7] = double(pntcov[4]);
  hcovel[8] = double(pntcov[5]);
  hcov.SetRotation(hcovel);
  hcov.Multiply(&matL2Gor);
  hcov.MultiplyLeft(&matL2Gor.Inverse());
  hcov.Multiply(&matT2L);
  hcov.MultiplyLeft(&matT2L.Inverse());
  //
  Double_t *hcovscl = hcov.GetRotationMatrix();
  const double *sysE = sens->GetAddError(); // additional syst error
  pnt->SetYZErrTracking(hcovscl[4]+sysE[0]*sysE[0],hcovscl[5],hcovscl[8]+sysE[1]*sysE[1]);
  pnt->SetXYZTracking(tra[0],tra[1],tra[2]);
  pnt->SetAlphaSens(sens->GetAlpTracking());
  pnt->SetXSens(sens->GetXTracking());
  pnt->SetDetID(GetDetID());
  pnt->SetSID(sid);
  //
  pnt->SetContainsMeasurement();
  //
  pnt->Init();
  //
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
  fNPoints = 0;
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
    vol->PrepareMatrixL2G();
    // ideal global-local matrix
    vol->PrepareMatrixL2GOrig();
    //
  }
  // Now set tracking-local matrix (MUST be done after ALL L2G matrices are done!)
  // Attention: for sensor it is a real tracking matrix extracted from
  // the geometry but for container alignable volumes the tracking frame
  // is used for as the reference for the alignment parameters only,
  // see its definition in the AliAlgVol::PrepateMatrixT2L
  next.Reset();
  while ( (vol=(AliAlgVol*)next()) ) vol->PrepareMatrixT2L();
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
void AliAlgDet::InitGeom()
{
  // define hiearchy, initialize matrices
  if (GetInitGeomDone()) return;
  //
  DefineVolumes();
  SortSensors();    // VolID's must be in increasing order
  DefineMatrices();
  //
  SetInitGeomDone();
}

//_________________________________________________________
void AliAlgDet::InitDOFs()
{
  // initialize free parameters
  if (GetInitDOFsDone()) return;
  //
  int gloCount0(fAlgSteer->GetNDOFs()), gloCount(fAlgSteer->GetNDOFs());
  int nvol = GetNVolumes();
  for (int iv=0;iv<nvol;iv++) {
    AliAlgVol *vol = GetVolume(iv);
    // call init for root level volumes, they will take care of their children
    if (!vol->GetParent() && !vol->GetInitDOFsDone()) vol->InitDOFs(gloCount);
  }
  //
  fNDOFs = gloCount-gloCount0;
  //
  SetInitDOFsDone();
  return;
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
  printf("Detector:%5s %5d volumes %5d sensors {VolID: %5d-%5d} Def.Sys.Err: %.4e %.4e\n",
	 GetName(),GetNVolumes(),GetNSensors(),GetVolIDMin(),GetVolIDMax(), fAddError[0],fAddError[1]);
  if (opts.Contains("long")) for (int iv=0;iv<GetNVolumes();iv++) GetVolume(iv)->Print(opt);
  //
}

//____________________________________________
void AliAlgDet::SetDetID(UInt_t tp)
{
  // assign type
  if (tp>=AliAlgSteer::kNDetectors) AliFatalF("Detector typeID %d exceeds allowed range %d:%d",
					      tp,0,AliAlgSteer::kNDetectors-1);
  SetUniqueID(tp);
}

//____________________________________________
void AliAlgDet::SetAddError(double sigy, double sigz)
{
  // add syst error to all sensors
  AliInfoF("Adding sys.error %.4e %.4e to all sensors",sigy,sigz);
  fAddError[0] = sigy;
  fAddError[1] = sigz;
  for (int isn=GetNSensors();isn--;) GetSensor(isn)->SetAddError(sigy,sigz);
  //
}

//____________________________________________
void AliAlgDet::SetObligatory(Bool_t v)
{
  // mark detector presence obligatory in the track
  fObligatory = v;
  fAlgSteer->SetObligatoryDetector(GetDetID(),v);
}
