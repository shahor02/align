#include "AliAlgDet.h"
#include "AliESDtrack.h"
#include "AliAlgTrack.h"

ClassImp(AliAlgDet)


//____________________________________________
AliAlgDet::AliAlgDet()
  :fVolIDMin(0)
  ,fVolIDMax(0)
{
  // def c-tor
}

//____________________________________________
AliAlgDet::AliAlgDet(const char* name, const char* title) :
  TNamed(name,title)
{
  // def c-tor
  
}

//____________________________________________
Int_t AliAlgDet::ProcessPoints(const AliESDtrack* esdTr, AliAlgTrack* fAlgTrack)
{
  // extract the points corresponding to this detector, recalibrate/realign them to the
  // level of the "starting point" for the alignment/calibration session
  const AliESDfriendTrack* trF=0;
  const AliTrackPointArray* trp=0;
  if (!(trF=esdTr->GetFriendTrack()) || !(trp=trF->GetTrackPointArray())) return kFALSE;
  //
  UShort_t* vidArr = trPoints->GetVolumeID();
  int np = trp->GetNPoints();  
  int npSel = 0;
  for (int ip=0;ip<np;ip++) {
    if ( !VIDofDetector(vidArr[ip]) ) continue;
    npSel++;
  }
  //
  return npSel;
}


//____________________________________________
AliAlgPoint* AliAlgDet::TrackPoint2AlgPoint(int pntId, const AliTrackPoitArray* trpArr)
{
  // convert the pntId-th point to AliAlgPoint, detectors may override this method
  AliAlgPoint* pnt;
  //
  // convert to detector tracking frame
  UShort_t vid = trpArr->GetVolumeID()[pntID];
  //
  Int_t sid = GetSensorIndex(vid); // sensor index within the detector
  if (!sid) return 0;
  //
  TGeoHMatrix *sMatrix = GetSensorT2GMatrixBySID(sid);  // matrix for global - tracking frame translation
  double tra[3],glo[3] = {pnt->GetX()[pntId], pnt->GetY()[pntId], pnt->GetY()[pntId]};
  //
  

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
  for (int sid=0;i<fNSensors;sid++) {
    int volID = SID2VID(sid);
    TGeoHMatrix *mcurr = new TGeoHMatrix();
    
  }
  
}
