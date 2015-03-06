/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliAlgSteer.h"
#include "AliAlgDet.h"
#include "AliLog.h"
#include <TGeoMatrix.h>

const char* AliAlgSteer::fgkDetectorName[AliAlgSteer::kNDetectors] = {"ITS", "TPC", "TRD", "TOF", "HMPID" };
const int   AliAlgSteer::fgkSkipLayers[AliAlgSteer::kNLrSkip] = {AliGeomManager::kPHOS1,AliGeomManager::kPHOS2,
								 AliGeomManager::kMUON,AliGeomManager::kEMCAL};

ClassImp(AliAlgSteer)

//________________________________________________________________
AliAlgSteer::AliAlgSteer()
  :fNDet(0)
  ,fRunNumber(-1)
  ,fAlgTrack(0)

{
  // def c-tor
  for (int i=kNDetectors;i--;) {
    fDetectors[i] = 0;
    fDetVolMinMax[i][0] = fDetVolMinMax[i][1] = 0;
    fDetPos[i] = -1;
  }
}

//________________________________________________________________
AliAlgSteer::~AliAlgSteer()
{
  // d-tor
  delete fAlgTrack;
  for (int i=0;i<fNDet;i++) delete fDetectors[i];
  for (int i=0;i<kNGeoms;i++) {
    delete fMatrixT2G[i];
    delete fMatrixT2L[i];
  }
}

//________________________________________________________________
void AliAlgSteer::LocalInit()
{
  // init all structures
  //
  static Bool_t done = kFALSE;
  if (done) return;
  done = kTRUE;
  //
  fAlgTrack = new AliAlgTrack();

  //
  int nLrSkip = sizeof(fgkSkipLayers)/sizeof(int);
  //
  // fill start and end of matrices for each used layer
  int nmat = 0;
  for (int i=0;i<AliGeomManager::kLastLayer;i++) {
    fMatrixID[kMatStart][i] = nmat;
    fMatrixID[kNMat][i] = 0;
    if (i<AliGeomManager::kFirstLayer) continue;
    Bool_t skip = kFALSE;
    for (int j=nLrSkip;j--;) if (fgkSkipLayers[j]==i) {skip=kTRUE; break;}
    if (skip) continue;
    nmat += fMatrixID[kNMat][i] = AliGeomManager::LayerSize(i);
  }
  for (int i=0;i<kNGeoms;i++) {
    fMatrixT2L[i] = new TClonesArray("TGeoHMatrix",nmat); // book the space for matrices
    fMatrixT2G[i] = new TClonesArray("TGeoHMatrix",nmat); // book the space for matrices
  }
  //
}

//________________________________________________________________
void AliAlgSteer::LoadMatrices(Int_t geomTyp)
{
  // fetch matrices for given geom type
  if (geomTyp<0||geomTyp>=kNGeoms) return;
  if (!AliGeomManager::GetGeometry()) AliFatal("No geometry is loaded");
  //
  TClonesArray &arrT2G  = *fMatrixT2G[geomTyp];
  TClonesArray &arrT2L  = *fMatrixT2L[geomTyp];
  arrT2G.Clear();
  arrT2L.Clear();
  //
  TGeoHMatrix mtt;
  //
  for (int il=AliGeomManager::kFirstLayer;il<AliGeomManager::kLastLayer;il++) {
    int nmat = fMatrixID[kNMat][il];
    AliInfo(Form("Loading %4d matrices%d for Layer: %s",nmat,geomTyp,AliGeomManager::LayerName(il)));
    if (!nmat) continue; // skip layer
    for (int imd=0;imd<nmat;imd++) {
      int vid = AliGeomManager::LayerToVolUID(il,imd);
      const TGeoHMatrix *matL2G=0,*matT2L=0;
      // attention: some layer's matrices are somewhat special
      switch (il) {
	/*
      case AliGeomManager::kSPD1 : 
      case AliGeomManager::kSPD2 : 
      case AliGeomManager::kSDD1 : 
      case AliGeomManager::kSDD2 : 
      case AliGeomManager::kSSD1 : 
      case AliGeomManager::kSSD2 : 
	mat = GetITSSensVolMatrix(vid); 
	break;
	//
	*/
      default:	
	matL2G = AliGeomManager::GetMatrix(vid); // g2l
	matT2L = AliGeomManager::GetTracking2LocalMatrix(vid); // t2l
      }
      if (!matL2G || !matT2L) {AliDebug(1,Form("No matrix for module %d",vid)); continue;}
      mtt = *matL2G;
      mtt.Multiply(matT2L);
      int id = GetMatrixID(il,imd);
      if (id<0) AliFatal(Form("MatrixID for VolID=%d should have been defined",vid));
      printf("Add %s at %d\n", AliGeomManager::SymName(vid) ,id);
      TGeoHMatrix* m0 = new( arrT2G[id] ) TGeoHMatrix(mtt);
      TGeoHMatrix* m1 = new( arrT2L[id] ) TGeoHMatrix(*matT2L);
      m0->SetName(AliGeomManager::SymName(vid));
      m1->SetName(AliGeomManager::SymName(vid));
      //
    }
  }
  //
}


//________________________________________________________________
void AliAlgSteer::AddDetector(const char* name)
{
  // add detector participating in the alignment
  int id = -1;
  TString names = name;
  names.ToUpper();
  for (int i=kNDetectors;i--;) if (names==fgkDetectorName[i]) {id=i;break;}
  if (id<0) AliFatal(Form("Detector %s is not known to alignment framework",name));
  AliAlgDet* det = 0;
  switch(id) {
  case kITS: det = new AliAlgITS(); break;
    //  case kTPC: det = new AliAlgTPC(); break;
    //  case kTRD: det = new AliAlgTRD(); break;
    //  case kTOF: det = new AliAlgTOF(); break;
  default: AliErrorF("%s not implemented yet",name); break;
  };
  //
  fDetectors[fNDet] = det;
  fDetVolMinMax[fNDet][0] = det->GetVolIDMin();
  fDetVolMinMax[fNDet][1] = det->GetVolIDMax();
  fDetPos[id] = fNDet;
  //
  fNDet++;
  //
}

//________________________________________________________________
/*
TGeoHMatrix* AliAlgSteer::GetITSSensVolMatrix(Int_t vid)
{
  // special matrix extraction for its
  static TGeoHMatrix mat;
  AliGeomManager::ELayerID lay = AliGeomManager::VolUIDToLayer(vid);
  if (lay<1|| lay>6) return 0;
  Int_t idx=Int_t(vid)-2048*lay;
  if (idx>=AliGeomManager::LayerSize(lay)) return -1;
  for (Int_t ilay=1; ilay<lay; ilay++) idx += AliGeomManager::LayerSize(ilay);
  Double_t rot[9];
  if (!AliITSgeomTGeo::GetRotation(idx,rot)) return 0;
  mat.SetRotation(rot);
  Double_t oLoc[3]={0,0,0}, oGlo[3]={0,0,0};
  if (!AliITSgeomTGeo::LocalToGlobal(idx,oLoc,oGlo)) return -3;
  mat.SetTranslation(oGlo);
  return &mat;
}
*/

//_________________________________________________________
Bool_t AliAlgSteer::AcceptTrack(const AliESDtrack* esdTr) const
{
  // decide if the track should be processed
  return kTRUE;
  //
}

//_________________________________________________________
Bool_t AliAlgSteer::ProcessTrack(const AliESDtrack* esdTr)
{
  // process single track
  //
  int nPnt = 0;
  const AliESDfriendTrack* trF = esdTr->GetFriendTrack();
  if (!trF) return kFALSE;
  const AliTrackPointArray* trPoints = trF->GetTrackPointArray();
  if (!trPoints || (nPnt=trPoints->GetNPoints())<1) return kFALSE;
  //
  if (!AcceptTrack(esdTr)) return kFALSE;
  //
  fAlgTrack->Clear();
  //
  // process the track points for each detector, fill needed points in the fAlgTrack
  for (int idet=0;idet<fNDet;idet++) {
    GetDetector(idet)->ProcessPoints(esdTr, fAlgTrack);
  }
  //
  return kTRUE;
}

//_________________________________________________________
Bool_t AliAlgSteer::AddDetector(AliAlgDet* det)
{
  // add new detector to alignment
  if (!det) return kFALSE;
  if (fDetectors.FindObject(det->GetName())) {
    AliError(Form("Detector %s was already added",det->GetName()));
    return kFALSE;
  }
  fDetectors.AddLast(det);
  fNDet++;
  return kTRUE;
}

//_________________________________________________________
void AliAlgSteer::SetRunNumber(Int_t run)
{
  if (run==fRunNumber) return;  // nothing to do
  //
  AcknowledgeNewRun(run); 
}

//_________________________________________________________
void AliAlgSteer::AcknowledgeNewRun(Int_t run)
{
  // load needed info for new run
  fRunNumber = run;
  for (int idet=0;idet<fNDet;idet++) GetDetector(idet)->AcknowledgeNewRun(run);
  //
}
