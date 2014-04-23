#include "AliGloAlg.h"
#include "AliLog.h"
#include <TGeoMatrix.h>

const char* AliGloAlg::fgkDetectorName[AliGloAlg::kNDetectors] = {"ITS", "TPC", "TRD", "TOF", "HMPID" };
const int   AliGloAlg::fgkSkipLayers[AliGloAlg::kNLrSkip] = {AliGeomManager::kPHOS1,AliGeomManager::kPHOS2,AliGeomManager::kMUON,AliGeomManager::kEMCAL};

ClassImp(AliGloAlg)

//________________________________________________________________
AliGloAlg::AliGloAlg() :
fAlgTrack(0)
{
  // def c-tor
  //  for (int i=AliGloAlg::kNDetectors;i--;) fAlignerID[i] = 0;



}

//________________________________________________________________
AliGloAlg::~AliGloAlg()
{
  // d-tor
  delete fAlgTrack;

  for (int i=0;i<kNGeoms;i++) {
    delete fMatrixT2G[i];
    delete fMatrixT2L[i];
  }
}

//________________________________________________________________
void AliGloAlg::LocalInit()
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
void AliGloAlg::LoadMatrices(Int_t geomTyp)
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
void AliGloAlg::AddDetector(const char* name)
{
  // add detector participating in the alignment
  int id = -1;
  TString names = name;
  names.ToUpper();
  for (int i=kNDetectors;i--;) if (names==fgkDetectorName[i]) {id=i;break;}
  if (id<0) AliFatal(Form("Alignment of %s is not foreseen",name));
  
}

//________________________________________________________________
/*
TGeoHMatrix* AliGloAlg::GetITSSensVolMatrix(Int_t vid)
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
Bool_t AliGloAlg::AcceptTrack(const AliESDtrack* esdTr) const
{
  // decide if the track should be processed
  AliESDfriendTrack* trF = esdTr->GetFriendTrack();
  if (!trF) {
    AliError("No friend track found");
    return kFALSE;
  }
  const AliTrackPointArray* trPoints = trF->GetTrackPointArray();
  if (!trPoints || trPoints.GetNPoints()<1) return kFALSE;
  //
  // do other checks
  //
  return kTRUE;
  //
}

//_________________________________________________________
Bool_t AliGloAlg::ProcessTrack(const AliESDtrack* esdTr)
{
  // process single track
  //
  if (!AcceptTrack(esdTr)) return kFALSE;
  fAlgTrack->Clear();
  //
  UInt_t detStat = 0; // status of track processing by each detector
  // process the track points for each detector, fill needed points in the fAlgTrack
  for (int idet=0;idet<fNDet;idet++) if (GetDetector(idet)->ProcessTrack(esdTr, fAlgTrack)) detStat |= 0x1<<idet;
  
}

//_________________________________________________________
Bool_t AliGloAlg::AddDetector(AliAlgDetector* det)
{
  // add new detector to alignment
  if (!det) return kFALSE;
  if (fDetectors->FindObject(det->GetName())) {
    AliError(Form("Detector %d was already added",det->GetName()));
    return kFALSE;
  }
  fDetectors->AddLast(det);
  fNDet++;
  return kTRUE;
}
