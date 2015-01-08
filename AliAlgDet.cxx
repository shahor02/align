#include "AliAlgDet.h"
#include "AliESDtrack.h"
#include "AliAlgTrack.h"

ClassImp(AliAlgDet)


//____________________________________________
AliAlgDet::AliAlgDet()
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
Bool_t AliAlgDet::ProcessTrack(const AliESDtrack* esdTr, AliAlgTrack* fAlgTrack)
{
  // extract the points corresponding to this detector, recalibrate/realign them to the
  // level of the "starting point" for the alignment/calibration session
  const AliTrackPointArray* trp = esdTr->GetFriendTrack()->GetTrackPointArray();
  if (!trp) return kFALSE;
  //
  int np = trp->GetNPoints();
  int npSel = 0;
  for (int ip=0;ip<np;ip++) {
    int vID = trp->GetVolumeID()[ip];
    if (vID<fVolIDMin || vID>fVolIDMax) continue;
    npSel++;
  }
  //
  return kTRUE;
}

//_________________________________________________________
void AliAlgDet::AcknowledgeNewRun(Int_t run)
{
  // update parameters needed to process this run
}
