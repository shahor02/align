#include "AliAlgDetTRD.h"
#include "AliAlgVol.h"
#include "AliAlgSensTRD.h"
#include "AliAlgSteer.h"
#include "AliGeomManager.h"
#include "AliESDtrack.h"
#include "AliTRDgeometry.h"
#include "TGeoManager.h"

ClassImp(AliAlgDetTRD);

//____________________________________________
AliAlgDetTRD::AliAlgDetTRD()
{
  // default c-tor
  SetDetID(AliAlgSteer::kTRD);
}

//____________________________________________
AliAlgDetTRD::AliAlgDetTRD(const char* name, const char* title)
  :AliAlgDet(name,title)
{
  SetDetID(AliAlgSteer::kTRD);
}

//____________________________________________
AliAlgDetTRD::~AliAlgDetTRD()
{
  // d-tor
}

//____________________________________________
void AliAlgDetTRD::DefineVolumes()
{
  // define TRD volumes
  //
  const int kNSect = 18, kNStacks = 5, kNLayers = 6;
  AliAlgSensTRD *chamb=0;
  //
  //  AddVolume( volTRD = new AliAlgVol("TRD") ); // no main volume, why?
  AliAlgVol *sect[kNSect] = {0};
  //
  for (int ilr=0;ilr<kNLayers;ilr++) { // layer
    for (int ich=0;ich<kNStacks*kNSect;ich++) { // chamber
      Int_t isector   = ich/AliTRDgeometry::Nstack();
      Int_t istack    = ich%AliTRDgeometry::Nstack();
      Int_t lid       = AliTRDgeometry::GetDetector(ilr,istack,isector);
      const char *symname = Form("TRD/sm%02d/st%d/pl%d",isector,istack,ilr);
      if (!gGeoManager->GetAlignableEntry(symname)) continue;
      UShort_t vid    = AliGeomManager::LayerToVolUID(AliGeomManager::kTRD1+ilr,ich);
      AddVolume(chamb = new AliAlgSensTRD(symname,vid,lid) );
      if (!sect[isector]) sect[isector] = new AliAlgVol(Form("TRD/sm%02d",isector));
      chamb->SetParent(sect[isector]);
    } // chamber
  } // layer
  //
  for (int isc=0;isc<kNSect;isc++) {
    if (sect[isc]) AddVolume(sect[isc]);
  }
  //
}

//____________________________________________
Bool_t AliAlgDetTRD::PresentInTrack(const AliESDtrack* trc) const 
{
  // test if detector had seed this track
  return trc->IsOn(AliESDtrack::kTRDout);
}
