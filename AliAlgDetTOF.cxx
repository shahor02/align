#include "AliAlgDetTOF.h"
#include "AliAlgVol.h"
#include "AliAlgSensTOF.h"
#include "AliAlgSteer.h"
#include "AliGeomManager.h"
#include "AliTOFGeometry.h"
#include "AliESDtrack.h"
#include <TGeoManager.h>

ClassImp(AliAlgDetTOF);

//____________________________________________
AliAlgDetTOF::AliAlgDetTOF(const char* title)
{
  // default c-tor
  SetNameTitle(AliAlgSteer::GetDetNameByDetID(AliAlgSteer::kTOF),title);
  SetDetID(AliAlgSteer::kTOF);
}

//____________________________________________
AliAlgDetTOF::~AliAlgDetTOF()
{
  // d-tor
}

//____________________________________________
void AliAlgDetTOF::DefineVolumes()
{
  // define TOF volumes
  //
  const int kNSect = 18,kNStrips = AliTOFGeometry::NStripA()+2*AliTOFGeometry::NStripB()+2*AliTOFGeometry::NStripC();
  AliAlgSensTOF *strip=0;
  //
  //  AddVolume( volTOF = new AliAlgVol("TOF") ); // no main volume, why?
  AliAlgVol *sect[kNSect] = {0};
  //
  for (int isc=0;isc<kNSect;isc++) AddVolume(sect[isc]=new AliAlgVol(Form("TOF/sm%02d",isc)));
  
  int iid = -1;
  for (int isc=0;isc<kNSect;isc++) {
    for (int istr=1;istr<=kNStrips;istr++) { // strip
      int vid = AliGeomManager::LayerToVolUID(AliGeomManager::kTOF, ++iid);
      const char *symname = Form("TOF/sm%02d/strip%02d",isc,istr);
      if (!gGeoManager->GetAlignableEntry(symname)) continue;
      AddVolume(strip = new AliAlgSensTOF(symname,vid,iid,isc) );
      strip->SetParent(sect[isc]);
    } // strip
  } // layer
  //
}

//____________________________________________
Bool_t AliAlgDetTOF::AcceptTrack(const AliESDtrack* trc,Int_t trtype) const 
{
  // test if detector had seed this track
  return CheckFlags(trc,trtype);
}
