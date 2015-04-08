#include "AliAlgDetTPC.h"
#include "AliAlgVol.h"
#include "AliAlgSensTPC.h"
#include "AliAlgSteer.h"
#include "AliGeomManager.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include <TGeoManager.h>

ClassImp(AliAlgDetTPC);

//____________________________________________
AliAlgDetTPC::AliAlgDetTPC()
{
  // default c-tor
  SetDetID(AliAlgSteer::kTPC);
}

//____________________________________________
AliAlgDetTPC::AliAlgDetTPC(const char* name, const char* title)
  :AliAlgDet(name,title)
{
  SetDetID(AliAlgSteer::kTPC);
}

//____________________________________________
AliAlgDetTPC::~AliAlgDetTPC()
{
  // d-tor
}

//____________________________________________
void AliAlgDetTPC::DefineVolumes()
{
  // define TPC volumes
  //
  const int kNSect = 18, kAC=2, kIOROC=2;
  const char* kSide[kAC] = {"A","C"};
  const char* kROC[kIOROC] = {"Inner","Outer"};
  //  AliAlgSensTPC *chamb=0;
  //
  AliAlgVol* volTPC = new AliAlgVol("ALIC_1/TPC_M_1");
  AddVolume( volTPC ); 
  //
  for (int roc=0;roc<kIOROC;roc++) { // inner/outer
    for (int side=0;side<kAC;side++) { // A/C
      for (int isc=0;isc<kNSect;isc++) { // sector ID
	const char *symname = Form("TPC/Endcap%s/Sector%d/%sChamber",kSide[side],isc+1,kROC[roc]);
	if (!gGeoManager->GetAlignableEntry(symname)) {
	  AliErrorF("Did not find alignable %s",symname);
	  continue;
	}
	Int_t iid = side*kNSect+isc;
	UShort_t vid = AliGeomManager::LayerToVolUID(AliGeomManager::kTPC1+roc,iid);
	AliAlgSensTPC* sens = new AliAlgSensTPC(symname,vid,iid,isc);
	sens->SetParent(volTPC);
	AddVolume(sens);
      } // sector ID
    } // A/C
  } // inner/outer
  //
}

//____________________________________________
Bool_t AliAlgDetTPC::PresentInTrack(const AliESDtrack* trc) const 
{
  // test if detector had seed this track
  return trc->IsOn(AliESDtrack::kTPCrefit);
}
