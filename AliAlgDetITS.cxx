#include "AliAlgDetITS.h"
#include "AliAlgVol.h"
#include "AliAlgSens.h"
#include "AliITSgeomTGeo.h"
#include "AliGeomManager.h"

ClassImp(AliAlgDetITS);

//____________________________________________
AliAlgDetITS::AliAlgDetITS()
{
  // default c-tor
}

//____________________________________________
AliAlgDetITS::AliAlgDetITS(const char* name, const char* title)
  :AliAlgDet(name,title)
{

}

//____________________________________________
AliAlgDetITS::~AliAlgDetITS()
{
  // d-tor
}

//____________________________________________
void AliAlgDetITS::DefineVolumes()
{
  // define ITS volumes
  //
  const int kNSPDSect = 10;
  AliAlgVol *volITS=0,*hstave=0,*ladd=0;
  AliAlgSens *sens=0;
  //
  AddVolume( volITS = new AliAlgVol("ITS") );
  int cntVID=0;
  //
  // SPD
  AliAlgVol *sect[kNSPDSect] = {0};
  for (int isc=0;isc<kNSPDSect;isc++) { // sectors
    AddVolume( sect[isc] = new AliAlgVol(Form("ITS/SPD0/Sector%d",isc)) );
    sect[isc]->SetParent(volITS);
  }
  for (int ilr=0;ilr<=1;ilr++) { // SPD layers
    //
    cntVID = 0;
    int nst = AliITSgeomTGeo::GetNLadders(ilr+1)/kNSPDSect; // 2 or 4 staves per sector
    for (int isc=0;isc<kNSPDSect;isc++) { // sectors
      for (int ist=0;ist<nst;ist++) { // staves of SPDi
	for (int ihst=0;ihst<2;ihst++) { // halfstave
	  AddVolume ( hstave = new AliAlgVol(Form("ITS/SPD%d/Sector%d/Stave%d/HalfStave%d",ilr,isc,ist,ihst)) );
	  hstave->SetParent(sect[isc]);
	  for (int isn=0;isn<2;isn++) { // "ladder" (sensor)	    
	    AddVolume( sens = new AliAlgSens(Form("ITS/SPD%d/Sector%d/Stave%d/HalfStave%d/Ladder%d",ilr,isc,ist,ihst,isn), 
					     AliGeomManager::LayerToVolUID(ilr+1,cntVID++)) );
	    sens->SetParent(hstave);
	  }
	}
      } // staves of SPDi
    } // sectors
  } // SPD layers
  //
  // SDD
  for (int ilr=2;ilr<=3;ilr++) { // layer
    cntVID = 0;
    for (int ist=0;ist<AliITSgeomTGeo::GetNLadders(ilr+1);ist++) { // ladder
      AddVolume( ladd = new AliAlgVol(Form("ITS/SDD%d/Ladder%d",ilr,ist)) );
      ladd->SetParent(volITS);
      for (int isn=0;isn<AliITSgeomTGeo::GetNDetectors(ilr+1);isn++) { // sensor
	AddVolume( sens = new AliAlgSens(Form("ITS/SDD%d/Ladder%d/Sensor%d",ilr,ist,isn), 
					 AliGeomManager::LayerToVolUID(ilr+1,cntVID++)) );
	sens->SetParent(ladd); 
      }
    } // ladder
  } // layer
  //
  // SSD
  for (int ilr=4;ilr<=5;ilr++) { // layer
    cntVID = 0;
    for (int ist=0;ist<AliITSgeomTGeo::GetNLadders(ilr+1);ist++) { // ladder
      AddVolume( ladd = new AliAlgVol(Form("ITS/SSD%d/Ladder%d",ilr,ist)) );
      ladd->SetParent(volITS);
      for (int isn=0;isn<AliITSgeomTGeo::GetNDetectors(ilr+1);isn++) { // sensor
	AddVolume( sens = new AliAlgSens(Form("ITS/SDD%d/Ladder%d/Sensor%d",ilr,ist,isn),
					 AliGeomManager::LayerToVolUID(ilr+1,cntVID++)) );
	sens->SetParent(ladd); 
      }
    } // ladder
  } // layer
  //
  //
}

//____________________________________________
void AliAlgDetITS::PrintHierarchy()
{
  // print ITS volumes
  for (int iv=0;iv<GetNVolumes();iv++) {
    const AliAlgVol* vol = GetVolume(iv);
    int offs = vol->CountParents();
    for (int i=offs;i--;) printf("  ");
    if (vol->IsSensor()) printf(" VId: %6d ", ((AliAlgSens*)vol)->GetVolID());
    printf("%s\n",vol->GetName());
  }
  //
}
