#include "AliAlgAux.h"
#include "AliCDBId.h"
#include "AliCDBManager.h"
#include "AliLog.h"
#include <TList.h>
#include <TMap.h>
#include <TObjString.h>
#include <TPRegexp.h>
#include <TGrid.h>

//_______________________________________________________________
void AliAlgAux::PrintBits(ULong64_t patt, Int_t maxBits)
{
  // print maxBits of the pattern
  maxBits = Min(64,maxBits);
  for (int i=0;i<maxBits;i++) printf("%c",((patt>>i)&0x1) ? '+':'-');
}

//_______________________________________________________________
AliCDBId* AliAlgAux::FindCDBId(const TList* cdbList,const TString& key)
{
  // Find enty for the key in the cdbList and create its CDBId
  // User must take care of deleting created CDBId
  TIter next(cdbList);
  TObjString* entry;
  while ( (entry=(TObjString*)next()) ) if (entry->GetString().Contains(key)) break;    
  if (!entry) return 0;
  return AliCDBId::MakeFromString(entry->GetString());
  //
}

//_________________________________________________________
void AliAlgAux::RectifyOCDBUri(TString& inp)
{
  // render URI from cdbMap to usable form
  TString uri = "";
  int ind;
  if (inp.BeginsWith("alien:/")) { // alien folder
    TPRegexp fr("[Ff]older=/");
    if ( (ind=inp.Index(fr))>0 ) inp.Remove(0,ind);
    ind = inp.First('?');
    if (ind>0) inp.Resize(ind);
    inp.Prepend("alien://");
  }
  else if (inp.BeginsWith("local:/")) {
    ind = inp.First('?');
    if (ind>0) inp.Resize(ind);
  }
  else {
    AliFatalGeneralF("::RectifyOCDBUri","Failed to extract OCDB URI from %s",inp.Data());
  }
  //
}

//_________________________________________________________
Bool_t AliAlgAux::PreloadOCDB(int run, const TMap* cdbMap, const TList* cdbList)
{
  // Load OCDB paths for given run from pair of cdbMap / cdbList
  // as they are usually stored in the UserInfo list of esdTree
  // In order to avoid unnecessary uploads, the objects are not actually 
  // loaded/cached but just added as specific paths with version
  //  
  TObjString *ostr,*okey;
  TString uriDef,uri,key;
  //
  CleanOCDB();
  //
  ostr = (TObjString*)cdbMap->GetValue("default");
  RectifyOCDBUri( uriDef=ostr->GetString() );
  AliInfoGeneralF("","Default storage %s",uriDef.Data());
  //
  AliCDBManager* man = AliCDBManager::Instance();  
  man->SetDefaultStorage(uriDef.Data());
  man->SetRun(run);
  //
  TIter nextM(cdbMap);
  while ( (okey=(TObjString*)nextM()) ) {
    if ( (key=okey->GetString())=="default") continue;
    ostr = (TObjString*)cdbMap->GetValue(okey);
    RectifyOCDBUri( uri=ostr->GetString() );
    // fetch object from the list 
    AliCDBId* cdbID = FindCDBId(cdbList,key);
    int ver=-1,sver=-1;
    if (cdbID) {
      ver = cdbID->GetVersion();
      sver= cdbID->GetSubVersion();
      delete cdbID;
    }
    else {
      AliWarningGeneralF("::PreloadOCDB","Key %s has special storage %s but absent in the cdbList",
			 key.Data(),uri.Data());
    }
    AliInfoGeneralF("::PreloadOCDB","Setting storage for %s to %s",key.Data(),uri.Data());
    man->SetSpecificStorage(key.Data(),uri.Data(),ver,sver);
  }
  //
  // now set remaining objects from the list as specific storages
  TIter nextL(cdbList);
  while ( (ostr=(TObjString*)nextL()) ) {
    AliCDBId *cdbIDF = AliCDBId::MakeFromString(ostr->GetString());
    cdbIDF->Print();
    man->SetSpecificStorage(cdbIDF->GetPath().Data(),uriDef.Data(),
			    cdbIDF->GetVersion(),cdbIDF->GetSubVersion());
    delete cdbIDF;
  }
  //
  return kTRUE;
}

//_________________________________________________________
void AliAlgAux::CleanOCDB()
{
  // brings OCDB to virgin state
  Bool_t isGrid = gGrid!=0;
  AliCDBManager::Destroy();
  if (isGrid && !gGrid) TGrid::Connect("alien://");
}
