#include "AliAlgVol.h"

ClassImp(AliAlgVol)

//------------------------------------------------------------
AliAlgVol::AliAlgVol(const char* name,const char* title) :
  TNamed(name,title)
  ,fFirstParOffs(-1)
  ,fParOffs(0)
  ,fDOF(0)
  ,fNDOF(0)
  ,fNDOFGeomFree(0)
  ,fNDOFFree(0)

  ,fParent(0)
  ,fChildren(0)

  ,fNProcPoints(0)
  ,fParVals(0)
  ,fParErrs(0)
  ,fParCstr(0)
{
  // def c-tor
}

//------------------------------------------------------------
AliAlgVol::~AliAlgVol()
{
  // d-tor
  delete[] fParOffs;
  delete[] fParVals;
  delete[] fParErrs;
  delete[] fParCstr;
  delete fChildren;
}
