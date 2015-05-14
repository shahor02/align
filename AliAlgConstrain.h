#ifndef ALIALGCONSTRAINT_H
#define ALIALGCONSTRAINT_H

#include <TNamed.h>
#include <TObjArray.h>
#include "AliAlgVol.h"

/*--------------------------------------------------------
  Base class for detector: wrapper for set of volumes
  -------------------------------------------------------*/


class AliAlgConstraint : public TNamed
{
 public:
  enum {kNDOFGeom=AliAlgVol::kNDOFGeom};

  AliAlgConstraint(const char* name=0,const char* title=0);
  virtual ~AliAlgConstraint();
  //
  void        SetParent(const AliAlgVol* par)        {fParent = par;}
  AliAlgVol*  GetParent()                      const {return fParent;}
  //
  Int_t       GetNChildren()                   const {return fChildren.GetEntriesFast();}
  AliAlgVol*  GetChild(int i)                  const {return (AliAlgVol*)fChildren[i];}
  void        AddChild(const AliAlgVol* v)           {fChildren.AddLast(v);}
  //
  Bool_t     IsDOFConstrained(Int_t dof)       const {return fConstraint&0x1<<dof;}
  UChar_t    GetConstraintPattern()            const {return fConstraint;}
  void       ConstrainDOF(Int_t dof)                 {fConstraint |= 0x1<<dof;}
  void       UContrainDOF(Int_t dof)                 {fConstraint &=~(0x1<<dof);}
  void       SetConstrainPattern(UInt_t pat)         {fConstraint = pat;}
  Bool_t     HasConstraint()                   const {return  fConstraint;}
  //
  virtual void   WriteChildrenConstraints(FILE* conOut) const;
  virtual void   CheckConstraints()                     const;
  virtual const char* GetDOFName(int i)                 const {AliAlgVol::GetGeomDOFName(i);}
  //
 protected:
  UInt_t      fConstraint;          // bit pattern of constraint
  AliAlgVol*  fParent;              // parent volume for contraint, lab if 0
  TObjArray   fChildren;            // volumes subjected to constraints
  //
  ClassDef(AliAlgConstraint,1);
};

#endif
