#ifndef ALIGLOALGMODULE_H
#define ALIGLOALGMODULE_H



class AliAlgSteerModule: public TNamed
{
  
 public:

  
 protected:

  TObjArray *fMatrices[];

  ClassDef(AliAlgSteerModule,1)  // base class for module global alignment
};


#endif
