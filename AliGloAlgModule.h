#ifndef ALIGLOALGMODULE_H
#define ALIGLOALGMODULE_H



class AliGloAlgModule: public TNamed
{
  
 public:

  
 protected:

  TObjArray *fMatrices[];

  ClassDef(AliGloAlgModule,1)  // base class for module global alignment
};


#endif
