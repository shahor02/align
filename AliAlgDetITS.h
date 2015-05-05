#ifndef ALIALGDETITS_H
#define ALIALGDETITS_H

#include "AliAlgDet.h"

class AliAlgDetITS : public AliAlgDet
{
 public:
  AliAlgDetITS(const char* title="");
  virtual ~AliAlgDetITS();
  //
  virtual void DefineVolumes();  
  //
  Bool_t AcceptTrack(const AliESDtrack* trc) const;

  void   SetAddErrorLr(int ilr, double sigY, double sigZ);
  void   SetSkipLr(int ilr);
  //
  virtual void  UpdatePointByTrackInfo(AliAlgPoint* pnt, const AliExternalTrackParam* t) const;
  virtual void  SetUseErrorParam(Int_t v);
  //
 protected:

  void GetErrorParamAngle(int layer,double tgl,double tgphitr,double &erry,double &errz) const;

 protected:
  //
  ClassDef(AliAlgDetITS,1);
};

#endif
