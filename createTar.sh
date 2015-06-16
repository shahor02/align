#!/bin/bash

tar cvzf algGlo.tar.gz \
AliAlgAux.{h,cxx}  AliAlgDetITS.{h,cxx} AliAlgDetTPC.{h,cxx}  AliAlgMPRecord.{h,cxx}  \
AliAlgRes.{h,cxx}  AliAlgSensITS.{h,cxx} AliAlgSensTPC.{h,cxx}  AliAlgSteer.{h,cxx}  \
AliAlgVol.{h,cxx}  AliAlgDet.{h,cxx}  AliAlgDetTOF.{h,cxx}  AliAlgDetTRD.{h,cxx}  \
AliAlgPoint.{h,cxx}  AliAlgSens.{h,cxx} AliAlgSensTOF.{h,cxx}  AliAlgSensTRD.{h,cxx} \
AliAlgTrack.{h,cxx}  AliAlgVtx.{h,cxx} AliAlgConstraint.{h,cxx} AliAlgResFast.{h,cxx} \
AliAlgDOFStat.{h,cxx} AlgLinkDef.h \
Mille.{h,cxx} \
Makefile

