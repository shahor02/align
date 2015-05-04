#ifndef ALIALGVTX_H
#define ALIALGVTX_H

/*--------------------------------------------------------
  Special fake "sensor" for event vertex.
  It is needed to allow adjustement of the global IP position
  if the event event is used as a measured point.
  Its degrees of freedom of LOCAL X,Y,Z, coinciding with
  GLOBAL X,Y,Z. 
  Since the vertex added to the track as a mesured point must be
  defined in the frame with X axis along the tracks, the T2L
  matrix of this sensor need to be recalculated for each track!
  -------------------------------------------------------*/
#include "AliAlgSens.h"

class AliAlgVtx : public AliAlgSens
{
 public:
  AliAlgVtx();
  //
  void SetAlpha(double alp)              {fAlp=alp; PrepareMatrixT2L();}
  virtual void   PrepareMatrixL2G()      {fMatL2G.Clear();} // unit matrix
  virtual void   PrepareMatrixL2GIdeal()  {fMatL2GIdeal.Clear();} // unit matrix
  virtual void   PrepareMatrixT2L();
  
 protected:
  //
  ClassDef(AliAlgVtx,1);
};


#endif
