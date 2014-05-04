#ifndef ALIALGAUX_H
#define ALIALGAUX_H

#include <TMath.h>

using namespace TMath;

namespace AliAlgAux {
  const double kAlmostZeroD = 1e-15;
  const double kAlmostZeroF = 1e-11;
  const double kAlmostOneD = 1.-kAlmostZeroD;
  const double kAlmostOneF = 1.-kAlmostZeroF;
  //
  void   PrintBits(ULong64_t patt, Int_t maxBits);
  Bool_t SmallerAbs(double d, double tolD)    {return Abs(d)<tolD;}
  Bool_t SmallerAbs(float  f, double tolF)    {return Abs(f)<tolF;}
  Bool_t Smaller(double d, double tolD)       {return d<tolD;}
  Bool_t Smaller(float  f, double tolF)       {return f<tolF;}
  Bool_t IsZeroAbs(double d)                  {return SmallerAbs(d,kAlmostZeroD);}
  Bool_t IsZeroAbs(float  f)                  {return SmallerAbs(f,kAlmostZeroF);}
  Bool_t IsZeroPos(double d)                  {return Smaller(d,kAlmostZeroD);}
  Bool_t IsZeroPos(float  f)                  {return Smaller(f,kAlmostZeroF);}
}

#endif
