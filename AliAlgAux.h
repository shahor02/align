#ifndef ALIALGAUX_H
#define ALIALGAUX_H

#include <TMath.h>

using namespace TMath;

namespace AliAlgAux {
  const double kAlmostZeroD = 1e-15;
  const double kAlmostZeroF = 1e-11;
  const double kAlmostOneD = 1.-kAlmostZeroD;
  const double kAlmostOneF = 1.-kAlmostZeroF;
  const double kTinyDist   = 1.e-7; // ignore distances less that this
  //
  template<typename F> void   BringTo02Pi(F &phi);
  template<typename F> void   BringToPiPM(F &phi);
  void   PrintBits(ULong64_t patt, Int_t maxBits);
  Int_t  NumberOfBitsSet(UInt_t x);
  Bool_t SmallerAbs(double d, double tolD)    {return Abs(d)<tolD;}
  Bool_t SmallerAbs(float  f, double tolF)    {return Abs(f)<tolF;}
  Bool_t Smaller(double d, double tolD)       {return d<tolD;}
  Bool_t Smaller(float  f, double tolF)       {return f<tolF;}
  Bool_t IsZeroAbs(double d)                  {return SmallerAbs(d,kAlmostZeroD);}
  Bool_t IsZeroAbs(float  f)                  {return SmallerAbs(f,kAlmostZeroF);}
  Bool_t IsZeroPos(double d)                  {return Smaller(d,kAlmostZeroD);}
  Bool_t IsZeroPos(float  f)                  {return Smaller(f,kAlmostZeroF);}
}

//_________________________________________________________________________________
template<typename F>
inline void AliAlgAux::BringTo02Pi(F &phi) {
  // bring phi to 0-2pi range
  if (phi<0) phi+=TwoPi(); else if (phi>TwoPi()) phi-=TwoPi();
}

//_________________________________________________________________________________
template<typename F>
inline void AliAlgAux::BringToPiPM(F &phi) {
  // bring phi to -pi:pi range
  if (phi>TMath::Pi()) phi-=TwoPi();
}

//_________________________________________________________________________________
inline Int_t AliAlgAux::NumberOfBitsSet(UInt_t x) {
  // count number of non-0 bits in 32bit word
  x = x - ((x >> 1) & 0x55555555);
  x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
  return (((x + (x >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

#endif
