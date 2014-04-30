void load()
{
  gROOT->ProcessLine(".L AliAlgAux.cxx+");
  gROOT->ProcessLine(".L AliAlgPoint.cxx+");
  gROOT->ProcessLine(".L AliAlgTrack.cxx+");
}
