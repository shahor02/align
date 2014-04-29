
void Load();

AliAlgTrack* algTrack=0;

const int kNITS = 6;
double rITS[kNITS] = {3.9,7.6,15.0,23.9,38.0,43.0};

Bool_t TestTrack(const AliExternalTrackParam& trSrc)
{
  Load();
  AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  Double_t bz = fld->SolenoidField();
  //
  AliExternalTrackParam tr0(trSrc);
  //
  algTrack = new AliAlgTrack();
  //
  algTrack.AliExternalTrackParam::operator=(tr0);
  //
  // add points
  for (int i=0;i<kNITS;i++) {
    if (!tr0->PropagateTo(rITS,bz)) return kFALSE;
    //
    AliAlgPoint* pnt = new AliAlgPoint();
    //
    
  }

}


//________________________________________________________________________________
void Load()
{
  if (!AliGeomManager::GetGeometry()) AliGeomManager::LoadGeometry("geomery.root");
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    AliMagF* fld = new AliMagF("bmap","bmap");
    TGeoGlobalMagField::Instance()->SetField( fld );
    TGeoGlobalMagField::Instance()->Lock();
  }
}
