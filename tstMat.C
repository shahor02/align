void tstMat(AliAlgSens* sens, double locOr[3])
{
  //
  sens->SetMatrixL2GReco(sens->GetMatrixL2G());

  //
  printf("LocOR  : ");for (int i=0;i<3;i++) printf("%+9.4f ",locOr[i]); printf("\n");
  double tmp[3];
  sens->GetMatrixT2L().MasterToLocal(locOr,tmp);
  double hit[3]={0};
  hit[0]=tmp[1];
  hit[1]=tmp[2];
  double fXYZ[3] = {0,hit[0],hit[1]};
  Double_t fXYZLoc[3] = {0, 0, 0};
  sens->GetMatrixT2L().LocalToMaster(fXYZ,fXYZLoc);
  fXYZLoc[1] = 0;

  printf("TraCLId: ");for (int i=0;i<3;i++) printf("%+9.4f ",fXYZ[i]); printf("\n");
  printf("LocCLId: ");for (int i=0;i<3;i++) printf("%+9.4f ",fXYZLoc[i]); printf("\n");

  double traOrId[3];
  sens->GetMatrixT2L().MasterToLocal(locOr,traOrId);
  traOrId[0] = 0;
  printf("TraORId: ");for (int i=0;i<3;i++) printf("%+9.4f ",traOrId[i]); printf("\n");
  //
  double traOrMl[3];
  sens->GetMatrixClAlg().LocalToMaster(traOrId,traOrMl);
  printf("TraORMl: ");for (int i=0;i<3;i++) printf("%+9.4f ",traOrMl[i]); printf("\n");
  //
  Double_t lxyz[3] = {0, 0, 0};
  double gloOrMlW[3],gloOrMl[3];
  sens->GetMatrixT2L().LocalToMaster(traOrMl,lxyz);
  printf("LocORMl: ");for (int i=0;i<3;i++) printf("%+9.4f ",lxyz[i]); printf("\n");
  sens->GetMatrixL2GReco().LocalToMaster(lxyz,gloOrMlW);
  sens->GetMatrixL2GIdeal().LocalToMaster(lxyz,gloOrMl);
  //
  printf("Wrong transform to global:\n");  
  printf("GloORMl: ");for (int i=0;i<3;i++) printf("%+9.4f ",gloOrMlW[i]); printf("\n");
  printf("Correct transform to global:\n");  
  printf("GloORMl: ");for (int i=0;i<3;i++) printf("%+9.4f ",gloOrMl[i]); printf("\n");

  //

  // backward
  printf("\nREVERT\n");
  //
  const TGeoHMatrix& matL2Grec = sens->GetMatrixL2GReco(); // local to global matrix used for reconstruction
  //const TGeoHMatrix& matL2G    = sens->GetMatrixL2G();     // local to global orig matrix used as a reference 
  const TGeoHMatrix& matT2L    = sens->GetMatrixT2L();     // matrix for tracking to local frame translation
  
  double loc[3],locId[3];
  matL2Grec.MasterToLocal(gloOrMl,locId); // go to local frame using reco-time matrix 
  printf("LocId  : ");for (int i=0;i<3;i++) printf("%+9.4f ",   locId[i]); printf("\n");  
  sens->GetMatrixL2GIdeal().MasterToLocal(gloOrMl,loc);
  printf("LocMl  : ");for (int i=0;i<3;i++) printf("%+9.4f ",   loc[i]); printf("\n");  

  double traId[3];
  matT2L.MasterToLocal(locId,traId); // go to tracking frame 
  printf("TraId  : ");for (int i=0;i<3;i++) printf("%+9.4f ", traId[i]); printf("\n");    
  double tra[3];
  sens->GetMatrixClAlg().LocalToMaster(traId,tra);   // apply alignment
  printf("Tra    : ");for (int i=0;i<3;i++) printf("%+9.4f ", tra[i]); printf("\n");    
  //
  double glo[3];
  double fAlphaSens = sens->GetAlpTracking();
  double xtra = sens->GetXTracking();
  //
  double cs=TMath::Cos(fAlphaSens);
  double sn=TMath::Sin(fAlphaSens);
  double x= xtra + tra[0];
  double r[3];
  r[0] = x*cs - tra[1]*sn; 
  r[1] = x*sn + tra[1]*cs;
  r[2] = tra[2];
  printf("Glo1   : ");for (int i=0;i<3;i++) printf("%+9.4f ", r[i]); printf("\n");   
  //
  TGeoHMatrix t2g;
  sens->GetMatrixT2G(t2g);
  t2g.LocalToMaster(tra,glo);
  printf("GloT2G : ");for (int i=0;i<3;i++) printf("%+9.4f ", glo[i]); printf("\n");
  //
}
