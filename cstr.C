#include "TGeoMatrix.h"

void PrintMat(const char* name, const double mt[4][4])
{
  printf("\nD/D%s\n",name);
  for (int i=0;i<4;i++) {
    for (int j=0;j<4;j++) printf("%+e ",mt[i][j]); printf("\n");
  }
}

void cstr(TGeoHMatrix &matRD, double jac[6][6])
{
  // If the transformation R brings the vector from "local" frame to "master" frame as V=R*v
  // then application of the small LOCAL correction tau to vector v is equivalent to
  // aplication of correction TAU in MASTER framce V' = R*tau*v = TAU*R*v
  // with TAU = R*tau*R^-1
  // Constraining the LOCAL modifications of child volumes to have 0 total movement in their parent
  // frame is equivalent to request that sum of all TAU matrices is unity matrix, or TAU-I = 0.
  //
  // This routine calculates derivatives of the TAU-I matrix over local corrections x,y,z, psi,tht,phi
  // defining matrix TAU. In small corrections approximation the constraint is equivalent to
  // Sum_over_child_volumes{ [dTAU/dParam]_ij * deltaParam } = 0
  // for all elements ij of derivative matrices. Since only 6 out of 16 matrix params are independent,
  // we request the constraint only for  [01],[02],[13],[30],[31],[32]
  //
  TGeoHMatrix matRI = matRD.Inverse();
  const int ij[6][2] = {{0,1},{0,2},{1,2},{3,0},{3,1},{3,2}}; 
  //
  const double *rd=matRD.GetRotationMatrix(),*ri=matRI.GetRotationMatrix();  
  const double *td=matRD.GetTranslation(),   *ti=matRI.GetTranslation();
  //
  double dDPar[6][4][4] = {
    // dDX[4][4] 
    {{0,0,0,0},{0,0,0,0},{0,0,0,0},{rd[0],rd[3],rd[6],0}},
    // dDY[4][4]
    {{0,0,0,0},{0,0,0,0},{0,0,0,0},{rd[1],rd[4],rd[7],0}},
    // dDZ[4][4]
    {{0,0,0,0},{0,0,0,0},{0,0,0,0},{rd[2],rd[5],rd[8],0}},
    // dDPSI[4][4]
    {{rd[2]*ri[3]-rd[1]*ri[6],rd[2]*ri[4]-rd[1]*ri[7],rd[2]*ri[5]-rd[1]*ri[8],0},
     {rd[5]*ri[3]-rd[4]*ri[6],rd[5]*ri[4]-rd[4]*ri[7],rd[5]*ri[5]-rd[4]*ri[8],0},
     {rd[8]*ri[3]-rd[7]*ri[6],rd[8]*ri[4]-rd[7]*ri[7],rd[8]*ri[5]-rd[7]*ri[8],0},
     {rd[2]*ti[1]-rd[1]*ti[2],rd[5]*ti[1]-rd[4]*ti[2],rd[8]*ti[1]-rd[7]*ti[2],0}},
    // dDTHT[4][4]
    {{rd[0]*ri[6]-rd[2]*ri[0], rd[0]*ri[7]-rd[2]*ri[1], rd[0]*ri[8]-rd[2]*ri[2],0},
     {rd[3]*ri[6]-rd[5]*ri[0], rd[3]*ri[7]-rd[5]*ri[1], rd[3]*ri[8]-rd[5]*ri[2],0},
     {rd[6]*ri[6]-rd[8]*ri[0], rd[6]*ri[7]-rd[8]*ri[1], rd[6]*ri[8]-rd[8]*ri[2],0},
     {rd[0]*ti[2]-rd[2]*ti[0], rd[3]*ti[2]-rd[5]*ti[0], rd[6]*ti[2]-rd[8]*ti[0],0}},
    // dDPHI[4][4]
    {{rd[1]*ri[0]-rd[0]*ri[3],rd[1]*ri[1]-rd[0]*ri[4],rd[1]*ri[2]-rd[0]*ri[5],0},
     {rd[4]*ri[0]-rd[3]*ri[3],rd[4]*ri[1]-rd[3]*ri[4],rd[4]*ri[2]-rd[3]*ri[5],0},
     {rd[7]*ri[0]-rd[6]*ri[3],rd[7]*ri[1]-rd[6]*ri[4],rd[7]*ri[2]-rd[6]*ri[5],0},
     {rd[1]*ti[0]-rd[0]*ti[1],rd[4]*ti[0]-rd[3]*ti[1],rd[7]*ti[0]-rd[6]*ti[1],0}},
  };
  //
  for (int cs=0;cs<6;cs++) {
    int i=ij[cs][0],j=ij[cs][1];
    for (int ip=0;ip<6;ip++) jac[cs][ip] = dDPar[ip][i][j];
  }
}
