/**************************************************************

    Name:    ProtonPathComputation.cpp

    Purpose: This program is used to generate phantoms for
             simulation purposes in pCT. It's strengths lie
             in its ability to generate custom phantoms of
             varying complexity.

    Credits: Micah ______, Howard Heaton
             
    Micah originally wrote this code during his work for
    his Master's thesis. Then extensive commenting and
    appropriate modifications have been made by Howard
    Heaton to incorporate blob basis functions in the
    system matrix.

    Revision Date: July 16, 2014

**************************************************************/

#ifndef   __IncludedLibraries__
#define   __IncludedLibraries__

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <string.h>
#include <time.h>
#include <fstream>
#include <sstream>

#endif // __IncludedLibraries__

/************************************************************************/

/***********************************PROTONPATH*******************************/

class PointVec {
 public:
  PointVec(float xx,float yy) : x(xx), y(yy) {};
  PointVec() {};
  void set(float xx,float yy) {
    x=xx;
    y=yy;
  }
  float x;
  float y;
};

class PointVec3 {
public:
  PointVec3(float xx,float yy,float zz) : x(xx),y(yy),z(zz) {};
  float x;
  float y;
  float z;
};

//Returns the Euclidean Norm of Vector connecting
//P1 to P2
float normEu(PointVec P1,PointVec P2) {
  return sqrt(pow(P2.x-P1.x,2)+pow(P2.y-P1.y,2));
}

// each PointVec3  (varLat, cov, VarAng)
vector<PointVec3> scatteringParameters() {
  vector<PointVec3> output;
  output.push_back(PointVec3(0,0,0));
  output.push_back(PointVec3(0.00112,  0.0001686, 3.397*pow(10,-5)));
  output.push_back(PointVec3(0.009335, 0.0007052, 7.154*pow(10,-5)));
  output.push_back(PointVec3(0.0324,   0.001638,  0.0001117));
  output.push_back(PointVec3(0.07861,  0.002992,  0.0001542));
  output.push_back(PointVec3(0.1567,   0.004793,  0.0001994));
  output.push_back(PointVec3(0.2761,   0.007067,  0.0002472));
  output.push_back(PointVec3(0.4466,   0.009843,  0.0002979));
  output.push_back(PointVec3(0.6786,   0.01315,   0.0003519));
  output.push_back(PointVec3(0.9833,   0.01703,   0.0004094));
  output.push_back(PointVec3(1.372,    0.0215,    0.0004709));
  output.push_back(PointVec3(1.859,    0.02663,   0.0005368));
  output.push_back(PointVec3(2.456,    0.03245,   0.0006078));
  output.push_back(PointVec3(3.178,    0.03902,   0.0006847));
  output.push_back(PointVec3(4.041,    0.0464,    0.0007683));
  output.push_back(PointVec3(5.063,    0.05467,   0.0008599));
  output.push_back(PointVec3(6.261,    0.06392,   0.0009611));
  output.push_back(PointVec3(7.658,    0.07425,   0.001074));
  output.push_back(PointVec3(9.275,    0.08579,   0.001201));
  output.push_back(PointVec3(11.14,    0.09871,   0.001347));
  output.push_back(PointVec3(13.28,    0.1132,    0.001518));
  return output;
}

// return uniform random number from [0,1]
float ranf() {
  return rand()%1000001/(float)1000000;
}

// return uniform rv from [-125,125]
float MyUniformRand() {
  float rv=rand()%1000001/(float)1000000;
  rv=rv*250.0; //scale rv to be in [0,250]
  rv=rv-125.0; //shift rv to be in [-125,125]
  return rv;
}

//for vertical displacement
float MyUniformRand2() {
  float rv=rand()%1000001/(float)1000000;
  rv=rv*100.0; //scale rv to be in [0,250]
  rv=rv-50.0; //shift rv to be in [-125,125]
  return rv;
}

// Marsaglia Polar Method ~ modified from Box-Muller
// to take 2 uniform rvs and transform them into
// 2 independent standard normal rvs
PointVec generate2RandStdNorm() {
  float u,v,x,y,s,t;
  s=1.0;
  while (s>=1.0 or equalf(s,0)) {
    u=2.0*ranf()-1.0;
    v=2.0*ranf()-1.0;
    s=u*u+v*v;
  }
  t=sqrt((-2.0*log(s))/s);
  x=u*t;
  y=v*t;
  return PointVec(x,y);
}

// hard coded calculation of determinant
float Determinant3x3(float* A) {
  float d=A[0]*(A[4]*A[8]+A[5]*A[7])
         -A[1]*(A[3]*A[8]+A[5]*A[6])
         +A[2]*(A[3]*A[7]+A[4]*A[6]);
  return d;
}

// solve linear 3x3 system of equations Ax=b
float* CramerSolve3x3(float* A,float* b) {
  float* output=(float*)malloc(3*sizeof(float));
  float detA=Determinant3x3(A);
  if (detA==0) {
    cout<<"error: zero determinant for J"<<endl;
    output[0]=0;
    output[1]=0;
    output[2]=0;
    return output;
  }
  float* tempM=(float*)malloc(9*sizeof(float));
  for (int i=0; i<9; i++)
    tempM[i]=A[i];
  //calculate the x of x=(x,y,z)
  tempM[0]=b[0];
  tempM[3]=b[1];
  tempM[6]=b[2];
  float x=Determinant3x3(tempM)/detA;
  tempM[0]=A[0];
  tempM[3]=A[3];
  tempM[6]=A[6];
  //calculate y
  tempM[1]=b[0];
  tempM[4]=b[1];
  tempM[7]=b[2];
  float y=Determinant3x3(tempM)/detA;
  tempM[1]=A[1];
  tempM[4]=A[4];
  tempM[7]=A[7];
  //calculate z
  tempM[2]=b[0];
  tempM[5]=b[1];
  tempM[8]=b[2];
  float z=Determinant3x3(tempM)/detA;
  output[0]=x;
  output[1]=y;
  output[2]=z;
  return output;
}

// F is vector function for nonlinear system
// of equations, this function returns -F
// which is needed in newton method
// v = (var1,cov,var2)
float* MyF(float*x, float* v) {
  float* F=(float*)malloc(3*sizeof(float));
  F[0]=-(x[0]*x[0]+x[1]*x[1]-v[0]);
  F[1]=-(x[1]*x[1]+x[2]*x[2]-v[2]);
  F[2]=-(x[0]*x[1]+x[1]*x[2]-v[1]);
  return F;
}

// hard coded calculation of the jacobian
// for the nonlinear system of equation to
// solve for the entries in matrix to
// multiply to the vector of std norms to
// get the joint norms
float* MyJacobian(float* x) {
  float* J=(float*)malloc(9*sizeof(float));
  J[0]=2*x[0];
  J[1]=2*x[1];
  J[2]=0;
  J[3]=0;
  J[4]=2*x[1];
  J[5]=2*x[2];
  J[6]=x[1];
  J[7]=x[0]+x[2];
  J[8]=x[1];
  return J;
}

// J~jacobian of F, v is needed to calculate F
// x~solution vector
float* NewtonNonlinearSysMethod(float* v) {
  float *J,*F,*s;
  float* x=(float*)malloc(3*sizeof(float));
  float* xprevious=(float*)malloc(3*sizeof(float));
  x[0]=10; x[1]=0.1; x[2]=0.01; //intitial x guess
  for (int k=0; k<20; k++) { //20 iterations
    J=MyJacobian(x); // 3x3 jacobian
    F=MyF(x,v); //returns 3x1 vector
    //solve Js=-f for s using Cramer
    s=CramerSolve3x3(J,F);
    x[0]=x[0]+s[0];
    x[1]=x[1]+s[1];
    x[2]=x[2]+s[2];
    free(J); free(F); free(s);
  }
  return x;
}

PointVec generateJointRandNorm(float var1,float cov,float var2) {
  PointVec U=generate2RandStdNorm();;
  float* v=(float*)malloc(3*sizeof(float));
  v[0]=var1; v[1]=cov; v[2]=var2;
  float* temp=NewtonNonlinearSysMethod(v);
  PointVec X;
  X.x=temp[0]*U.x+temp[1]*U.y;
  X.y=temp[1]*U.x+temp[2]*U.y;
  free(v); free(temp);
  return X;
}

// m~slopeLine,b~yinterceptLine,xa~lowerXboundBox,
// ya~lowerYboundBox,xb~UpperXboundBox,yb~upperYboundBox
// function returns the 2 points of interesection
// of the line and the box
vector<PointVec> LineBoxHits(float m,float b,float xa,
 float xb,float ya,float yb) {
  // calculate intersection of proton entry line with
  // reconstruction box
  float xdot,ydot;
  vector<PointVec> v; //entry and exit of the box
  xdot=(1.0/m)*(yb-b);
  if (xa<=xdot and xdot<=xb) //hit top line
    v.push_back(PointVec(xdot,yb));
  xdot=(1.0/m)*(ya-b);
  if (xa<=xdot and xdot<=xb) //hit bottom line
    v.push_back(PointVec(xdot,ya));
  ydot=m*xb+b;
  if (ya<=ydot and ydot<=yb) //hit right line
    v.push_back(PointVec(xb,ydot));
  ydot=m*xa+b;
  if (ya<=ydot and ydot<=yb) //hit left line
    v.push_back(PointVec(xa,ydot));
  if (v.size()!=2)
    cout<<"number of interesections="<<v.size()<<endl;
  return v;
}
// Point (x0,y0) is for the upper left corner of each pixel
float CalculateX0(float x,float xstep) {
  float tempDivide=x/xstep;
  float remainPercent=tempDivide-floor(tempDivide);
  float xmove=remainPercent*xstep;
  float x0=x-xmove;
  return x0;
}
float CalculateY0(float y,float ystep) {
  float tempDivide=y/ystep;
  float remainPercent=tempDivide-floor(tempDivide);
  float ymove=(1-remainPercent)*ystep;
  float y0=y+ymove;
  return y0;
}



