/**************************************************************

    Name:    ProtonHistoryComputation.cpp

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

#endif // __IncludedLibraries_

#ifndef Path_Computation
#define Path_Computation
        #include "ProtonPathComputation.cpp"
#endif  // Path_Computation


/**************************************************************

    Name:    class History2d

    Purpose: This object contains the WEPL value for a proton
             history and colInd contains a record of each of
             the basis functions (voxels) intersected by
             the given proton

**************************************************************/
class History2d {
public:
  History2d(vector<int> c,float w)
   : colInd(c),wepl(w) {};
  vector<int> colInd;
  float wepl;
};




// KnownHull indicates that we intersected the proton path with the ellipse boundary of the neo phantom
History2d* generateProtonPath2d_KnownHull(float theta,float* Xtrue,
        int yparts,int xparts,float boxLength,float boxWidth,neo MyNeo,
        int pathOption,int sliceIndex)
{
  // Multiply # of voxels per slice by the # of the current slice.
  // This will give the current starting index for ___ entries.
  int startCol=sliceIndex*yparts*xparts;

  vector<int> colInd,nullvec;
  float wepl=0;
  int col,localCol;
  PointVec Pin,Pt,Pout;
  float thetaOut;
  float a=MyNeo.E[MyNeo.E.size()-1].a; // semi major xaxis length of outer ellipse
  float b=MyNeo.E[MyNeo.E.size()-1].b; // semi minor yaxis ...
  float ystep=boxLength/(float)yparts;
  float xstep=boxWidth/(float)xparts;
  // the shifts to center the recon space
  float yshift=boxLength/2;
  float xshift=boxWidth/2;
  float gantryRad=3000;
  vector<PointVec3> scatParam=scatteringParameters();
  PointVec P1(gantryRad*cos(theta),gantryRad*sin(theta));
  PointVec v1(cos(theta),sin(theta));
  PointVec v2(-sin(theta),cos(theta));
  float varLat,covariance,varAng;
  float d1=MyUniformRand();
  //float d1=0; /**************random********/
  //calculate intersection of proton entry line with ellipse
  PointVec P3(P1.x+d1*v2.x,P1.y+d1*v2.y);
  float m=v1.y/v1.x; //slope of proton entry line
  float B=P3.y-m*P3.x; //y-intercept of proton entry line
  //to calculate the entry point Pin:
  //the intersection test of line and ellipse --> solve for x in
  //the quadratic equation:  alpha*x^2+beta*x+gamma=0
  float alpha=1/pow(a,2)+pow(m/b,2);
  float beta=2*m*B/pow(b,2);
  float gamma=pow(B/b,2)-1;
  float disc=pow(beta,2)-4*alpha*gamma; //discriminant
  if ( disc<0 or equalf(disc,0) ) {
    // proton path missed the object -> return blank history
    History2d* T=new History2d(nullvec,0);
    return T;
  }
  else {
    float x1=(-beta+sqrt(disc))/(2*alpha);
    float x2=(-beta-sqrt(disc))/(2*alpha);
    float y1=m*x1+B;
    float y2=m*x2+B;
    float y11=b*sqrt(1-pow(x1/a,2));
    float y22=-b*sqrt(1-pow(x2/a,2));
    float norm1=normEu(PointVec(x1,y1),P3);
    float norm2=normEu(PointVec(x2,y2),P3);
    if (norm1<norm2) {
      Pin.set(x1,y1);
      Pt.set(x2,y2);
    }
    else {
      Pin.set(x2,y2);
      Pt.set(x1,y1);
    }
    //depth is the distance proton traveled through object
    //assuming it moves in a straight line
    float depth=normEu(Pin,Pt);
    //the lookup table is indexed by cm so divide by 10
    //to convert mm to cm
    int index=ceil(depth/10);
    if (index>20) index=20; //avoid seg faults
    varLat=scatParam[index].x;
    covariance=scatParam[index].y;
    varAng=scatParam[index].z;
    PointVec R2=generateJointRandNorm(varLat,covariance,varAng);
    float d2=R2.x; //exiting lateral displacement
    float psi=R2.y; //exiting angular displacement
    //generate temporary line to intersect with ellipse to find Pout
    PointVec P4(Pt.x+d2*v2.x,Pt.y+d2*v2.y);
    float m2=m; //slope
    float B2=P4.y-m2*P4.x; //intercept
    //line-ellipse interesection quadratic coefficients
    float alpha2=1/pow(a,2)+pow(m2/b,2);
    float beta2=2*m2*B2/pow(b,2);
    float gamma2=pow(B2/b,2)-1;
    float disc2=pow(beta2,2)-4*alpha2*gamma2; //discriminant
    if (disc2<0) { // very rare circumstance
      History2d* T=new History2d(nullvec,0);
      return T;
    }
    else if (equalf(disc2,0)) {
      float x21=-beta2/(2*alpha2);
      float y21=m2*x21+B2;
      Pout.set(x21,y21);
    }
    else {
      float x21=(-beta2+sqrt(disc2))/(2*alpha2);
      float x22=(-beta2-sqrt(disc2))/(2*alpha2);
      float y21=m2*x21+B2;
      float y22=m2*x22+B2;
      float norm21=normEu(PointVec(x21,y21),P4);
      float norm22=normEu(PointVec(x22,y22),P4);
      if (norm21<norm22)
        Pout.set(x21,y21);
      else
        Pout.set(x22,y22);
    }
    //exiting angle relative to xy axis
    thetaOut=theta+psi;
    //calculate upper left corner coordinates of the entry voxel
    //and exit voxel to find xmin and ymin
    float xin0=CalculateX0(Pin.x,xstep);
    float xout0=CalculateX0(Pout.x,xstep);
    float xmin=min(xin0,xout0);
    float xmax=max(xin0,xout0);
    if (pathOption==1) { //straight-line between Pin and Pout
      float ms=(Pout.y-Pin.y)/(Pout.x-Pin.x); //slope
      float Bs=Pout.y-ms*Pout.x; //y-intercept
      for (float xx=xmin; xx<xmax; xx+=xstep) {
        float ys=ms*xx+Bs;
        float y0=CalculateY0(ys,ystep);
        int i=(int)((yshift-y0)/ystep);
        int j=(int)((xshift+xx)/xstep);
        //if the slope of the line is greater than 1 in abs val
        //then the line will intersect more voxels above i (if
        //positive slope) before the next xx so need to mark all
        //those voxels by calculating y at the next xx
        float ysnext=ms*(xx+xstep)+Bs;
        float y0next=CalculateY0(ysnext,ystep);
        int inext=(int)((yshift-y0next)/ystep);
        if (i>=0 and inext>=0) {
          for (int ii=min(i,inext); ii<=max(i,inext); ii++) {
            localCol=ii*xparts+j;
            if (0<localCol and localCol<xparts*yparts)
            {    // LocalCol gives the voxel in the plane
                 // Need to add the startCol to get the appropriate
                 // value taking into account the # of the slice
                 col=startCol+localCol;
                
                 colInd.push_back(col); //Save the # of the voxel intersected
                 wepl+=Xtrue[localCol]; //Sum the WEPL values for each pixel entered by proton
                                        //NOTE: THE LENGTH OF THE INTERSECTION IS NOT INCLUDED
                                        //      AND IS THUS ASSUMED TO BE 1!!!
            }
          }
        }
      }
    }
    else if (pathOption==2) //cubic spline
    { 
      float thetaIn=theta;
      float c=thetaIn*(Pout.x-Pin.x)-(Pout.y-Pin.x);
      float d=-thetaOut*(Pout.x-Pin.x)-(Pout.y-Pin.y);
      for (float xx=xmin; xx<xmax; xx+=xstep)
      {
        float t=(xx-Pin.x)/(Pout.x-Pin.x);
        float q=(1-t)*Pin.y+t*Pout.y+t*(1-t)*(c*(1-t)+d*t);
        float y0=CalculateY0(q,ystep);
       
        int i=(int)((yshift-y0)/ystep);
        int j=(int)((xshift+xx)/xstep);
        
        float tnext=(xx+xstep-Pin.x)/(Pout.x-Pin.x);
        float qnext=(1-tnext)*Pin.y+tnext*Pout.y+tnext*(1-tnext)*(c*(1-tnext)+d*tnext);
        float y0next=CalculateY0(qnext,ystep);

        int inext=(int)((yshift-y0next)/ystep);

        if (i>=0 and inext>=0)
        {
          for (int ii=min(i,inext); ii<=max(i,inext); ii++)
          {
            localCol=ii*xparts+j;
            if (0<=localCol and localCol<xparts*yparts)
            {
              col=startCol+localCol;
              colInd.push_back(col);
              wepl+=Xtrue[localCol];
            }
          }
        }
      }
    }
  }
  History2d* T=new History2d(colInd,wepl);
  return T;
}
