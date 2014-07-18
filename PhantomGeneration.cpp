/**************************************************************

    Name:    PhantomGeneration.cpp

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

/**************************************************************
***************************************************************

    THE REMAINING CODE IS DEDICATED TO PHANTOM GENERATION

**************************************************************
*************************************************************/



/**************************************************************

    Name:    class Ellipse

    Purpose: Object to give the well-known properties related
             to ellipses.

    Notes:   The prescribed object gives an ellipse centered
             at (xc, yc) with major axis length 2a and semi-
             major axis length 2b.

**************************************************************/
class Ellipse {
    public:
        //Initialize Ellipse values
        Ellipse(float xxc,float yyc,float aa,float bb) :
        xc(xxc),yc(yyc),a(aa),b(bb) {};

        float xc;     // X coordinate center
        float yc;     // Y coordinate center
        float a;      // Major Axis
        float b;      // Semi-Major Axis
};


/**************************************************************

    Name:    class NEO

    Purpose: NEO contains vectors of ellipses to describe
             convex objects

    Notes:   NEO = Nonuniform Ellipse Object phantom
             Last element of "V?" must be the boundary ellipse
             Inner ellipses must come before larger/outer
             ellipses in the E vector

**************************************************************/

class neo {
    public:
      neo(vector<Ellipse> e,vector<float> r) : E(e),rsp(r) {};
      vector<Ellipse> E;    // E   = Ellipses
      vector<float> rsp;    // rsp = Relative Stopping Power
};


/**************************************************************

    Name:    new CreateDefaultNeo

    Purpose: This function generates a standard NEO phantom

    Notes:   The float 's' is used to give.... RELATIVE SIZING?

**************************************************************/

neo CreateDefaultNeo(float s) {
    vector<Ellipse> E;
    vector<float>   rsp;

    // Ellipse Parameters -> (x, y, a, b)
    Ellipse e1(s*0  ,s*0 , s*70, s*90);
    Ellipse e2(s*0  ,s*85, s*10, s*2.5);
    Ellipse e3(s*0  ,s*0 , s*60, s*80);
    Ellipse e4(s*20 ,s*0 , s*10, s*20);
    Ellipse e5(s*-20,s*0 , s*10, s*20);

    E.push_back(e5);
    E.push_back(e4);
    E.push_back(e3);
    E.push_back(e2);
    E.push_back(e1);

    rsp.push_back(0.9); // ventricles
    rsp.push_back(0.9); // ventricles
    rsp.push_back(1.04);// brain tissue
    rsp.push_back(0.0); // air pocket
    rsp.push_back(1.6); // skull

    return neo(E,rsp);
}

/**************************************************************

      Name:    class LPoint

      Purpose: Generates a float (x,y) point in the plane

      Notes:   What is 'l' useful for???
.x
**************************************************************/
class LPoint {
 public:
  LPoint() {
    x=0;
    y=0;
    l=0;
  }
  LPoint(float xx,float yy,int ll) :
    x(xx),y(yy),l(ll) {};
  float x;
  float y;
  int l;
};


/**************************************************************

        Name:    float calculateTriangleArea

        Purpose: This one is rather self-explanatory

        Notes:   See the link below for more info on
                 how this computation works.
                 http://www.mathopenref.com/coordtrianglearea.html

**************************************************************/
float calculateTriangleArea(LPoint A,LPoint B, LPoint C) 
{
  float output=fabs((A.x*(B.y-C.y)+B.x*(C.y-A.y)+C.x*(A.y-B.y))/2.0);
  return output;
}


/* Simple equality function based upon set tolerance 'tol'.  */
bool equalf(float a,float b) {
  float tol=0.00001;
  if (a<b+tol && a>b-tol)
    return true;
  else
    return false;
}


/**************************************************************

          Name:    bool IsPointInEllipse

          Purpose: This one is rather self-explanatory

          Notes:   (k,h) gives the center of Ellipse of
                   interest. a and b are the respective
                   values for  the major and semi-major axis.
                   (xc, yc) is the point under investigation.

**************************************************************/
bool IsPointInEllipse(float xc,float yc,float a,float b, float k,float h)
{
  /* The core of this function is the equation for an ellipse:
     (x-k)^2 / a^2 + (y-h)^2 / b^2 = 1.
     Note that this is the reverse  of the normal notation where
     (h,k) denotes the center.Regardless, rearraging yields
     y = yc +- b sqrt(1 - (x-xc)^2 / a^2).
     So, checking to see if the passed value for yc is between
     the upper and lower bounds  of y would determine whether
     (xc,yc) is in the ellipse. (Note that in order to avoid
     a negative square root we must have |xc-k| <= a.
     Currently, no error handling is present to enfore this
     condition...   
   */
  float yup=h+b*sqrt(1.0-pow(xc-k,2)/pow(a,2)); //upper bound
  float ylo=h-b*sqrt(1.0-pow(xc-k,2)/pow(a,2)); //lower bound
 
  if (yc<=yup && yc>=ylo) //check for yc between bounds
    return true;
  else
    return false;
}

/**************************************************************

    Name:    float FindRSPatPoint

    Purpose: This function determines the relative
             stopping power RSP at a given point
             (x,y) inside the NEO phantom

    Notes:   The incorporated 'for' loop is where
             it is evident that the ellipses must
             be listed in the phantom from the most
             inner to outermost ellipse.
                     
***************************************************************/

float FindRSPatPoint(float x,float y,neo phant)
{
  //Temporary variables for iterations through phantom  
  float RSP=0;     
  float xcPhan,ycPhan,aPhan,bPhan;

  //Cycle through array of ellipses composing the phantom
  for (int i=0; i<phant.E.size(); i++)
  {  
    xcPhan=phant.E[i].xc;
    ycPhan=phant.E[i].yc;
    aPhan=phant.E[i].a;
    bPhan=phant.E[i].b;
    /*Note: instead of making new temporary variables
            and passing their values into the following
            'if' comparison, time could be saved by
            just referencing each ellipse "phant.E[i]".
    */

    if (IsPointInEllipse(x,y,aPhan,bPhan,xcPhan,ycPhan))
    {
      RSP=phant.rsp[i];
      return RSP;
    }
  }
  // If all of the Ellipses have been cycled through,
  // then return the initialized value 0 (hence the
  // point is outside the object).
  return RSP;
}

/**************************************************************

      Name:    float FindRSPatPoint_DefaultNeo

      Purpose: This is simply a hardwiring of the above
               function in the defualt

**************************************************************/

float FindRSPatPoint_DefaultNeo(float xc,float yc) {
  float RSP;
  if (IsPointInEllipse(xc,yc,10,20,-20,0) ||
  (IsPointInEllipse(xc,yc,10,20,20,0)))
    RSP=0.9;
  else if (IsPointInEllipse(xc,yc,60,80,0,0)) {
    RSP=1.04;
  }
  else if (IsPointInEllipse(xc,yc,10,2.5,0,85))
    RSP=0;
  else if (IsPointInEllipse(xc,yc,70,90,0,0))
    RSP=1.6;
  else
    RSP=0;
  return RSP;
}


/**********************************************************

    Name:       float* CreateNeoSlice

    Purpose:    This function is used to generate a float
                version of the true X vector?

**********************************************************/

// method==0 --> centerPointMethod      vs
// method==1 --> Corner Point Averaging
float* CreateNeoSlice(neo phantom,int ypartitions,int xpartitions,
                      float boxLength,float boxWidth,int method)
{
  int sizeXtrue=xpartitions*ypartitions*sizeof(float);
  float *Xtrue=(float*)malloc(sizeXtrue);
  float yshift=boxLength/2.0;
  float xshift=boxWidth/2.0;
  float ystep=boxLength/ypartitions;
  float xstep=boxWidth/xpartitions;
  float x0,x1,x2,x3,x4,y0,y1,y2,y3,y4;
  float xc,yc,RSP,temp0,temp1,temp2,temp3,temp4;
  float xdot1l0,xdot2l0,xdot1l2,xdot2l2,
        ydot1l1,ydot2l1,ydot1l3,ydot2l3;
  float base,rbase,height,rheight,area1,area2;
  float a,b;
  float voxelArea=xstep*ystep;
  int counter1=0;

  //Center Point Method
  if (method==0)
  { 
    for (int i=0; i < ypartitions; i++)
    { for (int j=0; j < xpartitions; j++)
      {
        x1=j*xstep-xshift;      // set to left-hand side of pixel
        x2=x1+xstep;            // set to right side of pixel
        y1=yshift-i*ystep;      // top of pixel
        y2=y1-ystep;            // bottom of pixel
        xc=(x1+x2)/2.0;         // horizontal center of pixel
        yc=(y1+y2)/2.0;         // vertical center of pixel

        //Find RSP at center point (xc,yc)
        RSP=FindRSPatPoint(xc,yc,phantom);

        //Set true vector to centered RSP value
        Xtrue[i*xpartitions+j]=RSP;
      }                            
    }
  }
  //Corner Point Averaging Method
  else
      if (method==1)
      {
        for (int i=0; i<ypartitions; i++)
        { for (int j=0; j<xpartitions; j++)
          {
            x1=j*xstep-xshift;  // set to left-hand side of pixel
            y1=yshift-i*ystep;  // top    of pixel
            x2=x1+xstep;        // right  of pixel
            y2=y1;              // top    of pixel
            x3=x2;              // right  of pixel
            y3=y1-ystep;        // bottom of pixel
            x4=x1;              // left   of pixel
            y4=y3;              // bottom of pixel
        
            //Find RSP Values at Each Corner
            temp1=FindRSPatPoint(x1,y1,phantom);//top-left
            temp2=FindRSPatPoint(x2,y2,phantom);//top-right
            temp3=FindRSPatPoint(x3,y3,phantom);//bottom-right
            temp4=FindRSPatPoint(x4,y4,phantom);//bottom-left
            
            //Average RSP Values
            RSP=(temp1+temp2+temp3+temp4)/4.0;
            Xtrue[i*xpartitions+j]=RSP;
          }
        }
      }
      else 
      {   cout << "ERROR: Phantom Option (variable name: 'method') is not a valid choice (i.e. 1 o r2)!\n";
          system("pause");
          exit(1);
      }

  return Xtrue;
}


float* CreateXtrueNeo3d(neo MyNeo,int zparts,int yparts,int xparts,
 float boxLength,float boxWidth,int method) {
  float* Neo2d=CreateNeoSlice(MyNeo,yparts,xparts,boxLength,boxWidth,method);
  float* XtrueNeo3d=(float*)malloc(zparts*yparts*xparts*sizeof(float));
  for (int k=0; k<zparts; k++) {
    for (int ii=0; ii<yparts*xparts; ii++) {
      XtrueNeo3d[k*yparts*xparts+ii]=Neo2d[ii];
    }
  }
  free(Neo2d);
  return XtrueNeo3d;
}



void WriteNeoToFile(float* Xtrue,int zs,int ys,int xs,string filename) {
  fstream outfile;
  outfile.open(filename.c_str(),ios::out);
  outfile<<zs<<' '<<ys<<' '<<xs<<endl;
  for (int ii=0; ii<zs*ys*xs; ii++)
    outfile<<Xtrue[ii]<<' ';
  outfile.close();
}
