/**************************************************************

    Name:    ConfigSetup.cpp

    Purpose: This file is used for initially reading the
             configuration parameters from a file.

    Credits: Micah ______, Howard Heaton
             
    Micah originally wrote this code during his work for
    his Master's thesis. Then extensive commenting and
    appropriate modifications have been made by Howard
    Heaton to incorporate blob basis functions in the
    system matrix.

    Revision Date: July 17, 2014

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

#ifndef  Data_Export
#define  Data_Export
    #include "MatrixExport.h"  // This file is used for writing
#endif  // Data_Export         // information to a .txt file.

using namespace std;

/**************************************************************

    Name:    class ConfigObj

    Purpose: This is an object for reading data from the
             config files??

**************************************************************/

class ConfigObj
{
    public:
        // There are many variables used to initialize this object,
        // which are explained in more detail below.
        ConfigObj(bool d3,int xp,int yp,int zp,float bw,float bl,
            float bh,int h,int tpa,int sa,int as,string sfn,int po1,
            int po2) :
            dimensional3(d3),xpartitions(xp),ypartitions(yp),
            zpartitions(zp),boxwidth(bw),boxlength(bl),
            boxheight(bh),histories(h),totalprojectionangles(tpa),
            startangle(sa),anglespace(as),simulationfoldername(sfn),
            pathoption(po1),phantomoption(po2) {};

        bool    dimensional3;
        int     xpartitions;            // # of partitions in each of
        int     ypartitions;            // the 3 Cartesian directions
        int     zpartitions;
        float   boxwidth;               // Dimensions of each voxel/box
        float   boxlength;
        float   boxheight;
        int     histories;              // # of proton histories
        int     totalprojectionangles;  // # of angles used in scan
        int     startangle;             // Initial angle of gantry
        int     anglespace;
        string  simulationfoldername;   // Folder name to put data in after
        int     phantomoption;          // ? What type of phantom?
        int     pathoption;             // 1 = straight-line approx OR
                                        // 2 = Choose cubic spline
};

/**************************************************************

    Name:    ConfigObj ReadConfigObjFromFile()

    Purpose: This function takes inputs from the "configfile"
             and passes them appropriately into a Config Object.

**************************************************************/
ConfigObj ReadConfigObjFromFile()
{
    fstream infile;
    infile.open("configfile",ios::in);    //Open config file

    /* The following temporary variables are defined
       to be passed into the configuration object. */

    bool d3;
    int xp,yp,zp,h,tpa,sa,as,po1,po2;
    float bw,bl,bh;

    string sfn,tempS;

    //The following reads successive strings from file??

    infile>>tempS>>d3;    // Boolean: True if 3 dimensional
    infile>>tempS>>xp;    // # of x partitions
    infile>>tempS>>yp;    // # of y partitions
    infile>>tempS>>zp;    // # of z partitions
    infile>>tempS>>bw;    // box width
    infile>>tempS>>bl;    // box length
    infile>>tempS>>bh;    // box height
    infile>>tempS>>h;     // # of histories
    infile>>tempS>>tpa;   // # of projection angles during scan
    infile>>tempS>>sa;    // Starting angle of scan
    infile>>tempS>>as;
    infile>>tempS>>sfn;   // Simulation folder name
    infile>>tempS>>po1;   // Path option
    infile>>tempS>>po2;   // Phantom option

    infile.close();       // Close file

    //Pass parameters into a configuration object and return.
    return ConfigObj(d3,xp,yp,zp,bw,bl,bh,h,tpa,sa,as,sfn,po1,po2);
}


/**************************************************************

    Name:    ConfigObj ReadConfigObjFromFile()

    Purpose: This function takes inputs from the "configfile"
             and passes them appropriately into a Config Object.

**************************************************************/
string ReadEntryStringFromConfigFile(int entryId)
{
    fstream infile;
    infile.open("configfile",ios::in);
    string temp1,temp2;

    for (int i=0; i<entryId; i++)
        infile>>temp1>>temp2; //move to correct entry

    infile>>temp1>>temp2;
    return temp2;
}

