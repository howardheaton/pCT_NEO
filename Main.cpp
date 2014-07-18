/**************************************************************

    Name:    Main.cpp

    Purpose: This program is used to generate phantoms for
             simulation purposes in pCT. It's strengths lie
             in its ability to generate custom phantoms of
             varying complexity.

    Notes:   In the context of this program the term "voxel"
             and "box" are synonymous and the later will
             typically be used.

    Credits: Micah ______, Howard Heaton
             
    Micah originally wrote this code during his work for
    his Master's thesis. Then extensive commenting and
    appropriate modifications have been made by Howard
    Heaton to incorporate blob basis functions in the
    system matrix.

    Revision Date: July 17, 2014

**************************************************************/

#ifndef Configuration_Setup
#define Configuration_Setup
    #include "ConfigSetup.cpp"
#endif

#ifndef Phantom_Generation
#define Phantom_Generation
    #include "PhantomGeneration.cpp" // As indicated, this file is used for
#endif                               // generating the NEO phantom.


#ifndef Path_Computation
#define Path_Computation
        #include "ProtonPathComputation.cpp"
#endif  // Path_Computation
   
#ifndef Proton_Histories
#define Proton_Histories
    #include "GenerateHistories.cpp"
#endif

#ifndef Simulate
#define Simulate
    #include "SimulatePct.cpp"
#endif

int main()
{    
  cout << "\nOpening SimulatePct.cpp\n";
  ConfigObj  C1             = ReadConfigObjFromFile();
  bool       threeDim       = C1.dimensional3; // 1
  int        xparts         = C1.xpartitions; // 160
  int        yparts         = C1.ypartitions;  // 200
  int        zparts         = C1.zpartitions; // 200
  float      boxWidth       = C1.boxwidth; // 160
  float      boxLength      = C1.boxlength; // 200
  float      boxHeight      = C1.boxheight; // 200
  int        histories      = C1.histories; // 64000000
  int        voxels         = C1.xpartitions*C1.ypartitions; //6400000
  int        totalProjAngles= C1.totalprojectionangles; // 180
  int        startAngle     = C1.startangle; // 1
  int        angleSpace     = C1.anglespace; // 2
  string     SimName        = C1.simulationfoldername;
  int        pathOption     = C1.pathoption; // 0 (straight-line)
  int        phantomOption  = C1.phantomoption; // 1 (corner-points)

  string     SimulationFolderName   = "SimulationData/" + SimName;
  string     MatrixFolderName       = SimulationFolderName + "/MatrixData/";
  string     NeoXtrueFilename       = SimulationFolderName + "/NEOXtrue.txt";
  string     MeasureVecFolder       = SimulationFolderName + "/DetectData/";
  neo        MyNeo                  = CreateDefaultNeo(1.0); // 1.0 is scale factor

  if (!threeDim)
  { // for 2d simulation
     
    cout << "\nExecuting 'SimulatePct2d'\n";
    string MatrixFilename=MatrixFolderName+"MicahFormat"; //This will be changed later...
    SimulatePct2d(MyNeo,histories,yparts,xparts,boxLength,
                  boxWidth,totalProjAngles,startAngle,angleSpace,
                  pathOption,phantomOption,MatrixFilename,
                  NeoXtrueFilename);
  }
  else
  { cout << "\nExecuting 'SimulatePct_KnownHull_2dHistories\n";
    SimulatePct_KnownHull_2dHistories(MyNeo,histories,zparts,yparts,xparts,boxLength,
                                      boxWidth,totalProjAngles,startAngle,angleSpace,
                                      pathOption,phantomOption,MatrixFolderName,
                                      NeoXtrueFilename);
  }
  cout << "\nSimulation Completed\n\n";
}
