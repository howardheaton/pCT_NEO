/**************************************************************

    Name:    SimulatePct.cpp

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

    Revision Date: July 16, 2014

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
    #include "ProtonHistoryComputation.cpp"
#endif



/**************************************************************

    Name:    bool fileExists

    Purpose: This function checks to see if a given file
             already exists. If so, it returns true (and false
             otherwise).

**************************************************************/
bool fileExists(string fileName)
{   
    // Convert inputted string to a character array
    char* test = (char*)fileName.c_str(); 

    ifstream infile(test);  // Create input file
    
    return infile.good();   // Check and Return whether file exists
}



/**************************************************************

    Name:    string Concatenate

    Purpose: This function appends an integer to the end of
             a string and returns the formed concatenation.

**************************************************************/
string concatenate(string name, int num)
{
    stringstream newString;    

    newString << name << num;

    return newString.str();
}










// 2dHistories -> no proton scattering in the vertical direction
void SimulatePct_KnownHull_2dHistories(neo MyNeo,int histories,
 int zparts,int yparts,int xparts,float boxLength,
 float boxWidth,int totalProjAngles,float startAngle,
 float angleSpace,int pathOption,int phantomOption,
 string MatrixFolderName,string NeoXtrueFilename) {
  fstream outfileX;
  outfileX.open(NeoXtrueFilename.c_str(),ios::out); //clobber old
  float* Xtrue=CreateNeoSlice(MyNeo,yparts,xparts,boxLength,
   boxWidth,phantomOption);
  for (int i=0; i<yparts*xparts; i++)
    outfileX<<Xtrue[i]<<' ';
  outfileX<<endl;
  outfileX.close();
  int rowsPerAngle=ceil(histories/(float)totalProjAngles);
  int rowsPerAnglePerSlice=ceil(rowsPerAngle/(float)zparts);
  float endAngle=startAngle+totalProjAngles*angleSpace;
  for (float theta=startAngle; theta<endAngle; theta+=angleSpace) {
    fstream outfile;
    int tempTheta=(int)theta;
    ostringstream o;
    o<<tempTheta;
    string num=o.str();
    string filenameMat=MatrixFolderName+"AcsrB_"+num+".txt";
    outfile.open(filenameMat.c_str(),ios::out|ios::app);
    float thetaRad=theta*PI/180.0;
    for (int sliceId=0; sliceId<zparts; sliceId++ ) {
      for (int localRow=0; localRow<rowsPerAnglePerSlice; localRow++) {
        History2d* T=generateProtonPath2d_KnownHull(thetaRad,Xtrue,
         yparts,xparts,boxLength,boxWidth,MyNeo,pathOption,sliceId);
        if (T->colInd.size()!=0) {
          outfile<<tempTheta<<' ';
          outfile<<localRow<<' ';
          outfile<<T->wepl<<' ';
          outfile<<T->colInd.size()<<' ';
          for (int j=0; j<T->colInd.size(); j++)
            outfile<<T->colInd[j]<<' ';
          outfile<<endl;
          //if (localRow%100000==0)
          //cout<<"angle="<<theta<<" || localRow="<<localRow<<endl;
        }
        else {
          localRow--;
          //cout<<"localRow--"<<localRow<<endl;
        }
        delete T;
      }
    }
    cout<<"file finished for angle "<<theta<<endl;
    outfile.close();
  }
}

// File Format:
// rowIndex angle B hits colInd[0] colInd[1] ... colInd[hits-1] endl
// startAngle,angleSpace in degrees
void SimulatePct2d(neo MyNeo,int histories,int yparts,int xparts,
                   float boxLength,float boxWidth,int totalProjAngles,
                   float startAngle,float angleSpace,int pathOption,
                   int phantomOption,string MatrixFilename,
                   string NeoXtrueFilename)
{
  int sliceIndex=0; // 2d sim has only 1 slice

  /***** FOLLOWING SECTION IS FOR WRITING TRUE X VECTOR TO FILE*******/
      fstream outfileX; // Variable used to write the true X vector to file

      outfileX.open(NeoXtrueFilename.c_str(),ios::out); // Open X vector file
      if(not outfileX.is_open())
            cout << "ERROR 1: SimulatePct2d unable to open NeoXtrueFilename.c_str()";
            

      float* Xtrue = CreateNeoSlice(MyNeo,yparts,xparts,boxLength,
                                    boxWidth,phantomOption);
      //Output the true image vector X
      for (int i=0; i<yparts*xparts; i++)
        outfileX<<Xtrue[i]<<' ';
     
      outfileX<<endl;
     
      outfileX.close(); //Close file for true X vector

  /*************** WRITING X VECTOR TO FILE COMPLETE *************/


  /************* BEGIN WRITING SYSTEM MATRIX TO FILE *************/
 
  fstream outfile;  // Variable used to write System Matrix to file
  
  int version = 0; // Used for determining the version of output,
                   // in case the program has already been executed.

  string fileToOpen;
  
  cout << "\nChecking for latest output version...\n";
    
  do{   // Loop until a nonexisting filename is found.
        version++;  // Increment # of version for each new file.

        //Concatenate the file name with the version # and appropriate suffix.
        fileToOpen = MatrixFilename +  concatenate("Version", version) + ".txt";

  }while(fileExists(fileToOpen)); // Loop until a new filename is made.

  outfile.open(fileToOpen.c_str(),ios::out|ios::app);   //Open Matrix file

  cout << "\nSaving data to "<<fileToOpen << ".\n";

  int rowIndex=0;

  int rowsPerAngle=ceil(histories/(float)totalProjAngles);
  
  float endAngle=startAngle+totalProjAngles*angleSpace;
 
  cout << "\nBeginning Virtual Proton Projections\n";

  //Cycle through each of the angles of the gantry
  for (float theta=startAngle; theta<endAngle; theta+=angleSpace)
  {
    
      float thetaRad=theta*PI/180.0;    //Convert to Radians
    
      if (equalf(thetaRad,0.0) or equalf(thetaRad,PI) or
          equalf(thetaRad,PI/2) or equalf(thetaRad,3/2*PI))
             thetaRad+=1*PI/180.0;  // to prevent division by zero (this will be fixed)

      // The # of proton histories per gantry angle is given 'rowsPerAngle'
      for (int localRow=0; localRow<rowsPerAngle; localRow++)
      {
         if (rowIndex<histories)    //1 row is dedicated to each history
         {
             History2d* T=generateProtonPath2d_KnownHull(thetaRad,Xtrue,yparts,
                                                         xparts,boxLength,boxWidth,
                                                         MyNeo,pathOption,sliceIndex);

             if (T->colInd.size()!=0) //If a history was generated?
             {
                 //Identify Matrix Row & Gantry Angle
                 outfile<<rowIndex<<' '<<thetaRad*180/PI<<' ';

                 //Identify WEPL Value and the label # of voxels that were intersected
                 outfile<<T->wepl<<' '<<T->colInd.size()<<' ';

                 //Cycle  through and give #s of each voxel intersected
                 for (int j=0; j<T->colInd.size(); j++)
                     outfile<<T->colInd[j]<<' ';
        
                 outfile<<endl; //End the current row of the matrix
                 rowIndex++;    //Index Row
            }
            else 
              localRow--;
        }
    }
  }
  cout << "\nFinished Virtual Proton Projections\n";
  outfile.close();
}

