/**************************************************************

    Name:    MatrixExport.cpp

    Purpose: This header contains the code for writing
             the data out to file.

    Notes: 

    Credits: Howard Heaton
             
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

    float PI = 3.1415926536;
    using namespace std;
#endif // __IncludedLibraries__

   
#ifndef Proton_Histories
#define Proton_Histories
    #include "GenerateHistories.cpp"
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



/**************************************************************

    Name:    string findFileName

    Purpose: This function searches for the any existing .txt
             files with a certain file name. If they exist,
             then a suffixed version number is incremented
             until a 'newest version' of the file can be
             created without overwriting existing data.

    Notes:   This saves from the tedious task of having to
             delete or rename files upon every execution so
             that data is not being overwritten.
             
**************************************************************/
string findFileName(string FolderName, string  FileName)
{   
    int  version = 0; 
    string fileToOpen;

    do{   // Loop until a nonexisting filename is found.
          version++;  // Increment # of version for each new file.

          //Concatenate the file name with the version # and appropriate suffix.
          fileToOpen = FolderName  + FileName +  concatenate("Version", version) + ".txt";

    }while(fileExists(fileToOpen)); // Loop until a new filename is made.
    
    cout << "\nSaving data to "<<fileToOpen<<".\n";
    return fileToOpen;
}


/**************************************************************

    Name:    void ExportMicahFormat

    Purpose: This function will output the system matrix data
             in the format originally proposed by Micah.

**************************************************************/

void ExportMicahFormat(string MatrixFolderName, int   histories, int    totalProjAngles,
                float  startAngle,     float angleSpace,float* Xtrue,
                int    yparts,         int   xparts,    float  boxLength,
                float  boxWidth,       neo   MyNeo,     int    pathOption,
                int    sliceIndex)
{
  /************* BEGIN WRITING SYSTEM MATRIX TO FILE *************/
 
  fstream outfile;  // Variable used to write System Matrix to file
 
  string fileName = findFileName(MatrixFolderName, "MicahFormat");
  
  outfile.open(fileName.c_str(),ios::out|ios::app);   //Open Matrix file

  int rowIndex=0;

  int rowsPerAngle=ceil(histories/(float)totalProjAngles);
  
  float endAngle=startAngle+totalProjAngles*angleSpace;
 
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
  outfile.close();
}





/**************************************************************

    Name:    void exportCSRandWEPL

    Purpose: This function will output the system matrix data
             in CSR format and the WEPL vector in a separate
             text file.

**************************************************************/

void    ExportCSRandWEPL(float   thetaStart,       float deltaTheta,    string MatrixFolderName,
                         string  MeasureVecFolder, string WEPLfileName, string CSRfileName,
                         int     totalProjAngles,  int    histories,    int    xparts,     
                         int     yparts,           int    pathOption,   float  boxWidth,   
                         neo     MyNeo,            float* Xtrue,        float  boxLength,  
                         int     sliceIndex)
{
    int rowIndex = 0;   // Variable used to indicate the current row
                        // of the system matrix.

    // The nonzero matrix entries will be stored in the 'nonzero' vector
    // and the 'WEPLvec' will store the Measurement Vector
    vector<float> nonzero, WEPLvec;

    // These vectors are used for giving the column indices of nonzero entries
    // and a pointer to the start of each row.
    vector<int>   columns, rowPointer;
                        
    int rowsPerAngle = ceil(histories/(float)totalProjAngles);

    float thetaStop = thetaStart + totalProjAngles * deltaTheta;
 
    //Convert the starting and ending gantry angles to radians
    thetaStart = thetaStart * PI / 180.0;
    thetaStop  = thetaStop  * PI / 180.0;
    deltaTheta = deltaTheta * PI / 180.0;

    for(float theta = thetaStart; theta < thetaStop; theta += deltaTheta)
    {
        // Apparently, there is some error if theta is a multiple of PI / 2.
        if(theta == 0 or theta == PI or theta == PI/2 or theta == (3*PI / 2))
            theta += PI / 180.0;

        // At each gantry angle multiple histories will occur. Here we loop
        // over each history generated at the current angle 'theta'.
        for(int r = 0; r < rowsPerAngle; r++)
        {   
            if(rowIndex < histories) // This provides a check to make sure the # of
            {                        // histories is not exceeded.
                
                // The following line generates a pointer for a proton history
                History2d* T = generateProtonPath2d_KnownHull(theta,Xtrue,yparts,
                                                                xparts,boxLength,boxWidth,
                                                                MyNeo,pathOption,sliceIndex);

                if(T->colInd.size() != 0)  //This is a check to ensure a history was generated
                {
                    // Loop over each of the basis functions intersected, which is
                    // given by the indices contained in the 'colInd' array.
                    for(int i = 0; i < T->colInd.size(); i++)
                    {
                        // Push each nonzero matrix entry onto the nonzero vector.
                       // nonzero.push_back(T->value[i]);
                          nonzero.push_back(1.00);

                            
                        // List the corresponding  # of the basis function intersected.
                        columns.push_back(T->colInd[i]);
                    }
                    
                    // Increment the rowIndex by the number of nonzero entries in
                    // the current history.
                    rowIndex += T->colInd.size();

                    // rowPointer's next entry indicateds where the next row of
                    // the matrix begins (according to the number of nonzero entries).
                    rowPointer.push_back(rowIndex);

                    // Add the current vector's WEPL value to the measurement vector.
                    WEPLvec.push_back(T->wepl);
                }
                else        // If a valid proton history was not generated, then the
                    r--;    // current iteration for the rows per angle must redone.

            }

        }

    }

    fstream CSRfile;    // CSR file is used to export the Syste matrix in CSR form.
   
        
    int version = 0; // Used for determining the version of output,
                       // in case the program has already been executed.

    string fileName = findFileName(MatrixFolderName, "CSRformat");
  
    CSRfile.open(fileName.c_str(),ios::out|ios::app);   //Open Matrix file

    //First, write out each of the nonzero matrix entry values
    for(int i = 0; i < nonzero.size(); i++)
        CSRfile << nonzero[i] << " ";

    CSRfile << endl;

    //Second, on the next line list the corresponding column entries
    for(int i = 0; i < columns.size(); i++)
        CSRfile << columns[i] << " ";

    CSRfile << endl;

    //Last, on the 3rd line list the indices corresponding to each row
    for(int i = 0; i < rowPointer.size(); i++)
        CSRfile << rowPointer[i] << " ";

    CSRfile << endl;

    CSRfile.close(); //Close the CSR file

    fstream WEPLfile;   // Generate variable for the WEPL measurement vector file.

    fileName = findFileName(MeasureVecFolder, "WEPLdata");
  
    WEPLfile.open(fileName.c_str(),ios::out|ios::app);   //Open Matrix file
    
    // List each WEPL history, separated by a space.
    for(int i = 0; i < WEPLvec.size(); i++)
        WEPLfile << WEPLvec[i] << " ";

    WEPLfile.close();   //Close WEPL file.

}



/**************************************************************

    Name:    void TempExport

    Purpose: This is a temporary function to transport the
             functionality of exporting the system matrix
             to file.

**************************************************************/

void TempExport(string MatrixFolder, string MeasureVecFolder,  int   histories, int    totalProjAngles,
                float  startAngle,     float angleSpace,float* Xtrue,
                int    yparts,         int   xparts,    float  boxLength,
                float  boxWidth,       neo   MyNeo,     int    pathOption,
                int    sliceIndex)
{

  ExportMicahFormat(MatrixFolder, histories, totalProjAngles, startAngle,
             angleSpace,     Xtrue,     yparts,          xparts,
             boxLength,      boxWidth,  MyNeo,           pathOption,
             sliceIndex);

  ExportCSRandWEPL(startAngle, angleSpace, MatrixFolder, MeasureVecFolder,
                         "WEPLdata",  "CSRmatrix",  totalProjAngles,
                         histories,    xparts,      yparts,
                         pathOption,   boxWidth,    MyNeo,
                         Xtrue,        boxLength,   sliceIndex);
}
