/**************************************************************

    Name:    MatrixExport.h

    Purpose: This header contains the code for writing
             the data out to file.

    Notes: 

    Credits: Howard Heaton
             
    Revision Date: July 17, 2014

**************************************************************/
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

using namespace std;

float PI = 3.1415926536;



/**************************************************************

    Name:    void exportCSRandWEPL

    Purpose: This function will output the system matrix data
             in CSR format and the WEPL vector in a separate
             text file.

**************************************************************

void    exportCSRandWEPL(float  thetaStart,   float  thetaStop,   float deltaTheta,
                         string WEPLfileName, string CSRfileName, int   rowsPerAngle,
                         int    histories,    int    xparts,      int   yparts,
                         int    pathOption,   float  boxWidth,    neo   MyNeo,
                         float* Xtrue,        )
{
    int rowIndex = 0;   // Variable used to indicate the current row
                        // of the system matrix.

    // The nonzero matrix entries will be stored in the 'nonzero' vector
    // and the 'WEPLvec' will store the Measurement Vector
    vector<float> nonzero, WEPLvec;

    // These vectors are used for giving the column indices of nonzero entries
    // and a pointer to the start of each row.
    vector<int>   columns, rowPointer;
                        
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
                        nonzero.push_back(T->value[i]);

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
                    WEPLvec.push_back(T->WEPL);
                }
                else        // If a valid proton history was not generated, then the
                    r--;    // current iteration for the rows per angle must redone.

            }

        }

    }

    fstream CSRfile;    // CSR file is used to export the Syste matrix in CSR form.
    
    CSRfile.open(CSRfileName.c_str(), ios::out|ios::app);   //Open the CSR file.

    if(not CSRfile.open())      //Check to make sure the file was able to open.
        cout << "ERROR in exportCSRandWEPL: Unable to open CSR file.";
    
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

    cout << "\nCSR format System matrix created.\n"; // Output confirmation.
    
    fstream WEPLfile;   // Generate variable for the WEPL measurement vector file.

    WEPLfile.open(WEPLfileName.c_str(), ios::out|ios::app); //Open WEPL file.
    
    // List each WEPL history, separated by a space.
    for(int i = 0; i < WEPLvec.size(); i++)
        WEPLfile << WEPLvec[i] << " ";

    WEPLfile.close();   //Close WEPL file.

    cout << "\nWEPL Measurement Vector created.\n"; // Output confirmation.
}


*/
