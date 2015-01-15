/* Serial Code to re-arrange multi-file output from projection code into a single large file */

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

using namespace std;


double pcTOcm = 3.08567758e18;
double totMass = 0.0;
double totMom[3] = {0,0,0};


void ReadCell(FILE* readMe, double *Array)
{
    fread(&Array[0],sizeof(double),8,readMe);
}

void WriteCell(FILE* writeMe, double *Array)
{
    fwrite(&Array[0],sizeof(double),8,writeMe);
}


int main(int argc, char **argv)
{
    
    // Ready... Set... GO!
    clock_t timeMe;
    timeMe = clock();
    
    int n,i,j;
    
    if(argc < 2) // Double check on inputs
    {
        printf("WATERLOO: 1 input needed: (int) numFiles\nYou gave %d.\n",argc-1);
        exit(-1);
    }
    
    int numFiles = atoi(argv[1]);
    char InFileHeader[25] = "./gadget2orion_";
    
    // Open the set of files
    FILE *InFiles[numFiles];
    
    FILE *OutFile;
    char out_file[50];
    sprintf(out_file,"%s%s",InFileHeader,"all");

    for(n=0;n<numFiles;n++)
    {
        char cur_file[50];
        sprintf(cur_file,"%s%d",InFileHeader,n);
        
        std::ifstream myfile(cur_file);
        if(!myfile)
        {
            printf("WATERLOO: Trying to read from file that doesn't exist!\n");
            printf("Cannot find file %s\n",cur_file);
            exit(-1);
        }
        
        InFiles[n] = fopen(cur_file,"r");
        printf("Opened file %s\n",cur_file);
    }
    
    
    // Reads in Header information
    double width[numFiles];
    int ref_lev[numFiles];
    int totProc[numFiles];
    int Cords[numFiles][3];
    
    for(n=0;n<numFiles;n++)
    {
        fread(&totProc[n],sizeof(int),1,InFiles[n]);
        fread(&ref_lev[n],sizeof(int),1,InFiles[n]);
        fread(&width[n],sizeof(double),1,InFiles[n]);
        fread(&Cords[n],sizeof(int),3,InFiles[n]);
    }
    
    // Checks to make sure headers all agree
    for(n=0;n<numFiles;n++)
    {
        if(totProc[n] != numFiles)
        {
            printf("WATERLOO: Data files thinks there are %d processors. You think there are %d.\n",totProc[n],numFiles);
            exit(-1);
        }
        
        if(ref_lev[n] != ref_lev[0])
        {
            printf("WATERLOO: Data files disagree on ref_level. File 0 thinks it is %d while File %d thinks it is %d.\n",ref_lev[0],n,ref_lev[n]);
            exit(-1);
        }
        
        if(width[n] != width[0])
        {
            printf("WATERLOO: Data files disagree on width. File 0 thinks it is %g while File %d thinks it is %g.\n",width[0],n,width[n]);
            exit(-1);
        }
        
        for(j=0;j<3;j++)
        {
            if(Cords[n][j] != Cords[0][j])
            {
                printf("WATERLOO: Data files disagree on Cordinate distribution. File 0 thinks it is %d,%d,%d while File %d thinks it is %d,%d,%d.\n",Cords[0][0],Cords[0][1],Cords[0][2],n,Cords[n][0],Cords[n][1],Cords[n][2]);
                exit(-1);
            }
        }
    }
    
    printf("Grid Cells (ref_lev): %d.\n",ref_lev[0]);
    printf("Width: %g pc.\n",width[0]);
    printf("Cordinate distribution: %d,%d,%d .\n",Cords[0][0],Cords[0][1],Cords[0][2],n);
    
    // DeltaX
    double DeltaX = pcTOcm*width[0]/((double)ref_lev[0]);
    
    // Number of cells per file
    int xCells = ref_lev[0] / Cords[0][0];
    int yCells = ref_lev[0] / Cords[0][1];
    int zCells = ref_lev[0] / Cords[0][2];
    cout << "iCells : " << xCells << "," << yCells << "," << zCells << endl;
    
    int CellsPP = xCells*yCells*zCells ;
    int totCells = CellsPP * totProc[0];
    
    int CRanges[numFiles][3][2];
    for(i=0;i<numFiles;i++)
    {
        int xBlock = i % Cords[0][0];
        int yBlock = ((i-xBlock)/Cords[0][0]) % Cords[0][1];
        int zBlock = (i-xBlock-yBlock*Cords[0][0])/Cords[0][0]/Cords[0][1];
        
        CRanges[i][0][0] = xBlock*xCells;
        CRanges[i][0][1] = (1+xBlock)*xCells-1;
        CRanges[i][1][0] = yBlock*yCells;
        CRanges[i][1][1] = (1+yBlock)*yCells-1;
        CRanges[i][2][0] = zBlock*zCells;
        CRanges[i][2][1] = (1+zBlock)*zCells-1;
        printf("rank %d : %d,%d   %d,%d   %d,%d\n",i,CRanges[i][0][0],CRanges[i][0][1],CRanges[i][1][0],CRanges[i][1][1],CRanges[i][2][0],CRanges[i][2][1]);
       
    }
    
    // Open output file
    double curValues[8];
    OutFile = fopen(out_file,"w");
    
    // Print header info
    fwrite(&totProc[0], sizeof(int), 1, OutFile);
    fwrite(&ref_lev[0], sizeof(int), 1, OutFile);
    fwrite(&width[0], sizeof(double), 1, OutFile);
    
    int c=0;
    int ReadCount[numFiles];
    for(n=0;n<numFiles;n++) ReadCount[n] = 0;
    
    
    
    // This loop will read in all the data and calculate the center of mass velocity,
    // which it will subtract out.
    double vCOM[3] = {0,0,0};
    double pos[3] = {0,0,0};
    double angM[3] = {0,0,0};
    double RhoT = 0;
    
    while(1)
    {
        // Are we done?
        if(c==totCells) break;
        
        // Determine x,y,z cell values
        int CurCord[3];
        CurCord[0] = c % ref_lev[0];
        CurCord[1] = ((c-CurCord[0])/ref_lev[0]) % ref_lev[0];
        CurCord[2] = (c-CurCord[0]-CurCord[1]*ref_lev[0]) / pow(ref_lev[0],2);
        
      
        
        // Determines which cell range these values are in
        int PossibleFiles[numFiles];
        for(j=0;j<numFiles;j++) PossibleFiles[j] = j;
        
        // Narrows down
        for(j=0;j<3;j++)
            for(n=0;n<numFiles;n++)
            {
                if( CurCord[j] < CRanges[n][j][0] || CurCord[j] > CRanges[n][j][1] )
                {
                    PossibleFiles[n] = -1;
                }
            }
        
        // Now only one entry should be different from -1
        int fNumber=-1;
        for(n=0;n<numFiles;n++)
        {
            if(PossibleFiles[n] > -1 && fNumber==-1 ) fNumber = n;
            else if(PossibleFiles[n] > -1 && fNumber!=-1)
            {
                printf("WATERLOO: Cell %d exists in multiple files!\n",c);
                exit(0);
            }
            
        }
        if(fNumber==-1)
        {
            printf("WATERLOO: Can't find cell %d!\n",c);
            exit(0);
        }

        // Now that we know the file the cell is in, read a strip xCells long
        // Add to vCOM
        
        for(j=0;j<xCells;j++)
        {
            // Determines cell center position from x,y,z values
            pos[0] = -width[0]*pcTOcm/2.0 + DeltaX/2.0 + ((double)j + (double)CurCord[0])*DeltaX;
            pos[1] = -width[0]*pcTOcm/2.0 + DeltaX/2.0 + ((double)CurCord[1])*DeltaX;
            pos[2] = -width[0]*pcTOcm/2.0 + DeltaX/2.0 + ((double)CurCord[2])*DeltaX;
            
            cout << "Reading cell " << c+j << " from file " << fNumber << endl;
            ReadCell(InFiles[fNumber], curValues);
            for(i=0;i<8;i++) cout << curValues[i] << " ";
            cout << endl;
            
            // Increment Rho
            RhoT = RhoT + curValues[0];
            
            // Increment vCOM
            for(i=0;i<3;i++) vCOM[i] = vCOM[i] + curValues[0]*curValues[i+1];
            
            // Bookkeeping
            totMass = totMass + curValues[0]*pow(DeltaX,3);
            for(i=0;i<3;i++) totMom[i] = totMom[i] + curValues[0]*pow(DeltaX,3)*(curValues[i+1]);
            for(i=0;i<3;i++) angM[i] = angM[i] + curValues[0]*pow(DeltaX,3)*(pos[(i+1)%3]*curValues[1+(i+2)%3] - pos[(i+2)%3]*curValues[1+(i+1)%3]);
        }
        
        c = c+xCells;
        
    }
    
    
    // Before shift, the total mass and center of mass velocities
    cout << "Before shifting velocities, the total mass and momentums are:" << endl;
    cout << "Total mass (preshift) : " << totMass << endl;
    cout << "Momentums (preshift) : " ;
    for(i=0;i<3;i++) cout << totMom[i] << " " ;
    cout << endl;
    cout << "Momentum sum (preshift) : " << totMom[0]+totMom[1]+totMom[2] << endl;
    cout << "AngMomentums (preshift) : " ;
    for(i=0;i<3;i++) cout << angM[i] << " " ;
    cout << endl;
    
    
    totMass = 0;
    for(i=0;i<3;i++) totMom[i] = 0;
    for(i=0;i<3;i++) angM[i] = 0;
    
    // Now with these quantities, calculate COM velocities
    for(i=0;i<3;i++) vCOM[i] = vCOM[i]/RhoT;
    cout << "The center of mass velocities are: " ;
    for(i=0;i<3;i++) cout << vCOM[i] << " " ;
    cout << endl;
    
    
    // We will subtract this out from each cell that gets written
    c = 0;
    
    
    // Resets file positions
    for(n=0;n<numFiles;n++)
    {
        fseek ( InFiles[n] , 5*sizeof(int) + sizeof(double) , SEEK_SET );
    }
    
   
    
    while(1)
    {
        // Are we done?
        if(c==totCells) break;
        
        // Determine x,y,z cell values
        int CurCord[3];
        CurCord[0] = c % ref_lev[0];
        CurCord[1] = ((c-CurCord[0])/ref_lev[0]) % ref_lev[0];
        CurCord[2] = (c-CurCord[0]-CurCord[1]*ref_lev[0]) / pow(ref_lev[0],2);
        
        // Determines which cell range these values are in
        int PossibleFiles[numFiles];
        for(j=0;j<numFiles;j++) PossibleFiles[j] = j;
        
        // Narrows down
        for(j=0;j<3;j++)
        for(n=0;n<numFiles;n++)
        {
            if( CurCord[j] < CRanges[n][j][0] || CurCord[j] > CRanges[n][j][1] )
            {
                PossibleFiles[n] = -1;
            }
        }
        
        // Now only one entry should be different from -1
        int fNumber=-1;
        for(n=0;n<numFiles;n++)
        {
            if(PossibleFiles[n] > -1 && fNumber==-1 ) fNumber = n;
            else if(PossibleFiles[n] > -1 && fNumber!=-1)
            {
                printf("WATERLOO: Cell %d exists in multiple files!\n",c);
                exit(0);
            }
                
        }
        if(fNumber==-1)
        {
            printf("WATERLOO: Can't find cell %d!\n",c);
            exit(0);
        }
        cout << "cells " << c << " to " << c+xCells << " are in file " << fNumber << endl;
        
        
        // Now that we know the file the cell is in, read a strip xCells long
        // Also do some book keeping?
   
        for(j=0;j<xCells;j++)
        {
            // Determines cell center position from x,y,z values
            pos[0] = -width[0]*pcTOcm/2.0 + DeltaX/2.0 + ((double)j + (double)CurCord[0])*DeltaX;
            pos[1] = -width[0]*pcTOcm/2.0 + DeltaX/2.0 + ((double)CurCord[1])*DeltaX;
            pos[2] = -width[0]*pcTOcm/2.0 + DeltaX/2.0 + ((double)CurCord[2])*DeltaX;
            
            cout << "Reading and Writing cell " << c+j << " from file " << fNumber << endl;
            ReadCell(InFiles[fNumber], curValues);
            
            // Subtract out COM velocities
            for(i=0;i<3;i++) curValues[i+1] = curValues[i+1] - vCOM[i];
            
            for(i=0;i<8;i++) cout << curValues[i] << " ";
            cout << endl;
            
            WriteCell(OutFile, curValues);
            ReadCount[fNumber] = ReadCount[fNumber] + 1;
            
            totMass = totMass + curValues[0]*pow(DeltaX,3);
            for(i=0;i<3;i++) totMom[i] = totMom[i] + curValues[0]*pow(DeltaX,3)*(curValues[i+1]);
            for(i=0;i<3;i++) angM[i] = angM[i] + curValues[0]*pow(DeltaX,3)*(pos[(i+1)%3]*curValues[1+(i+2)%3] - pos[(i+2)%3]*curValues[1+(i+1)%3]);
        }
        
        c = c+xCells;
    }
    
    // Close the input files
    for(n=0;n<numFiles;n++)
    {
        fclose(InFiles[n]);
        printf("Closed file %d. Read from it %d times.\n",n,ReadCount[n]);
    }
    
 
    // Close to output file
    fclose(OutFile);
    printf("Closed the output file.\n");
    
    // Book keeping results
    cout << "Total mass (postshift)= " << totMass << endl;
    cout << "Momentums (postshift)= " << totMom[0] << " " << totMom[1] << " " << totMom[2] << endl;
    cout << "Momentum sum (postshift) : " << totMom[0]+totMom[1]+totMom[2] << endl;
    cout << "New center of mass velocity : "  << totMom[0]/totMass << " " << totMom[1]/totMass << " " << totMom[2]/totMass << endl;
    cout << "AngMomentums (postshift) : " ;
    for(i=0;i<3;i++) cout << angM[i] << " " ;
    cout << endl;
    
    // Tick tock, stop the clock
    timeMe = clock() - timeMe;
    printf("All done. Took %g wall-seconds.\n",((double)timeMe)/CLOCKS_PER_SEC);

    
    // Write out gadget_all file for debugging
    OutFile = fopen(out_file,"r");
    int temp;
    double tempd;
    fread(&temp,sizeof(int),1,OutFile);
    cout << temp << endl;
    fread(&temp,sizeof(int),1,OutFile);
    cout << temp << endl;
    fread(&tempd,sizeof(double),1,OutFile);
    cout << tempd << endl;
    for(j=0;j<CellsPP*totProc[0];j++)
    {
        double Array[8];
        fread(&Array[0],sizeof(double),8,OutFile);
        for(i=0;i<8;i++) cout << Array[i] << " " ;
        cout << endl;
    }
    fclose(OutFile);
    
    return 0;
}



