#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

using namespace std;



int main(int argc, char **argv)
{
    
    char InFileN[35] = "./gadget2orion_all";
    FILE *InFile;
    InFile = fopen(InFileN,"r");
    
    double Cell[8];
    int totProc,ref_lev;
    double width;
    int Cords[3];
    
    fread(&totProc,sizeof(int),1,InFile);
    fread(&ref_lev,sizeof(int),1,InFile);
    fread(&width,sizeof(double),1,InFile);
    //fread(&Cords[0],sizeof(int),3,InFile);
    
    int i,j;
    for(i=0;i<pow(ref_lev,3);i++)
    {
        fread(&Cell[0],sizeof(double),8,InFile);
        for(j=0;j<8;j++) cout << Cell[j] << " " ;
        cout << endl;
    }
    
    fclose(InFile);
    
    return 0;
}