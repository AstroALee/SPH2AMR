 /*
 * SPH2AMR.c
 *
 * Reads in Gadget2 data and projects the data onto a uniform grid suitable for
 * Orion2, which can be read in during initialization
 * Written by Athena Stacy and Aaron Lee, 2014
 *
 * Takes two inputs: the number of cells in each dimension (integer) and the width
 * of the box to carve out, in units of parsecs.
 *
 * To use, you need think about / set a few parameters before compiling:
 *
 *
 * preprocessor BREAD (line ~38ish): 0 or 1, whether you do not or do want to use
 *          the bread model for distributing data to the processors.
 *          0 : If no, then each processor allocates memory for the entire grid,
 *          and loops over a fraction of the particles. Uses a lot of memory for
 *          large grids.
 *          1 : Each processor gets a fraction of the grids and loops over all
 *          the particles.
 *
 * Filename to read in (line ~216ish): Basename of file, the number of files, etc.
 * Filename to print out (line ~351ish): Filename to print out...
 *
 */

#include "sph2amr_final.h"


/* Here we load a snapshot file. It can be distributed onto several files (for files>1). The particles are brought back into the order implied by their
 ID's. A unit conversion routine is called to do unit conversion, and to
 evaluate the gas temperature.
*/
int main(int argc, char **argv)
{
    // Ready... Set... GO!
    clock_t timeMe;
    timeMe = clock();

    // counters
    int i,j,n;

    // Initialize the MPI Universe
    int npes, myrank, ierr=0;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &npes); // number processors
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank); // labels for each processor

    /* A few important input values are needed from the command line
     *
     * ref_lev : An integer giving the number of grid cells in each dimension
     *           for the resulting box
     * width   : A double given the size of each box side, in units of parsecs
     */


    if(argc < 4) // Double check on inputs
    {
        printf("WATERLOO: 3 inputs needed: (int)ref_lev, (double)width [unit: pc], and int(restart)\nYou gave %d.\n",argc-1);
        exit(-1);
    }


    ref_lev = atoi(argv[1]);
    width = atof(argv[2]);
    RestartMe = atoi(argv[3]);
    ref_lev_doub = ((double)ref_lev);



    /* I/O Information */

    char path_in[200], path_out[200], basename[200], input_fname[200], input_fname2[200];

    //If we are taking a particular snapshot number, identify it here
    int snapshot_number = 2;

    //Mulitiple files to load Gadget data from?
    int files=1;

    sprintf(path_in, pathname_in);
    sprintf(path_out, pathname_out);
    sprintf(basename, "bin_HR9");
    sprintf(input_fname, "%s/%s_%03d", path_in, basename, snapshot_number);
#if(READB)
    sprintf(input_fname2, "%s_bfield.dat", path_in);
#endif

    /* Output file */
    char  output_fname[200] ;
    sprintf(output_fname, "%s/gadget2orion", path_out);




    /* Prints out some information */
    if(myrank==0) {
    printf("\n\nWelcome to the Gadget to Orion2 Projection Code!\n\n");
    printf("Here are the inputs:\n");
    printf("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n");
    printf("You grid will be \t\t%d cells.\n",ref_lev);
    printf("The width will be \t\t%g pc.\n",width);
    printf("Are you restarting? \t\t%d\n",RestartMe);
    printf("Your file input name: \t\t%s\n",input_fname);
#if(READB)
    printf("Your B-field input name: \t\t%s\n",input_fname2);
#endif
    printf("Your output file name: \t\t%s\n",output_fname);
    printf("MPI--Num of processes: \t\t%d\n",npes);
    printf("\n\n");
    }




    /* Checks and Balances */

    // Are you trying to use an old method yet restart? What's happening?!
    if(!SIMPBREAD && (RestartMe==1) )
    {
        printf("WATERLOO: Cannot restart if you are not using the Simpson Bread method.\n");
        exit(-1);
    }

    // Does the number of processors divide evenly into the z-coordinate?
    if(BREAD && (ref_lev%npes != 0) && !SIMPBREAD)
    {
        printf("WATERLOO: For the bread model, we need the value of ref_lev to be divided evenly by the number of processors \n");
        exit(-1);
    }

    if(SIMPBREAD)
    {
        BreakUpDomain(Cords,npes);
        int temp = Cords[0];
        Cords[0] = Cords[2];
        Cords[2] = temp;

        for(j=0;j<3;j++) if(ref_lev%Cords[j] != 0)
        {
            printf("WATERLOO: Need the distribution of the number of processors and ref_lev to play nice with each other. Ideally, you'll have both be powers of 2... \n");
            exit(-1);
        }

        // Tries to load in the z-dimension
        while(1)
        {
            if( (Cords[0]>1) && (ref_lev % (2*Cords[2]) == 0) && (Cords[0] % 2) == 0)
            {
                Cords[2] = Cords[2]*2;
                Cords[0] = Cords[0]/2;
            }
            else
            {
                break;
            }

        }

        while(1)
        {
            if( (Cords[1]>1) && (ref_lev % (2*Cords[2]) == 0) && (Cords[1] % 2) == 0)
            {
                Cords[2] = Cords[2]*2;
                Cords[1] = Cords[1]/2;
            }
            else
            {
                break;
            }

        }

        // Final effort to load the y-dimension
        while(1)
        {
            if( (Cords[0]>1) && (ref_lev % (2*Cords[1]) == 0) && (Cords[0] % 2) == 0)
            {
                Cords[1] = Cords[1]*2;
                Cords[0] = Cords[0]/2;
            }
            else
            {
                break;
            }

        }


        if(myrank==0) if(DEBUGGING || CHATTY) {
        printf("Divided the Grid: ");
        for(i=0;i<3;i++) printf("Cord[%d]=%d ",i,Cords[i]);
        printf("\n\n");
        }

    }


    /* Reads in Gadget Particle Data */

    // Each processor will have access to this data
    double DelCoord[3] = {0,0,0};
    double boxsize = 140.0; // I don't think the 140 value matters
    Ngas = read_snapshot(input_fname, files, DelCoord, boxsize,myrank);


    // If dealing with B-field from analytic calculations, read that in here.
    FILE *infile;
#if(READB)
    infile = fopen(input_fname2, "r");
    for(n=0;n<Ngas;n++)
        fread(&P[n].Bfield[0], sizeof(double), 3, infile);
    for(n=0;n<Ngas;n++)
        fread(&P[n].nh_test, sizeof(double), 1, infile);
    fclose(infile);
#endif



    // Pesky units
    unit_conversion();


    /* With reading done, some more constants and conversions */

    // To convert from co-moving to physical
    double CtoP = 1.e3*Time/(hubble_param);

    // Define physical smoothing length
    for(n=0;n<Ngas;n++) P[n].hsm_phys = CtoP*P[n].hsm;

    // Calculate total mass in entire simulation
    double Part_Mtot = 0.0;
    for(n=0;n<Ngas;n++) Part_Mtot = Part_Mtot + P[n].Mass;
    if(myrank==0 && CHATTY) printf("MASSTEST: The mass in the entire simulation box is = %g grams (%g solar masses)\n",mass_conv*Part_Mtot,mass_conv*Part_Mtot/solarMass);



    /* Calculate comoving COM position of entire domain  */
    double pCOM[3] = {0,0,0};
    for(n=0;n<Ngas;n++)
    {
        pCOM[0] = pCOM[0] + P[n].Pos[0]*P[n].Mass;
        pCOM[1] = pCOM[1] + P[n].Pos[1]*P[n].Mass;
        pCOM[2] = pCOM[2] + P[n].Pos[2]*P[n].Mass;
    }
    pCOM[0]=pCOM[0]/Part_Mtot; // Comoving units
    pCOM[1]=pCOM[1]/Part_Mtot; // mass unit irrelevant; divided out
    pCOM[2]=pCOM[2]/Part_Mtot;


    /* Calculate position of Density maximum */
    double pDMax[3];
    int nhmax = 0;
    for(n=0;n<Ngas;n++)
    {
        int nh = P[n].Density; //used to be P[n].nh
        if(nh > nhmax && P[n].sink > -1)
        {
            nhmax = nh;
            pDMax[0] = P[n].Pos[0];  //del's based upon location of maximum density
            pDMax[1] = P[n].Pos[1];
            pDMax[2] = P[n].Pos[2]; //  keep in comoving coordinates
        }
    }






    if(myrank==0 && (DEBUGGING || CHATTY)) printf("physical location of COM and densest location (physical units): (%g,%g,%g) , (%g,%g,%g)\n",CtoP*pCOM[0],CtoP*pCOM[1],CtoP*pCOM[2],CtoP*pDMax[0],CtoP*pDMax[1],CtoP*pDMax[2]);




    // Calculate the center of mass position & velocity for particles
    // in box of interest, centered about densest point
    //double InterestMtot = 0.0; // only mass in box we're projecting for
    int Interestcount = 0;
    double pTot[3] = {0,0,0};
    double vCOM[3] = {0,0,0};
    double pCOMlocal[3] = {0,0,0};
    for(n=0;n < Ngas; n++)
    {

        // 1D distances, takes into account smoothing lengths
        double pLefts[3]={0,0,0};
        double pRights[3]={0,0,0};
        for(i=0;i<3;i++)
        {
            pLefts[i] = CtoP*(P[n].Pos[i]-pDMax[i])+hfac*P[n].hsm_phys;
            pRights[i] = CtoP*(P[n].Pos[i]-pDMax[i])-hfac*P[n].hsm_phys;
        }
        int DoWeCare = 1;

        for(i=0;i<3;i++) if( pLefts[i] < -width/2.0) DoWeCare=0;
        for(i=0;i<3;i++) if( pRights[i] > width/2.0) DoWeCare=0;

        // Uncomment below lines if inside we want to only care about particles
        // 100% inside the box
        //for(i=0;i<3;i++) if( pRights[i] < -width/2.0 ) DoWeCare=0;
        //for(i=0;i<3;i++) if( pLefts[i]  >  width/2.0 ) DoWeCare=0;

        if(DoWeCare) // overlaps with box we want
        {
            for(j=0;j<3;j++) vCOM[j] = vCOM[j] + P[n].Vel[j]*P[n].Mass;
            for(j=0;j<3;j++) pTot[j] = pTot[j] + P[n].Vel[j]*P[n].Mass;
            for(j=0;j<3;j++) pCOMlocal[j] = pCOMlocal[j] + (P[n].Pos[j]-pDMax[i])*P[n].Mass;
            InterestMtot = InterestMtot + P[n].Mass;
            Interestcount = Interestcount + 1;
        }

    }
    for(j=0;j<3;j++) vCOM[j] = vCOM[j]/InterestMtot; // Kept in co-moving units
    for(j=0;j<3;j++) pCOMlocal[j] = pCOMlocal[j]/InterestMtot;

    double vel_conv = 1.e5*pow(Time,0.5);
    for(j=0;j<3;j++) pTot[j] = pTot[j]*mass_conv*vel_conv;

    InterestMtot = mass_conv*InterestMtot;
    if(myrank==0 && (CHATTY || DEBUGGING))
    {
        printf("TESTME : INSIDE SIM: Total gadget particles %d, m_gadget = %10.10g grams (%10.10g solar masses)\n",Interestcount,InterestMtot,InterestMtot/solarMass);

        for(n=0;n<3;n++) printf("local phy COM[%d] = %lg , ",n,CtoP*pCOMlocal[n]);
        printf("\n");
        for(n=0;n<3;n++) printf("local phy vCOM[%d] = %lg , ",n,vel_conv*vCOM[n]);
        printf("\n");
        printf("Total phy momentum:");
        for(n=0;n<3;n++) printf(" Mom[%d] = %g , ",n,pTot[n]);
        printf("\n");
    }




#if(SIMPBREAD)
    char restartfilename[200];
    sprintf(restartfilename,"%s/restart.log%d",path_out,myrank);
#endif

    // Creates file(S) if we're not restarting
    if(!SIMPBREAD && myrank==0) //  (overwrites if already exists!)
    {
        FILE *outfile;
        outfile=fopen(output_fname,"w");
        fclose(outfile);
    }
    else if(SIMPBREAD && RestartMe==0)
    {
        // Each processor has its own file, create them here.
        sprintf(output_fname, "%s/gadget2orion_%d", path_out,myrank);
        FILE *outfile;
        outfile=fopen(output_fname,"w");
        fclose(outfile);

        // Also creates restart log if necessary
#if(RESTARTFILE)
        FILE *restartfile;
        restartfile = fopen(restartfilename,"w");
        int initrestart = 0;
        fwrite(&initrestart,sizeof(int),1,restartfile);
        fclose(restartfile);
#endif

    }
    else if(SIMPBREAD && RestartMe==1)
    {
        // Files already exist, don't overwrite them! But do they exist...
        sprintf(output_fname, "%s/gadget2orion_%d", path_out,myrank);
        std::ifstream myfile(output_fname);
        if(!myfile)
        {
            printf("WATERLOO: Trying to restart from file that doesn't exist!\n");
            printf("Cannot find file %s\n",output_fname);
            exit(-1);
        }

        // What about the restart log?
        std::ifstream myrestartfile(restartfilename);
        if(!myrestartfile)
        {
            printf("WATERLOO: Trying to restart, but cannot find a restart.log!\n");
            printf("Cannot find file %s\n",restartfilename);
            exit(-1);
        }

    }
    else
    {
        printf("WATERLOO: Not sure what to do regarding output!\n");
        exit(-1);
    }


    // About to start the real work.
    /* What are you centering on: pDMax ?
        or pCOM (of the entire simulation)?
        or pCOMlocal?
       What are you shifting velocities on: vCOM?
     */



    int Ngas2;
#if SIMPBREAD //
    Ngas2 = Projection_SimpBread(output_fname, restartfilename, pDMax, vCOM, npes, myrank);
#elif BREAD // Less memory intensive bread model
    Ngas2 = Projection_Bread(output_fname, pDMax, vCOM, npes, myrank);
#else // The original memory intensive method.
    Ngas2 = write_snapshot(input_fname, files, output_fname, delx, dely, delz, vCOM, Ngas, myrank,ref_lev,width);
#endif



    //if(myrank==0 && (DEBUGGING || CHATTY)) printf("Particle Count = %d\n",ParticleCounts);

    // Tick tock, stop the clock
    timeMe = clock() - timeMe;
    printf("Processor %d took %g wall-seconds on %d processors (dim,width = %d,%g pc).\n",myrank,((double)timeMe)/CLOCKS_PER_SEC,npes,ref_lev,width);


    //double timeAvg = 0;
    //double timeTook = ((double)timeMe)/CLOCKS_PER_SEC;
    //MPI_Allreduce(&timeTook, &timeAvg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //timeAvg = timeAvg / ( (double) npes) ;
    //printf("Average wall-second time amongst all processors : %g \n",timeAvg);

    // Clean up
    ierr=MPI_Finalize();


    return 0;
} // End of main()




/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 Simpson Projection Code
 =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */

// New attempt at writing snapshot writing so that it is less memory intensive and arranges the data smartly
int Projection_SimpBread(char *outname, char *restartfilename, double pCenter[], double vREL[], int totProc, int myrank)
{

    /* Inputs:
     outname: name of output file
     pCENTER: co-moving coordinates of location to center on
     vREL: co-moving relative velocity to subtract out
     totProc: total number of processors
     myrank: Rank of working processor
     */


    //Iterators
    int i,j,k,v,n;


    // Some conversions that will be useful later
    double vel_conv = 1.e5*pow(Time,0.5);
    double kernel_conv =  pow(pcTOcm*1.e3*Time/(hubble_param),-3);
    double CtoP = 1.e3*Time/(hubble_param);


    // If restarting, get the current cell numbers
    FILE* restartfile;
    int CurCellNum = 0;
    if(RestartMe)
    {
        restartfile = fopen(restartfilename,"r");
        fread(&CurCellNum,sizeof(int),1,restartfile);
        fclose(restartfile);
    }

    if(RestartMe) printf("We are engaging with cell number %d (on proc %d)\n",CurCellNum,myrank);


    /* Each processor has a unique value for myrank (out of totProc). Use this
     * and the number of processors to determine which grid cells it searches
     * over
     */

    // Using width, break each dimension into rectangles
    double Rect[3] = {width/( (double)Cords[0] ), width/( (double)Cords[1] ), width/( (double)Cords[2] ) };

    // Grid cell length
    double DeltaX = width/ref_lev_doub;
    double DeltaXh = DeltaX/2.0;

    // Cell range for each processor (different for each processor)
    int xCRange[2],yCRange[2],zCRange[2];
    int xCells,yCells,zCells;
    xCells=ref_lev/Cords[0];
    yCells=ref_lev/Cords[1];
    zCells=ref_lev/Cords[2]; // These are same for each proc

    for(i=0;i<totProc;i++)
    {
        int xBlock = i % Cords[0];
        int yBlock = ((i-xBlock)/Cords[0]) % Cords[1];
        int zBlock = (i-xBlock-yBlock*Cords[0])/Cords[0]/Cords[1];

        if(myrank==i)
        {
            xCRange[0] = xBlock*xCells;
            xCRange[1] = (1+xBlock)*xCells-1;
            yCRange[0] = yBlock*yCells;
            yCRange[1] = (1+yBlock)*yCells-1;
            zCRange[0] = zBlock*zCells;
            zCRange[1] = (1+zBlock)*zCells-1;
            if(DEBUGGING) printf("rank %d : %d,%d   %d,%d   %d,%d\n",myrank,xCRange[0],xCRange[1],yCRange[0],yCRange[1],zCRange[0],zCRange[1]);
        }
    }
    //zCRange[0] = myrank*zCells; // These two lines are different for each proc
    //zCRange[1] = (1+myrank)*zCells-1;


    // Number of cells stored on each processor
    int CellsPP = xCells*yCells*zCells;


    double perDone = ((double)CurCellNum)/((double)CellsPP);
    int PerInc = round( 0.05*((double)CellsPP) );
    int NextPer = PerInc;
    while(1)
    {
        if(NextPer < CurCellNum)
        {
            NextPer = NextPer+PerInc;
        }
        else
        {
            break;
        }
    }

    /*
    // Physical Z distance range for each processor (different for each proc)
    // Centers on middle of cell
    double zRange[2];
    zRange[0] = ((double) myrank)*Rect[2] - width/2.0 + DeltaXh;
    zRange[1] = zRange[0] + Rect[2] - DeltaX;
    //xRange[0] = -width/2.0 + DeltaXh;
    //xRange[1] = xRange[0] + RecX - DeltaX;


    if(CHATTY) printf("myrank = %d, cell ranges for z: (%d,%d)\n",myrank, zCRange[0],zCRange[1]);
    if(CHATTY) printf("myrank = %d, distance ranges for z: (%g,%g)\n",myrank,zRange[0],zRange[1]);
    */


    /*  Here we subtract out velocities and shift particle positions so our
     center location is at the "center" of the Orion2 box
     ( "center" = (-1/2,-1/2,-1/2)*DeltaX )
     */
    double gMassTot=0.0;
    double gMomTot[3] = {0,0,0};
    int gNcount=0;


    // Instead going to recalculate the vCOM and pCenter values and compare
    double pCenterLocal[3] = { 0,0,0};

    int dMaxIndex = 0; double dMaxVal = 0.0;

    // Find max location of entire box (will be inside box of interest, presumably)
    for(n=0;n<Ngas;n++) if(P[n].Density > dMaxVal && P[n].sink > -1)
        {
            dMaxIndex = n;
            dMaxVal = P[n].Density;
        }

    // Subtract out position, make physical units
    for(n=0;n<Ngas;n++)
    {
        P[n].disx = CtoP*(P[n].Pos[0] - P[dMaxIndex].Pos[0]);
        P[n].disy = CtoP*(P[n].Pos[1] - P[dMaxIndex].Pos[1]);
        P[n].disz = CtoP*(P[n].Pos[2] - P[dMaxIndex].Pos[2]);
    }

    // Check on densest location
    int dMaxIdxTest = 0; double dMaxValTest = 0.0;
    for(n=0;n<Ngas;n++) if(P[n].Density > dMaxValTest && P[n].sink > -1)
        {
            dMaxIdxTest = n;
            dMaxValTest = P[n].Density;
        }
    if(DEBUGGING) printf("Postshift check on densest particle location : %g %g %g \n",P[dMaxIdxTest].disx,P[dMaxIdxTest].disy,P[dMaxIdxTest].disz);

    // Calculate angular momentum using only paricles of interest
    double angMLocal[3] = {0,0,0};
    for(n=0;n<Ngas;n++)
    {
        //Determines if particles are inside box of interest
        int DoWeCare = 1;

        // Is there any overlap?
        if( P[n].disx - hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        if( P[n].disy - hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        if( P[n].disz - hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        if( P[n].disx + hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        if( P[n].disy + hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        if( P[n].disz + hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;

        // Uncomment these lines if we only care about totally enclosed particles
        //if( P[n].disx - hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        //if( P[n].disy - hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        //if( P[n].disz - hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        //if( P[n].disx + hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        //if( P[n].disy + hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        //if( P[n].disz + hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;

        if(DoWeCare)
        {
            angMLocal[0] = angMLocal[0] + P[n].Mass*(P[n].disy*P[n].Vel[2] - P[n].disz*P[n].Vel[1]);
            angMLocal[1] = angMLocal[1] + P[n].Mass*(P[n].disz*P[n].Vel[0] - P[n].disx*P[n].Vel[2]);
            angMLocal[2] = angMLocal[2] + P[n].Mass*(P[n].disx*P[n].Vel[1] - P[n].disy*P[n].Vel[0]);
        }
    }

    if(DEBUGGING) printf("TESTME: Total AngMomentum is %g %g %g\n",angMLocal[0]*mass_conv*vel_conv,angMLocal[1]*mass_conv*vel_conv,angMLocal[2]*mass_conv*vel_conv);

    // Calculate Vcom using only particles in box of interest
    double vCOMLocal[3] = {0,0,0};
    double MassOfInterest = 0;
    for(n=0;n<Ngas;n++)
    {
        //Determines if particles are inside box of interest
        int DoWeCare = 1;

        // Is there any overlap?
        if( P[n].disx - hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        if( P[n].disy - hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        if( P[n].disz - hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        if( P[n].disx + hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        if( P[n].disy + hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        if( P[n].disz + hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;

        // Uncomment these lines if we only care about totally enclosed particles
        //if( P[n].disx - hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        //if( P[n].disy - hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        //if( P[n].disz - hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        //if( P[n].disx + hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        //if( P[n].disy + hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        //if( P[n].disz + hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;

        if(DoWeCare)
        {
            MassOfInterest = MassOfInterest + P[n].Mass;
            for(i=0;i<3;i++) vCOMLocal[i] = vCOMLocal[i] + P[n].Vel[i]*P[n].Mass;
        }

    }

    if(DEBUGGING)
    {
    printf("TESTING: Center is located %g %g %g\n",P[dMaxIndex].Pos[0],P[dMaxIndex].Pos[1],P[dMaxIndex].Pos[2]);
    printf("TESTING: Previous Center was at %g %g %g\n",pCenter[0],pCenter[1],pCenter[2]);
    printf("TESTME: Total Momentum is %g %g %g\n",vCOMLocal[0]*mass_conv*vel_conv,vCOMLocal[1]*mass_conv*vel_conv,vCOMLocal[2]*mass_conv*vel_conv);
    printf("TESTING: Local vCOM is %g %g %g\n",vCOMLocal[0]*vel_conv/MassOfInterest,vCOMLocal[1]*vel_conv/MassOfInterest,vCOMLocal[2]*vel_conv/MassOfInterest);
    }



    // Calculate angular momentum using only paricles of interest with velocity shift
    for(i=0;i<3;i++) angMLocal[i] = 0.0;
    for(n=0;n<Ngas;n++)
    {
        //Determines if particles are inside box of interest
        int DoWeCare = 1;

        // Is there any overlap?
        if( P[n].disx - hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        if( P[n].disy - hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        if( P[n].disz - hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        if( P[n].disx + hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        if( P[n].disy + hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        if( P[n].disz + hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;

        // Uncomment these lines if we only care about totally enclosed particles
        //if( P[n].disx - hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        //if( P[n].disy - hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        //if( P[n].disz - hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        //if( P[n].disx + hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        //if( P[n].disy + hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        //if( P[n].disz + hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;

        if(DoWeCare)
        {
            angMLocal[0] = angMLocal[0] + P[n].Mass*(P[n].disy*(P[n].Vel[2]-vCOMLocal[2]/MassOfInterest) - P[n].disz*(P[n].Vel[1]-vCOMLocal[1]/MassOfInterest));
            angMLocal[1] = angMLocal[1] + P[n].Mass*(P[n].disz*(P[n].Vel[0]-vCOMLocal[0]/MassOfInterest) - P[n].disx*(P[n].Vel[2]-vCOMLocal[2]/MassOfInterest));
            angMLocal[2] = angMLocal[2] + P[n].Mass*(P[n].disx*(P[n].Vel[1]-vCOMLocal[1]/MassOfInterest) - P[n].disy*(P[n].Vel[0]-vCOMLocal[0]/MassOfInterest));
        }
    }

    if(DEBUGGING) printf("TESTING: Postshift Total AngMomentum is %g %g %g\n",angMLocal[0]*mass_conv*vel_conv,angMLocal[1]*mass_conv*vel_conv,angMLocal[2]*mass_conv*vel_conv);





    // Now shift velocities (all in comoving)
    //for(n=0;n<Ngas;n++)for(i=0;i<3;i++) P[n].Vel[i] = P[n].Vel[i] - vCOMLocal[i]/MassOfInterest;


    // Now calculate momentum of particles of interest.
    for(i=0;i<3;i++) vCOMLocal[i] = 0;
    double CheckSum1 = 0;
    for(n=0;n<Ngas;n++)
    {
        //Determines if particles are inside box of interest
        int DoWeCare = 1;

        // Is there any overlap?
        if( P[n].disx - hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        if( P[n].disy - hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        if( P[n].disz - hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        if( P[n].disx + hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        if( P[n].disy + hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        if( P[n].disz + hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;

        // Uncomment these if we only care about totally enclosed particles
        //if( P[n].disx - hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        //if( P[n].disy - hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        //if( P[n].disz - hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
        //if( P[n].disx + hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        //if( P[n].disy + hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
        //if( P[n].disz + hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;

        if(DoWeCare)
        {
            for(i=0;i<3;i++) vCOMLocal[i] = vCOMLocal[i] + P[n].Vel[i]*P[n].Mass;
            CheckSum1 = CheckSum1 + n;
        }

    }

    if(DEBUGGING) printf("TESTING: Total Pre-Shift Momentum is %g %g %g\n",vCOMLocal[0]*mass_conv*vel_conv,vCOMLocal[1]*mass_conv*vel_conv,vCOMLocal[2]*mass_conv*vel_conv);




    //return 0;




    // Given the way the data will be layed out now, we'll need some header
    // information so we can reconstruct it correctly in Orion2.
    FILE *outfile;
    if(RestartMe==0)
    {
        outfile=fopen(outname,"w");
        fwrite(&totProc, sizeof(int), 1, outfile);
        fwrite(&ref_lev, sizeof(int), 1, outfile);
        fwrite(&width, sizeof(double), 1, outfile);
        fwrite(&Cords[0], sizeof(int), 3, outfile);
        fclose(outfile);
        printf("Header information for file %s (proc %d) written.\n",outname,myrank);
    }
    else if(RestartMe==1)
    {
        // Do nothing to the output file.
    }
    else
    {
        // Should never get to this
        printf("WATERLOO: Not sure what is happening with output!\n");
    }



    /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
                    Here is where the big loops start
       =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */


    // For each cell
    if(myrank==0) printf("Beginning grand loop over cells (proc %d).\n",myrank);

    double StoredData[OUTFREQ][varnum];
    int timeToPrint = 0;

    double totProjMass = 0;
    double totProjMom[3] = {0,0,0};
    double CheckSum2 = 0;
    while(1)
    {
        // Are we done?
        if(CurCellNum==CellsPP) break;


        // Data for current cell
        double curData[varnum];
        for(i=0;i<varnum;i++) curData[i]=0.0;

        // Clear Touch Data for particles
        for(i=0;i<Ngas;i++) P[i].TouchMe = 0;
        for(i=0;i<Ngas;i++) P[i].TouchRho = 0.0;

        /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
            Physical Extent of current cell
         =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */

        // given myrank, we know the block
        int xBlock = myrank % Cords[0];
        int yBlock = ((myrank-xBlock)/Cords[0]) % Cords[1];
        int zBlock = (myrank-xBlock-yBlock*Cords[0])/Cords[0]/Cords[1];

        // given CurCellNum, we can derive back x,y,z cell numbers
        int curX = CurCellNum % xCells;
        int curY = ((CurCellNum-curX)/xCells) % yCells;
        int curZ = (CurCellNum-curX-curY*xCells)/xCells/yCells;

        // and from these, physical extents of the cell boundaries
        double c_Ranges[3][2];
        c_Ranges[0][0] = ((double) xBlock)*Rect[0] -width/2.0 + ((double)curX)*DeltaX ;
        c_Ranges[0][1] = ((double) xBlock)*Rect[0] -width/2.0 + ((double)curX)*DeltaX+DeltaX;


        c_Ranges[1][0] = ((double) yBlock)*Rect[1] -width/2.0 + ((double)curY)*DeltaX;
        c_Ranges[1][1] = ((double) yBlock)*Rect[1] -width/2.0 + ((double)curY)*DeltaX+DeltaX;

        c_Ranges[2][0] = ((double) zBlock)*Rect[2] - width/2.0 + ((double)curZ)*DeltaX ;
        c_Ranges[2][1] = ((double) zBlock)*Rect[2] - width/2.0 + ((double)curZ)*DeltaX + DeltaX;

        // What's going on?
        if(DEBUGGING && CHATTY)
        {
	        printf("Proc %d at cell %d of %d\n",myrank,CurCellNum+1,CellsPP);
            //printf("Processor %d is dealing with cell %d of %d.\n",myrank,CurCellNum+1, CellsPP);
            //for(n=0;n<3;n++) printf("\tRanges %d = %g %g\n",n,c_Ranges[n][0],c_Ranges[n][1]);
        }
        //exit(0);

        // Loop over each particle
        for(n=0;n<Ngas;n++)
        {
            //if(n != 942485 && n != 1075919 && n != 1075951 && n != 1075955 && n != 1036372) continue;


            /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
                Determine if particle touches the cell
               =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */

            // Physical extent of each particle (includes smoothing length)
            double p_xRange[2] = {P[n].disx - hfac*P[n].hsm_phys , P[n].disx + hfac*P[n].hsm_phys };
            double p_yRange[2] = {P[n].disy - hfac*P[n].hsm_phys , P[n].disy + hfac*P[n].hsm_phys };
            double p_zRange[2] = {P[n].disz - hfac*P[n].hsm_phys , P[n].disz + hfac*P[n].hsm_phys };

            // Determines if this particle touches cell
            int DoWeCare = 1;

            // Is there any overlap?
            if( P[n].disx - hfac*P[n].hsm_phys > c_Ranges[0][1])  DoWeCare=0;
            if( P[n].disy - hfac*P[n].hsm_phys > c_Ranges[1][1])  DoWeCare=0;
            if( P[n].disz - hfac*P[n].hsm_phys > c_Ranges[2][1])  DoWeCare=0;
            if( P[n].disx + hfac*P[n].hsm_phys < c_Ranges[0][0]) DoWeCare=0;
            if( P[n].disy + hfac*P[n].hsm_phys < c_Ranges[1][0]) DoWeCare=0;
            if( P[n].disz + hfac*P[n].hsm_phys < c_Ranges[2][0]) DoWeCare=0;

            // Uncomment these if we only care about totally enclosed particles (in that cell)
            /*
            if( P[n].disx - hfac*P[n].hsm_phys < c_Ranges[0][0]) DoWeCare=0;
            if( P[n].disy - hfac*P[n].hsm_phys < c_Ranges[1][0]) DoWeCare=0;
            if( P[n].disz - hfac*P[n].hsm_phys < c_Ranges[2][0]) DoWeCare=0;
            if( P[n].disx + hfac*P[n].hsm_phys > c_Ranges[0][1])  DoWeCare=0;
            if( P[n].disy + hfac*P[n].hsm_phys > c_Ranges[1][1])  DoWeCare=0;
            if( P[n].disz + hfac*P[n].hsm_phys > c_Ranges[2][1])  DoWeCare=0;
            */

            // Uncomment these lines if we only care about particles inside the
            // entire simulation domain (but not necessarily totally inside a given
            // cell.
            //if( P[n].disx - hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
            //if( P[n].disy - hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
            //if( P[n].disz - hfac*P[n].hsm_phys < -width/2.0) DoWeCare=0;
            //if( P[n].disx + hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
            //if( P[n].disy + hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;
            //if( P[n].disz + hfac*P[n].hsm_phys > width/2.0)  DoWeCare=0;

            // If out of range for one coordinate, 'continue' back to top of loop
            if(!DoWeCare) continue;



            /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
              Particle touches cell, calculate volume fraction
             =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
            P[n].TouchMe = 1;

            // If so, evaluate fraction of particle that gets projected
            double frac = 1.0;
            //double CtoPvol = pow(CtoP,-3);


            //In the case where the particle is entirely inside the cell, project everything.
            int AllInside = 1;

            if( p_xRange[0] < c_Ranges[0][0] ) AllInside = 0;
            if( p_xRange[1] > c_Ranges[0][1] ) AllInside = 0;
            if( p_yRange[0] < c_Ranges[1][0] ) AllInside = 0;
            if( p_yRange[1] > c_Ranges[1][1] ) AllInside = 0;
            if( p_zRange[0] < c_Ranges[2][0] ) AllInside = 0;
            if( p_zRange[1] > c_Ranges[2][1] ) AllInside = 0;

            if(!AllInside)
            {
                //printf("%g %g %g\n",kernel_conv,kernel_conv * pow(pcTOcm,3),pow(CtoP,-3));
                //Quadrature_Simpson(frac,P[n],c_Ranges,DeltaX);
                Quadrature_MCarlo(frac,P[n],c_Ranges,DeltaX,myrank);
                //Quadrature_Centering(frac,P[n],c_Ranges,DeltaX);

            }
            //printf("frac = %g\n",frac);


            /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
             Fraction calculated, determine density
             =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */

            // Density (comoving)
            P[n].TouchRho = (P[n].Mass)*frac/pow(pcTOcm*DeltaX,3);
            curData[0] = curData[0] + P[n].TouchRho;

            // Book Keeping (comoving)
            totProjMass = totProjMass + (P[n].Mass)*frac;

        }

        //  This cell now has a total Rho (curData[0]). Use this to compute the other variables
        //double AnotherBook[3] = {0,0,0};
        //double AnotherBook2[3] = {0,0,0};

        // If nothing gets projected, we don't want to divide by zero!
        if(curData[0] <= 0.0)
        {
            for(n=0;n<8;n++) curData[n] = 0.0;

        }
        else
        {
            for(n=0;n<Ngas;n++)
            {
                if(P[n].TouchMe==1)
                {
                    double MassWeightedFrac = P[n].TouchRho / curData[0];

                    // Les Velocities
                    curData[1] = curData[1] + (P[n].Vel[0])*MassWeightedFrac;
                    curData[2] = curData[2] + (P[n].Vel[1])*MassWeightedFrac;
                    curData[3] = curData[3] + (P[n].Vel[2])*MassWeightedFrac;
                    // Energy Density
                    curData[4] = curData[4] + P[n].U*MassWeightedFrac;
                    // Chemical Species
                    curData[5] = curData[5] + P[n].H2I*MassWeightedFrac;
                    curData[6] = curData[6] + P[n].HDI*MassWeightedFrac;
                    curData[7] = curData[7] + P[n].HII*MassWeightedFrac;
                    // B-field
#if(READB)
                    curData[8] = curData[8] + P[n].Bfield[0]*MassWeightedFrac;
                    curData[9] = curData[9] + P[n].Bfield[1]*MassWeightedFrac;
                    curData[10] = curData[10] + P[n].Bfield[2]*MassWeightedFrac;
#endif


                    // Book Keeping
                    CheckSum2 = CheckSum2 + n;

                    //for(i=0;i<3;i++) AnotherBook[i] = AnotherBook[i] + P[n].Mass*P[n].Vel[i];
                    //for(i=0;i<3;i++) AnotherBook2[i] = AnotherBook2[i] + P[n].Mass*P[n].Vel[i]*mass_conv*vel_conv;

                }
            }

        }



        // All done, convert to physical units
        curData[0] = curData[0]*mass_conv;
        curData[1] = curData[1]*vel_conv;
        curData[2] = curData[2]*vel_conv;
        curData[3] = curData[3]*vel_conv;


        //printf("CHECKSUMS: %d %d (diff: %g)\n",CheckSum1,CheckSum2,CheckSum1-CheckSum2);
        //printf("ANOTHERBOOK: %g %g %g\n",AnotherBook[0]*mass_conv*vel_conv,AnotherBook[1]*mass_conv*vel_conv,AnotherBook[2]*mass_conv*vel_conv);
        //printf("ANOTHERBOOK2: %g %g %g\n",AnotherBook2[0],AnotherBook2[1],AnotherBook2[2]);

        // Book keeping (cgs units)
        //totProjMass = totProjMass + curData[0]*pow(pcTOcm*DeltaX,3);
        for(i=0;i<3;i++) totProjMom[i] = totProjMom[i] + curData[0]*pow(pcTOcm*DeltaX,3)*curData[i+1];


        // This particular cell now has its complete set of data.
        for(i=0;i<varnum;i++) StoredData[timeToPrint][i] = curData[i];
        timeToPrint = timeToPrint + 1;

        //Do we print yet?
        if(timeToPrint == OUTFREQ || CurCellNum == CellsPP-1)
        {
            printf("Printing on Processor %d (%d: %d to %d, %d)!\n",myrank,timeToPrint,CurCellNum-timeToPrint+1,CurCellNum,CellsPP);
            outfile=fopen(outname,"a"); // Appending
            for(i=0;i<timeToPrint;i++) fwrite(&StoredData[i][0],sizeof(double),varnum,outfile);
            fclose(outfile);
            timeToPrint = 0;
        }




        if(DEBUGGING) {
        printf("(%d) %g %g %g : ",myrank,c_Ranges[0][0],c_Ranges[1][0],c_Ranges[2][0]);
        for(i=0;i<varnum;i++) printf("%g ",curData[i]);
        printf("\n");
        }

        // Update Current Cell Number
        CurCellNum = CurCellNum + 1;

        if(CurCellNum >= NextPer)
        {
            perDone = 100*((double)CurCellNum)/((double)CellsPP);
            NextPer = NextPer+PerInc;
            printf("Processor %d is %g percent done.\n",myrank,perDone);
        }




        // Update restart file (overwrite, not append)
#if(RESTARTFILE)
        restartfile = fopen(restartfilename,"w");
        fwrite(&CurCellNum,sizeof(int),1,restartfile);
        fclose(restartfile);
#endif

    }


    // And we're done.
    //fclose(outfile);

    // Book keeping results
    if(CHATTY || DEBUGGING)
    {
        printf("Rank %d -- Total projected mass: %10.10g (%10.10g solar masses)\n",myrank,mass_conv*totProjMass,totProjMass*mass_conv/solarMass);
        for(i=0;i<3;i++) printf("Rank %d -- Total projected momentum %d: %g\n",myrank,i,totProjMom[i]);
    }

    double sumProjMass = 0;
    double sumProjMom[3]  = {0,0,0};
    int CheckSumSum = 0;

    //MPI_Allreduce(&CheckSum2, &CheckSumSum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //MPI_Allreduce(&totProjMass, &sumProjMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //MPI_Allreduce(&totProjMom[0], &sumProjMom[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if(CHATTY && myrank==0)
    {
        //printf("MPISUM: CheckSum: %d \n",CheckSumSum);
        //printf("MPISUM: Total mass projected: %10.10g (%10.10g solar)\n",sumProjMass*mass_conv,sumProjMass*mass_conv/solarMass);
        //printf("MPISUM: Total momentum projected: %g %g %g\n",sumProjMom[0],sumProjMom[1],sumProjMom[2]);
    }

    return 0;
}




/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
                                 Projection Code
 =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */

// New attempt at writing snapshot writing so that it is less memory intensive and arranges the data smartly
int Projection_Bread(char *outname, double pCenter[], double vREL[], int totProc, int myrank)
{

    /* Inputs:
        outname: name of output file
        pCENTER: co-moving coordinates of location to center on
        vREL: co-moving relative velocity to subtract out
        totProc: total number of processors
        myrank: Rank of working processor
    */


    //Iterators
    int i,j,k,v,n;

    // Some conversions that will be useful later
    double vel_conv = 1.e5*pow(Time,0.5);
    double dens_conv = P[10].Rho/P[10].Density; //assumes at least 10 particles...
    double kernel_conv =  pow(pcTOcm*1.e3*Time/(hubble_param),-3);
    double CtoP = 1.e3*Time/(hubble_param);

    // Other counters and such
    int ParticleCountsBread = 0;


    /* Each processor has a unique value for myrank (out of totProc). Use this
     * and the number of processors to determine which grid cells it searches
     * over
     */

    // Using width, break each dimension into rectangles
    double Rect[3] = {width, width, width/( (double)totProc ) };

    // Grid cell length
    double DeltaX = width/ref_lev_doub;
    double DeltaXh = DeltaX/2.0;

    // Cell range for each processor (different for each processor)
    int zCRange[2];
    int xCells,yCells,zCells;
    xCells=ref_lev;
    yCells=ref_lev;
    zCells=ref_lev/totProc; // These are same for each proc

    zCRange[0] = myrank*zCells; // These two lines are different for each proc
    zCRange[1] = (1+myrank)*zCells-1;

    // Number of cells stored on each processor
    int CellsPP = xCells*yCells*zCells;


    // Physical Z distance range for each processor (different for each proc)
    // Centers on middle of cell
    double zRange[2];
    zRange[0] = ((double) myrank)*Rect[2] - width/2.0 + DeltaXh;
    zRange[1] = zRange[0] + Rect[2] - DeltaX;
    //xRange[0] = -width/2.0 + DeltaXh;
    //xRange[1] = xRange[0] + RecX - DeltaX;


    if(CHATTY) printf("myrank = %d, cell ranges for z: (%d,%d)\n",myrank, zCRange[0],zCRange[1]);
    if(CHATTY) printf("myrank = %d, distance ranges for z: (%g,%g)\n",myrank,zRange[0],zRange[1]);



    /* Each processor will then seek over all particles and fill only the
     * cells that it owns. It will then compile all this information in the
     * output file.
     */


    /*  Here we subtract out velocities and shift particle positions so our
        center location is at the "center" of the Orion2 box
        ( "center" = (-1/2,-1/2,-1/2)*DeltaX )
     */
    double gMassTot=0.0;
    int gNcount=0;
    for(n=0;n<Ngas;n++)
    {
        P[n].disx = CtoP*(P[n].Pos[0] - pCenter[0]);// - DeltaXh; // Now in physical units
        P[n].disy = CtoP*(P[n].Pos[1] - pCenter[1]);// - DeltaXh;
        P[n].disz = CtoP*(P[n].Pos[2] - pCenter[2]);// - DeltaXh;
        P[n].Vel[0] = P[n].Vel[0] - vREL[0]; // Still comoving
        P[n].Vel[1] = P[n].Vel[1] - vREL[1];
        P[n].Vel[2] = P[n].Vel[2] - vREL[2];


        //Also add up the mass in the box of interest (why?)
        if(fabs(fabs(P[n].disx) ) < width/2.0 &&
           fabs(fabs(P[n].disy) ) < width/2.0 &&
           fabs(fabs(P[n].disz) ) < width/2.0)
        { gMassTot = gMassTot + P[n].Mass; gNcount = gNcount + 1;}
        // same for each processor, since each processor loops over all parts
    }
    if(CHATTY) printf("Number and Mass total of gadget particles, whose centers are inside the box of interest: N = %d, m_gadget = %g grams (%g solar masses)\n",gNcount,mass_conv*gMassTot,mass_conv*gMassTot/solarMass);


    // Given the way the data will be layed out now, we'll need some header
    // information so we can reconstruct it correctly in Orion2.
    if(myrank==0)
    {
        FILE *outfile;
        outfile=fopen(outname,"w");
        fwrite(&totProc, sizeof(int), 1, outfile);
        fwrite(&ref_lev, sizeof(int), 1, outfile);
        fwrite(&width, sizeof(double), 1, outfile);
        fclose(outfile);
    }


    /* Here is where all the work starts. (varnum is a global variable) */
    double *Gdum ;
    double *GrhoT;
    double *ProjM;
    double *ProjM2;



    // For each variable...
    for(v=-1;v<varnum;v++)      // -1 is the loop for the projected mass ratio
    {

        // If v=-1, we need a different array of data
        if(v==-1)
        {
            if(CHATTY) printf("Allocating array for particle mass projections (rank %d)\n",myrank);
            ProjM = (double *) malloc( Ngas * sizeof(double));
            ProjM2 = (double *) malloc( Ngas * sizeof(double));
            if(CHATTY) printf("Allocating array for particle mass projections DONE! (rank %d)\n",myrank);
            for(n=0;n<Ngas;n++) ProjM[n]=0.0;
            for(n=0;n<Ngas;n++) ProjM2[n]=0.0;
        }

        if(v==0)
        {
            // Allocate memory, a dummy array that holds the current element of interest
            if(CHATTY) printf("Allocating memory!\n");
            Gdum = (double *) malloc( CellsPP * sizeof(double));
            GrhoT = (double *) malloc( CellsPP * sizeof(double));
            if(CHATTY) printf("Allocating memory DONE!\n");

            for(n=0;n<CellsPP;n++) GrhoT[n]=0.0; // Only want to do this once.

        }

        // Clean slate for each processor
        if(v!=-1) for(n=0;n<CellsPP;n++) Gdum[n]=0.0;  // if() statement is b/c this isn't allocated yet when v=-1

            int TouchCount=0;
        // Loop over each particle, dropping them into buckets
        for(n=0;n<Ngas;n++)
        {


            // For each particle, get physical extent of smoothing length
            double hsm = P[n].hsm;
            double p_mass = P[n].Mass * mass_conv;

            // Physical extent of each particle (includes smoothing length)
            double p_xRange[2] = {P[n].disx - hfac*P[n].hsm_phys , P[n].disx + hfac*P[n].hsm_phys };
            double p_yRange[2] = {P[n].disy - hfac*P[n].hsm_phys , P[n].disy + hfac*P[n].hsm_phys };
            double p_zRange[2] = {P[n].disz - hfac*P[n].hsm_phys , P[n].disz + hfac*P[n].hsm_phys };

            // Determines if this particle's extent is on this processor
            int NotHere = 0;

            //if( p_xRange[0] > width/2.0 || p_xRange[1] < -width/2.0 ) NotHere = 1;
            //if( p_yRange[0] > width/2.0 || p_yRange[1] < -width/2.0 ) NotHere = 1;
            //if( p_zRange[0] > zRange[1]+DeltaXh || p_zRange[1] < zRange[0]-DeltaXh  ) NotHere = 1;

            if( p_xRange[0] > width/2.0 || p_xRange[1] < -width/2.0 ) NotHere = 1;
            if( p_yRange[0] > width/2.0 || p_yRange[1] < -width/2.0 ) NotHere = 1;
            if( p_zRange[0] > zRange[1]+DeltaXh || p_zRange[1] < zRange[0]-DeltaXh  ) NotHere = 1;

            // If out of range for one coordinate, 'continue' back to top of loop
            if(NotHere) continue;



            if(v==0) ParticleCountsBread = ParticleCountsBread + 1; // Will be different for each processor
            // Some particles might be double counted, also.

            /*
             OK, at this point we've determined this particular particle extends
             on the current processor. Let's find what cells it extends over
             */


            // This processor doesn't care about what this particle has
            // beyond the processor's domain if we're not calculating the projected mass

            // Z is more complicated since that is how the grid is divided amongst processors

            if(p_xRange[0] < -width/2.0+DeltaXh && v!=-1) p_xRange[0]= -width/2.0+DeltaXh;
            if(p_yRange[0] < -width/2.0+DeltaXh && v!=-1) p_yRange[0]= -width/2.0+DeltaXh;
            if(p_zRange[0] < zRange[0])
            {
                if(v==-1 && myrank==0) p_zRange[0] = p_zRange[0]; // do nothing
                else p_zRange[0] = zRange[0]; // clip it
            }

            if(p_xRange[1] > width/2.0-DeltaXh && v!=-1) p_xRange[1]=width/2.0-DeltaXh;
            if(p_yRange[1] > width/2.0-DeltaXh && v!=-1) p_yRange[1]=width/2.0-DeltaXh;
            if(p_zRange[1] > zRange[1])
            {
                if(v==-1 && myrank==totProc-1) p_zRange[1] = p_zRange[1]; // do nothing
                else p_zRange[1]=zRange[1]; // clip it
            }


            //if(p_xRange[1] > width/2.0-DeltaXh) p_xRange[1]=width/2.0-DeltaXh;
            //if(p_yRange[1] > width/2.0-DeltaXh) p_yRange[1]=width/2.0-DeltaXh;
            //if(p_zRange[1] > zRange[1]) p_zRange[1]=zRange[1];


            // Convert to cell ranges
            int p_xCRange[2]={0,0};
            int p_yCRange[2]={0,0};
            int p_zCRange[2]={0,0};



            if(v==-1)
            {
                // If we have made it to this point, then we know there is some overlap

                // The left side of the smoothing region is at most x=width/2. Start there and work back
                i=ref_lev-1;
                j=zCells-1;
                int xdone,ydone,zdone;
                xdone=ydone=zdone=0;
                while( 1 )
                {
                    if( xdone==0 && p_xRange[0] > -width/2.0 + ((double) i)*DeltaX ){ xdone=1; p_xCRange[0] = i;}
                    if( ydone==0 && p_yRange[0] > -width/2.0 + ((double) i)*DeltaX ){ ydone=1; p_yCRange[0] = i;}
                    if( zdone==0 && p_zRange[0] > zRange[0] - DeltaXh + ((double) j)*DeltaX ) { zdone=1; p_zCRange[0] = myrank*zCells + j;}
                    j = j-1;
                    i = i-1;
                    if(xdone==1 && ydone==1 && zdone==1) break;
                }
                // The right side of the smoothing region is at least x=-width/2. Start there and work back
                zdone=xdone=ydone=0;
                i=0;
                while( 1 )
                {
                    if( xdone==0 && p_xRange[1] < -width/2.0 + (1.0 + ((double) i))*DeltaX ){ xdone=1; p_xCRange[1] = i;}
                    if( ydone==0 && p_yRange[1] < -width/2.0 + (1.0 + ((double) i))*DeltaX ){ ydone=1; p_yCRange[1] = i;}
                    if( zdone==0 && p_zRange[1] < zRange[0] - DeltaXh + (1.0 + ((double)i))*DeltaX ) { zdone=1; p_zCRange[1] = myrank*zCells + i;}
                    i = i+1;
                    if(xdone==1 && ydone==1 && zdone==1) break;
                }

            }
            else
            {

                for(j=0;j<2;j++)
                {
                    for(i=0;i<=ref_lev;i++)
                    {   double i_d = (double)i;
                        if(p_xRange[j] < -width/2.0 + (1.0 + i_d)*DeltaX  )
                        { p_xCRange[j] = i; break; }}
                    for(i=0;i<=ref_lev;i++)
                    {   double i_d = (double)i;
                        if(p_yRange[j] < -width/2.0 + (1.0 + i_d)*DeltaX )
                        { p_yCRange[j] = i; break; }}
                    for(i=0;i<zCells;i++)
                    {   double i_d = (double)i;
                        if(p_zRange[j] < zRange[0] - DeltaXh + (1.0 + i_d)*DeltaX )
                        { p_zCRange[j] = myrank*zCells + i; break; }} // so z also ranges from 0 to ref_lev-1
                }

            }



            if( v!=-1 && (p_zCRange[0] < 0 || p_zCRange[1] > ref_lev-1) ) printf("MAYDAY: WHAT?!!?\n");



            // Loop over all the cells this particle touches and put some mass in that bin
            for(i=p_xCRange[0];i<=p_xCRange[1];i++)
                for(j=p_yCRange[0];j<=p_yCRange[1];j++)
                    for(k=p_zCRange[0];k<=p_zCRange[1];k++)
                    {


                        // We only want to deal with cells that get a non-zero amt. of mass

                        // Need physical location of the cell center for i,j,k
                        double cur_xgrid = -width/2.0 + DeltaXh + ((double)i)*DeltaX;
                        double cur_ygrid = -width/2.0 + DeltaXh + ((double)j)*DeltaX;
                        double cur_zgrid = -width/2.0 + DeltaXh + ((double)k)*DeltaX;

                        double kernel = calc_kernel_spline(n,cur_xgrid,cur_ygrid,cur_zgrid,hsm,DeltaX,DeltaXh);
                        double kernel_phys = kernel*kernel_conv;

                        // particle contributes to cell
                        if(v==0) TouchCount = TouchCount+1;

                        double rho = p_mass*kernel_phys; // rho_ij (g/cm^3)
                        if(rho<=0 && cur_zgrid >= P[n].disz) { k=p_zCRange[1]; continue; }



                        double fac = 0.0; //(P[n].Mass/P[n].Density)*kernel; //no conversions
                        if(v!=-1) fac = (p_mass/P[n].mapped_mass);

                        // The z range possibly needs to be reduced for storage
                        int ThisZ = k;
                        while(1)
                        {
                            if(ThisZ-zCells<0) break;
                            else ThisZ = ThisZ - zCells;
                        }


                        //1D index of current cell of interest for processor memory
                        int idx = i + j*xCells + ThisZ*xCells*yCells; // memory location
                        //assert(idx < xCells*yCells*zCells);

                        // If the cell gets zero density, skip it
                        if(rho>0)
                        {
                            if(v==-1) // How much mass would've been projected?
                            {
                                ProjM[n] = ProjM[n] + rho*pow(pcTOcm*DeltaX,3);
                                //if(P[n].edge_flag) ProjM[n] = ProjM[n] + ProjMCorrect();
                            }
                            else if(v==0) // Density
                            {
                                Gdum[idx] = Gdum[idx] + rho*fac;
                            }
                            else if(v==1 || v==2 || v==3) // Les Momentums
                            {
                                //Gdum[idx] = Gdum[idx] + P[n].Vel[v-1]*fac*vel_conv; //orig
                                Gdum[idx] = Gdum[idx] + P[n].Vel[v-1]*vel_conv*(rho*fac/GrhoT[idx]);
                            }
                            else if(v==4) // Energy Density
                            {
                                //Gdum[idx] = Gdum[idx] + P[n].Pres*fac;
                                Gdum[idx] = Gdum[idx] + P[n].U*(rho*fac/GrhoT[idx]);
                            }
                            else if(v==5) // Chem Species
                            {
                                //Gdum[idx] = Gdum[idx] + P[n].H2I*fac;
                                Gdum[idx] = Gdum[idx] + P[n].H2I*(rho*fac/GrhoT[idx]);
                            }
                            else if(v==6)
                            {
                                //Gdum[idx] = Gdum[idx] + P[n].HDI*fac;
                                Gdum[idx] = Gdum[idx] + P[n].HDI*(rho*fac/GrhoT[idx]);
                            }
                            else if(v==7)
                            {
                                //Gdum[idx] = Gdum[idx] + P[n].HII*fac;
                                Gdum[idx] = Gdum[idx] + P[n].HII*(rho*fac/GrhoT[idx]);
                            }

                        } // End of if(rho>0)


                    } // End of triple for() loop



        } // End of for(every particle) loop
        if(v==0) printf("Total touches = %d (proc %d)\n",TouchCount,myrank);
        MPI_Barrier(MPI_COMM_WORLD);
        // At this point, all particles have been put into buckets.


        if(v==0) {for(n=0;n<CellsPP;n++) GrhoT[n] = Gdum[n];} // Copy of rho for other variables

        // Loop over each processor and dump the data

        if(v==-1)
        {
            printf("Reducing projected mass (proc %d)\n",myrank);
            //Send projected mass data to all other processors
            MPI_Allreduce(&ProjM[0],&ProjM2[0],Ngas,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            //MPI_Reduce(&ProjM[0],&ProjM[0],Ngas,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            //MPI_Bcast(&ProjM[0], Ngas, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        }
        else
        {
            for(i=0;i<totProc;i++)
            {
                if(myrank == i && i != 0)
                {
                    //Send data if it's your turn, else do nothing.
                    //No need to send data if i=0, Gdum already ready to go.
                    printf("Sending data on variable %d from %d to 0\n",v,i);
                    //for(j=0;j<CellsPP;j++) printf("Sending j = %d of %g\n",j,Gdum[j]);
                    MPI_Send( &Gdum[0], CellsPP, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                }
                else if(myrank==0 && i != 0)
                {
                    // Receive yummy data (overwrites!)
                    printf("Receiving data regarding variable %d from %d!\n",v,i);
                    MPI_Recv( &Gdum[0], CellsPP, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    //for(j=0;j<CellsPP;j++) printf("Receiving j = %d of %g\n",j,Gdum[j]);
                }

                // Proc 0 now has the next set of data to print
                if(myrank==0) PrintAway(Gdum,CellsPP,outname,i,width,ref_lev);
                MPI_Barrier(MPI_COMM_WORLD);

            }
        }


        if(v==-1)
        {
            for(n=0;n<Ngas;n++) P[n].mapped_mass = ProjM2[n];
            //for(n=0;n<Ngas;n++) printf("P[%d].mapped_mass = %g\n",n,P[n].mapped_mass);
            printf("Deallocating mass projection arrays (rank %d)\n",myrank);
            //And we're done with these guys
            free(ProjM);
            free(ProjM2);
        }

        if(v==0)
        {
            MPI_Reduce(&ParticleCountsBread,&ParticleCounts,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD); // still might double count some edge particles...

            double totProjMass = 0.0;
            double masssum = 0.0;
            for(i=0;i<CellsPP;i++) totProjMass = totProjMass + GrhoT[i]*pow(pcTOcm*DeltaX,3);
            MPI_Allreduce(&totProjMass,&masssum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            printf("Total projected mass = %g (rank %d mass = %g)\n",masssum,myrank,totProjMass);
        }


    } // end of for(each variable) loop



    free(Gdum);
    free(GrhoT);

    // Everything is written to the file, we're done! (phew...)
    return 0;
}
