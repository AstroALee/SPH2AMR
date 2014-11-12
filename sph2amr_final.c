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


/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */
int main(int argc, char **argv)
{
    // Ready... Set... GO!
    clock_t timeMe;
    timeMe = clock();
    
    
    // Initialize the MPI Universe
    int npes, myrank, ierr=0;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &npes); // number processors
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank); // labels for each processor
    if(CHATTY || DEBUGGING) printf("npes = %d\n",npes);
    
    /* A few important input values are needed from the command line
     *
     * ref_lev : An integer giving the number of grid cells in each dimension
     *           for the resulting box
     * width   : A double given the size of each box side, in units of parsecs
     */
    
    
    if(argc < 3) // Double check on inputs
    {
        printf("WATERLOO: 2 inputs needed: (int)ref_lev, and (double)width [unit: pc]\nYou gave %d.\n",argc-1);
        exit(-1);
    }
    
    int ref_lev = atoi(argv[1]);
    double width = atof(argv[2]);
    
    if(CHATTY || DEBUGGING) printf("ref_lev = %d and width = %g pc, argc = %d\n",ref_lev,width,argc);
    
    
#if BREAD
    // Does the number of processors divide evenly into the z-coordinate?
    if(ref_lev%npes != 0)
    {
        printf("WATERLOO: For the bread model, we need the value of ref_lev to be divided evenly by the number of processors \n");
        exit(-1);
    }
#endif
    
    
    
    /* I/O Information */
    
    char path_in[200], path_out[200], basename[200], input_fname[200], input_fname2[200];
    
    //If we are taking a particular snapshot number, identify it here
    int snapshot_number = 0;
    
    //Mulitiple files to load Gadget data from?
    int files=1;

    sprintf(path_in, pathname_in);
    sprintf(path_out, pathname_out);
    sprintf(basename, "bin_HR10");
    sprintf(input_fname, "%s/%s_%03d", path_in, basename, snapshot_number);
    
    
    // Time to read in the Gadget2 data to assemble the list of all particles
    // Each processor will have access to this data
    double DelCoord[3] = {0,0,0};
    double boxsize = 0;
    Ngas = read_snapshot(input_fname, files, DelCoord, boxsize);
    

    // If dealing with B-field from analytic calculations, read that in here.
#if (readB)
    sprintf(input_fname2, "%s_bfield.dat", path_in);
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
    
    
    
    
    
    
    
    
    char  output_fname[200] ;
    int  j=0, n, type, Ngas, Ngas2, random;
    double x,y,z,x1,y1, nh, nhmax, ref_lev_doub, grid_size, grid_size_half;
    double delx, dely, delz, dis, disx, disy, disz; //, disAU;
    double Part_Mtot;
    double pCOM[3],vCOM[3]; // position and velocity of COM
    double pDEL[3];
    FILE *outfile, *infile;
    
    
 
    
    //Initializations, no need to touch these
    boxsize = 140.0;
    delx = dely = delz = 0;
    ref_lev_doub = (double)ref_lev;
    grid_size = width/ref_lev_doub;
    grid_size_half = grid_size/2.0;
    for(j=0;j<3;j++) pCOM[j]=0.0;
    for(j=0;j<3;j++) vCOM[j]=0.0;
    for(j=0;j<3;j++) pDEL[j]=0.0;
    
   
    
    
   
   
    
    
    
    
    
    
    
    
    // Calculate total mass
    Part_Mtot = 0.0;
    for(n=0;n<Ngas;n++) Part_Mtot = Part_Mtot + P[n].Mass;
    if(myrank==0 && CHATTY) printf("MASSTEST: The mass in the entire simulation box is = %g\n",mass_conv*Part_Mtot);
    
    /* Calculate comoving COM position of entire domain (might want to center
     * on this at some point?
     */
    
    for(n=0;n<Ngas;n++)
    {
        pCOM[0] = pCOM[0] + P[n].Pos[0]*P[n].Mass;
        pCOM[1] = pCOM[1] + P[n].Pos[1]*P[n].Mass;
        pCOM[2] = pCOM[2] + P[n].Pos[2]*P[n].Mass;
    }
    pCOM[0]=CtoP*pCOM[0]/Part_Mtot; // Converted to physical units
    pCOM[1]=CtoP*pCOM[1]/Part_Mtot; // (tecnically, the masses never were, but it was
    pCOM[2]=CtoP*pCOM[2]/Part_Mtot; //  divided out)
    
    
    // Calculate location of Density maximum for "re-centering" of Orion2 data
    nhmax = 0;
    for(n=0;n<Ngas;n++)
    {
        nh = P[n].nh;
        if(nh > nhmax && P[n].sink > -1)
        {
            nhmax = nh;
            delx = P[n].Pos[0];  //del's based upon location of maximum density
            dely = P[n].Pos[1];
            delz = P[n].Pos[2];
        }
    }
    pDEL[0] = delx; // These are used later, keep in comoving coordinates
    pDEL[1] = dely;
    pDEL[2] = delz;
    
    if(CHATTY) printf("physical location of COM and denest location: (%g,%g,%g) , (%g,%g,%g)\n",pCOM[0],pCOM[1],pCOM[2],CtoP*pDEL[0],CtoP*pDEL[1],CtoP*pDEL[2]);
    
    
    // Calculate the center of mass velocity for particles in box of interest
    //double InterestMtot = 0.0; // only mass in box we're projecting for
    int Interestcount = 0;
    for(n=0;n < Ngas; n++)
    {
        // Calculate distance from densest point, converting to physical units
        dis = pow(((P[n].Pos[0]-delx)*(P[n].Pos[0]-delx)
                   + (P[n].Pos[1]-dely)*(P[n].Pos[1]-dely)
                   + (P[n].Pos[2]-delz)*(P[n].Pos[2]-delz)), 0.5);
        dis=dis*CtoP;
        
        //double DeltaXh = width/((double)ref_lev);
        // 1D distances
        disx = fabs( (P[n].Pos[0]-delx)*CtoP );
        disy = fabs( (P[n].Pos[1]-dely)*CtoP );
        disz = fabs( (P[n].Pos[2]-delz)*CtoP );
        
        if(disx < width/2.0 && disy < width/2.0 && disz < width/2.0) // if in box we want
        {
            vCOM[0] = vCOM[0] + P[n].Vel[0]*P[n].Mass;
            vCOM[1] = vCOM[1] + P[n].Vel[1]*P[n].Mass;
            vCOM[2] = vCOM[2] + P[n].Vel[2]*P[n].Mass;
            InterestMtot = InterestMtot + P[n].Mass;
            Interestcount = Interestcount + 1;
        }
    }
    vCOM[0] = vCOM[0]/InterestMtot; // Kept in co-moving units
    vCOM[1] = vCOM[1]/InterestMtot;
    vCOM[2] = vCOM[2]/InterestMtot;
    
    InterestMtot = mass_conv*InterestMtot;
    
    if(CHATTY) printf("Total gadget particles %d, m_gadget = %g\n",Interestcount,InterestMtot);
    
    if(myrank==0)
    {
        if(CHATTY) for(n=0;n<3;n++) printf("comoving vCOM[%d] = %lg , ",n,vCOM[n]);
        if(CHATTY) printf("\n");
    }
    
    
    
    sprintf(output_fname, "%s/gadget2orion", path_out);
    if(myrank == 0) // creates file (overwrites!)
    {
        outfile=fopen(output_fname,"w");
        fclose(outfile);
    }
    
    
    
#if BREAD
    // Bread model
    Ngas2 = write_snapshotLessMemBread(input_fname, files, output_fname, pDEL, vCOM, Ngas, npes, myrank,ref_lev,width);
#else
    
    // A method of trying to use less memory, but breaks the data up into rectangular chunks. Output is not in any sensible arrangement, which poses an additional challenge at read-in. Not throughly tested, either. Use with caution.
    //Ngas2 = write_snapshotLessMem(input_fname, files, output_fname, pDEL, vCOM, Ngas, npes, myrank,ref_lev,width);
    
    // The original memory intensive method.
    Ngas2 = write_snapshot(input_fname, files, output_fname, delx, dely, delz, vCOM, Ngas, myrank,ref_lev,width);
#endif
    
    
    if(CHATTY) if(myrank==0) printf("Particle Count = %d\n",ParticleCounts);
    
    // Tick tock
    timeMe = clock() - timeMe;
    printf("Processor %d took %g wall-seconds on %d processors (dim,width = %d,%g pc).\n",myrank,((double)timeMe)/CLOCKS_PER_SEC,npes,ref_lev,width);
    
    
    
    
    
    
    
    // Clean up
    ierr=MPI_Finalize();
    
} // End of main()





