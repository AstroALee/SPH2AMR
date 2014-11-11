/*
 * SPH2AMR.H
 *
 * Reads in Gadget2 data and projects the data onto a uniform grid suitable for
 * Orion2, which can be read in during initialization
 * Written by Athena Stacy and Aaron Lee, 2014
 *
 */
 




/* PreProcessor Directives =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */

/* Debugging */
#define DEBUGGING 0

/* Diagnostic Info? (Automatically included if DEBUGGING) */
#define CHATTY 1

/* Read in B-field information? */
#define READB 0

/* Use less memory intensive Bread method? */
#define BREAD 1 


/* Where do I look for the Gadget2 output file? Where do I print the result
 file? Simply a period means the current working directory. */
#define pathname_in "."
#define pathname_out "."


/* Some other pre-processors */
#define MAXREF 20 // Gadget2-related, probably never have to touch this


/* End PreProcessor Directives =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */


/* Let me see all those libraries */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mpi.h"

#if(DEBUGGING)
#include <cassert>
#endif

/* Needed for cosmological situations H_0 = 100*h */
#define hubble_param 0.7
#define PI 3.14159265359


/* =-=-=-=-=-=-=-=-=-=-=-=- Global Constants =-=-=-=-=-=-=-=-=-=-=-=- */

// Global Unit Conversions
double pcTOcm = 3.08567758e18;


void PrintAway(double Gdum[],int CellsPP,char *outname, int curRank, double width, int ref_lev);
void BreakUpDomain(int Cord[], int npes);
int read_snapshot(char *fname, int files, char *outname, double delx,
                  double dely, double delz, double  boxsize);
int write_snapshotLessMemBread(char *fname, int files, char *outname, double pDEL[],
                               double vCOM[], int NumGas, int npes, int proc, int ref_lev,
                               double width);
int write_snapshot(char *fname, int files, char *outname, double delx,
                   double dely, double delz, double vCOM[], int NumGas,
                   int myrank, int ref_lev, double width);
int reordering(void);
int unit_conversion(void);
int allocate_memory(void);
double calc_kernel_spline(int n, double x, double y, double z, double hsm,
                          double grid_size, double grid_size_half);
double calc_kernel_tsc(int n, double x, double y, double z, double hsm,
                       double grid_size, double grid_size_half);


struct io_header_1
{
    int      npart[MAXREF];
    double   mass[MAXREF];
    double   time;
    double   redshift;
    int      flag_sfr;
    int      flag_feedback;
    int      npartTotal[MAXREF];
    int      flag_cooling;
    int      num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam;
    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8]; //fills to 256 Bytes
} header1;

struct io_header_old
{
    int      npart[MAXREF];
    double   mass[MAXREF];
    double   time;
    double   redshift;
    int      flag_sfr;
    int      flag_feedback;
    int      npartTotal[MAXREF];
    int      flag_cooling;
    int      num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam;
    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8]; //fills to 256 Bytes
} header_old;

struct particle_data
{
    double  Pos[3];
    double  Vel[3];
    double  Mass;
    int    Type;
    double disx, disy, disz;
    double  Rho, U, Pres, nh, Density, hsm, hsm_phys;
    double mapped_mass;
    //double Temp, sink;
    //double  elec, HI, HII, HeI, HeII,  HeIII, H2I, H2II, HM, hsm, DI, \
    // DII, HDI, DM, HDII, FosHII gam;
    double H2I, HII, DII, HDI, HeII, HeIII, gam, sink;
    double nh_test;
#if (READB)
    double nh_test, Bfieldx, Bfieldy, Bfieldz;
#endif
    double dummy;
    int edge_flag;
} *P;


// Global variables
int *Id;
int NumPart, Ngas, NumEPart;
double Time, zred;
int varnum = 8;
double InterestMtot = 0;
int ParticleCounts = 0;




/* =-=-=-=-=-=-=-=- Other functions =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */




// New attempt at writing snapshot writing so that it is less memory intensive and arranges the data smartly
int write_snapshotLessMemBread(char *fname, int files, char *outname, double pDEL[], double vCOM[], int NumGas, int totProc, int myrank, int ref_lev, double width)
{
    
    //Iterators
    int i,j,k,v,n;
    
    // Some numbers that will be useful later
    double vel_conv = 1.e5*pow(Time,0.5);
    double dens_conv = P[10].Rho/P[10].Density; //assumes at least 10 particles...
    double mass_conv = (1.e10/hubble_param)*1.989e33;
    double kernel_conv =  pow(pcTOcm*1.e3*Time/(hubble_param),-3);
    double CtoP = 1.e3*Time/(hubble_param);
    double hfac = 2.0;
    
    double ref_lev_doub = (double)ref_lev;
    
    double delvx = vCOM[0]; // assuming comoving coordinates
    double delvy = vCOM[1];
    double delvz = vCOM[2];
    double delx  = pDEL[0]; // assuming comoving coordinates
    double dely  = pDEL[1];
    double delz  = pDEL[2];
    
    
    /* Each processor has a unique value for myrank (out of totProc). Use this
     * and the number of processors to determine which grid cells it searches
     * over
     */
    
    // Using width, break each dimension into rectangles
    double RecX = width;
    double RecY = width;
    double RecZ = width/( (double)totProc );
    
    // Grid cell length (same as grid_size and grid_size_half)
    double DeltaX = width/ref_lev_doub;
    double DeltaXh = DeltaX/2.0;
    
    
    // Cell range for each processor (different for each processor)
    int zCRange[2];
    int xCells,yCells,zCells;
    xCells=ref_lev;
    yCells=ref_lev;
    zCells=ref_lev/totProc; // These are same for each proc
    
    zCRange[0] = myrank*zCells;
    zCRange[1] = (1+myrank)*zCells-1; // This line is different for each proc
    
    // Number of cells stored on each processor
    int CellsPP = xCells*yCells*zCells;
    
    
    // Physical Z distance range for each processor (different for each proc)
    // Centers on middle of cell
    double zRange[2];
    zRange[0] = ((double) myrank)*RecZ - width/2.0 + DeltaXh;
    zRange[1] = zRange[0] + RecZ - DeltaX;
    //xRange[0] = -width/2.0 + DeltaXh;
    //xRange[1] = xRange[0] + RecX - DeltaX;
    
    
    if(CHATTY) printf("myrank = %d, cell ranges for z: (%d,%d)\n",myrank, zCRange[0],zCRange[1]);
    if(CHATTY) printf("myrank = %d, distance ranges for z: (%g,%g)\n",myrank,zRange[0],zRange[1]);
    
    
    
    /* Each processor will then seek over all particles and fill only the
     * cells that it owns. It will then compile all this information in the
     * output file.
     */
    
    
    /*  Here we subtract out COM velocities and shift particle positions so
     the densest cell is at the "center" of the Orion2 box
     ( "center" = (-1/2,-1/2,-1/2)*DeltaX )
     */
    double gMassTot=0.0;
    int gNcount=0;
    for(n=0;n<NumGas;n++)
    {
        P[n].hsm_phys = P[n].hsm*CtoP;
        P[n].disx = CtoP*(P[n].Pos[0] - delx);// - DeltaXh; // Now in physical units
        P[n].disy = CtoP*(P[n].Pos[1] - dely);// - DeltaXh;
        P[n].disz = CtoP*(P[n].Pos[2] - delz);// - DeltaXh;
        P[n].Vel[0] = P[n].Vel[0] - delvx; // Still comoving
        P[n].Vel[1] = P[n].Vel[1] - delvy;
        P[n].Vel[2] = P[n].Vel[2] - delvz;
        
        
        //Also add up the mass in the box of interest
        if(fabs(fabs(P[n].disx) ) < width/2.0 &&
           fabs(fabs(P[n].disy) ) < width/2.0 &&
           fabs(fabs(P[n].disz) ) < width/2.0)
        { gMassTot = gMassTot + P[n].Mass; gNcount = gNcount + 1;}
        // same for each processor, since each processor loops over all parts
    }
    if(CHATTY) printf("Number and Mass total of gadget particles, summed at cell centers. N = %d, m_gadget = %g\n",gNcount,mass_conv*gMassTot);
    
    
    /*
     // Determines which particles reside on the edge of the domain. Since this is the same for all processors, no
     // need to sync
     
     for(n=0;n<NumGas;n++)
     {
     double smooth = P[n].hsm_phys;
     double x_left = P[n].disx - smooth;
     double x_cent = x_left + smooth;
     double x_right = x_cent + smooth;
     
     
     if( d_left < -width/2.0  && d_right > -width/2.0) {P[n].edge_flag=1;continue;}
     if( d_left < width/2.0  && d_right > width/2.0) {P[n].edge_flag=1;continue;}
     
     d_left = P[n].disy - P[n].hsm_phys;
     d_right = P[n].disy + P[n].hsm_phys;
     if( d_left < -width/2.0  && d_right > -width/2.0) {P[n].edge_flag=1;continue;}
     if( d_left < width/2.0  && d_right > width/2.0) {P[n].edge_flag=1;continue;}
     
     d_left = P[n].disz - P[n].hsm_phys;
     d_right = P[n].disz + P[n].hsm_phys;
     if( d_left < -width/2.0  && d_right > -width/2.0) {P[n].edge_flag=1;continue;}
     if( d_left < width/2.0  && d_right > width/2.0) {P[n].edge_flag=1;continue;}
     
     }
     for(n=0;n<NumGas;n++) if(P[n].edge_flag==1) NumEPart = NumEPart +1;
     if(CHATTY) printf("Number of edge particles = %d\n",NumEPart);
     */
    
    // Now each processor has the same information for which particles creep over edges
    
    
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
    
    int ParticleCountsBread = 0;
    
    
    
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
            ProjM = (double *) malloc( NumGas * sizeof(double));
            ProjM2 = (double *) malloc( NumGas * sizeof(double));
            if(CHATTY) printf("Allocating array for particle mass projections DONE! (rank %d)\n",myrank);
            for(n=0;n<NumGas;n++) ProjM[n]=0.0;
            for(n=0;n<NumGas;n++) ProjM2[n]=0.0;
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
        
        // Loop over each particle, dropping them into buckets
        for(n=0;n<NumGas;n++)
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
                        
                        double rho = p_mass*kernel_phys; // rho_ij (g/cm^3)
                        
                        
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
        MPI_Barrier(MPI_COMM_WORLD);
        // At this point, all particles have been put into buckets.
        
        
        if(v==0) {for(n=0;n<CellsPP;n++) GrhoT[n] = Gdum[n];} // Copy of rho for other variables
        
        // Loop over each processor and dump the data
        
        if(v==-1)
        {
            printf("Reducing projected mass (proc %d)\n",myrank);
            //Send projected mass data to all other processors
            MPI_Allreduce(&ProjM[0],&ProjM2[0],NumGas,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            //MPI_Reduce(&ProjM[0],&ProjM[0],NumGas,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            //MPI_Bcast(&ProjM[0], NumGas, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
            for(n=0;n<NumGas;n++) P[n].mapped_mass = ProjM2[n];
            //for(n=0;n<NumGas;n++) printf("P[%d].mapped_mass = %g\n",n,P[n].mapped_mass);
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
    //free(Grho_tot);
    // Everything is written to the file, we're done! (phew...)
    
    return 0;
}


void PrintAway(double Gdum[],int CellsPP, char *outname, int curRank,double width, int ref_lev)
{
    // Inputed a lot of information in case we need to reconstruct cell locations
    // (perhaps to double check the edges are not being double counted... might have
    //  been smart enough above so this is unnecessary.)
    FILE *outfile;
    int k=0;
    
    //printf("Printing\n");
    
    outfile=fopen(outname,"a");
    for(k=0;k<CellsPP;k++) fwrite(&Gdum[k], sizeof(double), 1, outfile);
    for(k=0;k<CellsPP;k++) printf("%d : %g\n",k,Gdum[k]);
    fclose(outfile);
    
    //for(k=0;k<CellsPP;k++) printf("Gdum printed = %g\n",Gdum[k]);
    
}


void BreakUpDomain(int Cord[], int curProc)
{
    //Loop over each dimension distributing 2's, 3's, 5's, etc.
    int i=0,j=0,k=0;
    int Primes[5]={2,3,5,7,11}; // Can have more if you want... but why?!
    
    int istart = 0;
    
    for(i=istart;i<3;i++)
    {
        int done=5;
        for(j=0;j<5;j++)
        {
            if(curProc%(Primes[j])==0)
            {
                Cord[i]=Cord[i]*Primes[j];
                curProc=curProc/Primes[j];
                break;
            }
            else
            {
                done = done-1;
            }
        }
        if(done==0)
        {
            // Higher order prime number left (why did you chose that?!)
            // Multiply what's left to the coordinate with the min value
            int idxMin=0;
            
            for(k=1;k<3;k++)
            {
                if( Cord[k] < Cord[idxMin])
                {
                    idxMin=k;
                }
            }
            Cord[idxMin]=Cord[idxMin]*curProc;
            curProc=curProc/curProc; // = 1
        }
        
    }
    
    if(curProc!=1) BreakUpDomain(Cord,curProc);
    
}






int write_snapshot(char *fname, int files, char *outname, double delx, double dely, double delz, double vCOM[], int NumGas, int myrank, int ref_lev, double width)
{
    FILE *outfile;
    char   buf[200];
    int i, nmin, nmax, j, k, l, n, send_int, send_size, tot_proc=1, tot_proc_sum;
    int imin, imax, jmin, jmax, kmin, kmax;
    int NumGas_per_proc, vartype, varnum=8;
    double ref_lev_doub, i_doub, vel_conv, dens_conv, mass_conv, kernel_conv, kernel_phys, mass, rho, hsm;
    double fac, kernel, rad, ratio, xgrid, ygrid, zgrid, grid_size, grid_size_half;
    double grid_arr[ref_lev], masstot, hfac = 2.0;
    double *Gdum, *Grho, *Grho_tot, *Gvelx, *Gvely, *Gvelz, *Gpres, *GH2I, *GHDI, *GHII, *Gtot; //, *Gnpart, *Gnpart_tot;
    MPI_Status status;
    double binnedMass = 0;
    
    double delvx = vCOM[0];
    double delvy = vCOM[1];
    double delvz = vCOM[2];
    
    if(READB == 1)
        varnum=12;
    
    
    //For NON-cosmological runs
    //Time = 1.0;
    
    Gdum =  (double *) malloc(pow(ref_lev,3) * sizeof(double));
    Grho =  (double *) malloc(pow(ref_lev,3) * sizeof(double));
    Gtot =  (double *) malloc(pow(ref_lev,3) * sizeof(double));
    //Gnpart =(double *) malloc(pow(ref_lev,3) * sizeof(double));
    //Gnpart_tot =(double *) malloc(pow(ref_lev,3) * sizeof(double));
    Grho_tot =(double *) malloc(pow(ref_lev,3) * sizeof(double));
    
    printf("myrank = %d, Memory allocation done\n", myrank);
    
    ref_lev_doub = (double)ref_lev;
    send_size = pow(ref_lev,3);
    vel_conv = 1.e5*pow(Time,0.5);
    dens_conv = P[100].Rho/P[100].Density;
    mass_conv = (1.e10/hubble_param)*1.989e33;
    kernel_conv =  pow(pcTOcm*1.e3*Time/(hubble_param),-3);
    grid_size = width/ref_lev_doub;
    grid_size_half = grid_size/2.0;
    
    MPI_Allreduce(&tot_proc, &tot_proc_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    printf("myrank = %d, tot_proc_sum = %d\n", myrank, tot_proc_sum);
    
    double CtoP =1.e3*Time/hubble_param;
    //subtract out COM velocities and shift particle positions to center of Orion2 box
    for(n=0;n<NumGas;n++)
    {
        P[n].hsm_phys = P[n].hsm*1.e3*Time/hubble_param;
        P[n].disx = 1.e3*Time/(hubble_param)*(P[n].Pos[0] - delx);// - grid_size_half;
        P[n].disy = 1.e3*Time/(hubble_param)*(P[n].Pos[1] - dely) ;//- grid_size_half;
        P[n].disz = 1.e3*Time/(hubble_param)*(P[n].Pos[2] - delz) ;//- grid_size_half;
        P[n].Vel[0] = P[n].Vel[0] - delvx;
        P[n].Vel[1] = P[n].Vel[1] - delvy;
        P[n].Vel[2] = P[n].Vel[2] - delvz;
        P[n].mapped_mass = 0;
        if(fabs(fabs(P[n].disx) - 0.0*P[n].hsm_phys) < width/2.0 && fabs(fabs(P[n].disy) - 0.0*P[n].hsm_phys) < width/2.0 && fabs(fabs(P[n].disz) - 0.0*P[n].hsm_phys) < width/2.0)
            masstot = masstot + P[n].Mass;
        
        //if(n==282575) printf("OINK Particle position: %g,%g,%g\n",P[n].disx,P[n].disy,P[n].disz);
        
    }
    
    
    for(i=0;i<ref_lev;i++)
    {
        i_doub = (double)i;
        grid_arr[i] = (width * i_doub/ref_lev_doub) - width/2.0;
    }
    
    NumGas_per_proc = (int) (NumGas/tot_proc_sum);
    nmin = NumGas_per_proc * (myrank);
    nmax = NumGas_per_proc * (myrank+1) -1;
    
    
    if(myrank == tot_proc_sum - 1)  //make sure highest-rank processor actually covers final SPH particle
        nmax = NumGas-1;
    if(nmin >= NumGas-1)            //make sure nmin and nmax for each processor does not fall out of range of 0 to NumGas-1
        nmin = nmax = NumGas-1;
    if(nmax >= NumGas-1)
        nmax = NumGas-1;
    
    
    printf("myrank = %d, nmin = %d, nmax = %d, ref_lev = %d, ref_lev_doub = %lg\n", myrank, nmin, nmax, ref_lev, ref_lev_doub);
    
    // For each processor, find particles whose smoothing lengths overlap with grid cell centers, and add their weighted values to the grid
    l=0;
    for(i=0;i<ref_lev;i++)
        for(j=0;j<ref_lev;j++)
            for(k=0;k<ref_lev;k++)
            {
                Grho[l] = 0.0;
                Gtot[l] = 0.0;
                //Gnpart[l] = 0.0;
                //Gnpart_tot[l] = 0.0;
                Grho_tot[l] = 0.0;
                l++;
            }
    
    for(vartype=-1;vartype<varnum;vartype++)
    {
        
        l=0;
        for(i=0;i<ref_lev;i++)
            for(j=0;j<ref_lev;j++)
                for(k=0;k<ref_lev;k++)
                {
                    Gdum[l] = 0.0;
                    l=l+1;            // I don't think it was incrementing before?
                }
        
        for(n=nmin;n<=nmax;n++)
        {
            hsm = P[n].hsm;
            mass = P[n].Mass*mass_conv;
            
            imin = (int)((P[n].disx - hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);
            imax = (int)((P[n].disx + hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);
            jmin = (int)((P[n].disy - hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);
            jmax = (int)((P[n].disy + hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);
            kmin = (int)((P[n].disz - hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);
            kmax = (int)((P[n].disz + hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);
            
            if(imin < 0)
                imin = 0;
            if(imax > ref_lev-1)
                imax = ref_lev-1;
            if(jmin < 0)
                jmin = 0;
            if(jmax > ref_lev-1)
                jmax = ref_lev-1;
            if(kmin < 0)
                kmin = 0;
            if(kmax > ref_lev-1)
                kmax = ref_lev-1;
            
            if(imin > ref_lev-1 || jmin > ref_lev-1 || kmin > ref_lev-1 || imax < 0 || jmax < 0 || kmax < 0)
                continue;
            
            //if(vartype==0) printf("ZOOP: part %d indices:  (%d,%d) (%d,%d) (%d,%d)\n",n,imin,imax,jmin,jmax,kmin,kmax);
            if(vartype==0) ParticleCounts = ParticleCounts + 1;
            
            for(i=imin;i<=imax;i++)
                for(j=jmin;j<=jmax;j++)
                    for(k=kmin;k<=kmax;k++)
                    {
                        //l = k + (j)*(ref_lev) + (i)*(ref_lev)*(ref_lev);
                        l = i + (j)*(ref_lev) + (k)*(ref_lev)*(ref_lev);
                        
                        if(vartype == -1)
                        {
                            xgrid = grid_arr[i]+grid_size_half;
                            ygrid = grid_arr[j]+grid_size_half;
                            zgrid = grid_arr[k]+grid_size_half;
                            
                            kernel = calc_kernel_spline(n, xgrid, ygrid, zgrid, hsm, grid_size, grid_size_half);
                            kernel_phys = kernel*kernel_conv;
                            
                            P[n].mapped_mass = P[n].mapped_mass + mass*kernel_phys*pow(grid_size*pcTOcm, 3);
                        }
                        
                        if(vartype == 0)
                        {
                            xgrid = grid_arr[i]+grid_size_half;
                            ygrid = grid_arr[j]+grid_size_half;
                            zgrid = grid_arr[k]+grid_size_half;
                            
                            kernel = calc_kernel_spline(n, xgrid, ygrid, zgrid, hsm, grid_size, grid_size_half);
                            //kernel = calc_kernel_tsc(n, xgrid, ygrid, zgrid, hsm, grid_size, grid_size_half);
                            kernel_phys = kernel*kernel_conv;
                            
                            rho = mass*kernel_phys;
                            
                            //fac = pow(grid_size*pcTOcm, -3) / kernel_phys;
                            
                            //fac = P[n].Mass/P[n].Density*kernel;
                            fac = mass / P[n].mapped_mass;
                            
                            //fac = 1.0;
                            
                            if(kernel_phys > 0 && i%1000 == 0)
                                printf("n = %d, mass_mapped = %lg, fac = %lg \n", n, P[n].mapped_mass, fac);
                            
                            
                            if(rho > 0)
                            {
                                //Grho[l]  = Grho[l] + rho;
                                //Grho[l]  = Grho[l] + P[n].Rho*fac;
                                Grho[l]  = Grho[l] + rho*fac;
                                //Gnpart[l] = Gnpart[l] + 1.0;
                            }
                        }
                        
                        
                        if(vartype > 0 && Grho[l] > 0)
                        {
                            xgrid = grid_arr[i]+grid_size_half;
                            ygrid = grid_arr[j]+grid_size_half;
                            zgrid = grid_arr[k]+grid_size_half;
                            
                            kernel = calc_kernel_spline(n, xgrid, ygrid, zgrid, hsm, grid_size, grid_size_half);
                            kernel_phys = kernel*kernel_conv;
                            
                            
                            rho = mass*kernel_phys;
                            
                            //fac = P[n].Mass/P[n].Density*kernel;
                            fac = mass / P[n].mapped_mass;
                            
                            if(vartype == 1 && rho > 0)
                                Gdum[l] = Gdum[l] + P[n].Vel[0]*vel_conv*rho*fac; //Go ahead and convert velocity to cgs units
                            if(vartype == 2 && rho > 0)
                                Gdum[l] = Gdum[l] + P[n].Vel[1]*vel_conv*rho*fac;
                            if(vartype == 3 && rho > 0)
                                Gdum[l] = Gdum[l] + P[n].Vel[2]*vel_conv*rho*fac;
                            if(vartype == 4 && rho > 0)
                                Gdum[l] = Gdum[l] + P[n].U*rho*fac; //P[n].Pres*rho*fac;
                            if(vartype == 5 && rho > 0)
                                Gdum[l] = Gdum[l] + P[n].H2I*rho*fac;
                            if(vartype == 6 && rho > 0)
                                Gdum[l] = Gdum[l] + P[n].HDI*rho*fac;
                            if(vartype == 7 && rho > 0)
                                Gdum[l] = Gdum[l] + P[n].HII*rho*fac;
#if(READB)
                            if(vartype == 8 && rho > 0)
                                Gdum[l] = Gdum[l] + P[n].Bfieldx*rho*fac;
                            if(vartype == 9 && rho > 0)
                                Gdum[l] = Gdum[l] + P[n].Bfieldy*rho*fac;
                            if(vartype == 10 && rho > 0)
                                Gdum[l] = Gdum[l] + P[n].Bfieldz*rho*fac;
                            if(vartype == 11 && rho > 0)
                                Gdum[l] = Gdum[l] + P[n].nh_test*rho*fac;
#endif
                        }
                    }
            
        }
        
        printf("myrank = %d, main loop finished\n", myrank);
        
        //printf("MASSTEST: The amount of mass put into bins = %g\n",binnedMass);
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        if(vartype == 0)
        {
            //MPI_Allreduce(&Gnpart[0], &Gnpart_tot[0], send_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&Grho[0], &Grho_tot[0], send_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
        
        //////////write mapped data to file for Orion2!!!!!!!!!!!
        
        //write out density for each grid cell
        if(vartype == 0)
            MPI_Reduce(&Grho[0], &Gtot[0], send_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        
        if(vartype > 0)
            MPI_Reduce(&Gdum[0], &Gtot[0], send_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        
        printf("Rank 0 has the data now, myrank=%d\n", myrank);  fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
        
        if(myrank==0 && vartype==0) {
            double gridVol = pow(pcTOcm*width/ref_lev_doub,3);
            double binnedMass = 0;
            l=0;
            for(k=0;k<ref_lev;k++)
                for(j=0;j<ref_lev;j++)
                    for(i=0;i<ref_lev;i++)
                    {
                        binnedMass = binnedMass + Gtot[l]*gridVol;
                        l=l+1;
                    }
            printf("MASSTEST: The total binned mass m_grid = %g\n",binnedMass);
            printf("MASSTEST: The fraction m_gadget / m_grid = %g\n",InterestMtot/binnedMass);
        }
        
        if(myrank == 0)
        {
            outfile=fopen(outname,"a");
            masstot = 0;
            l=0;
            for(i=0;i<ref_lev;i++)
                for(j=0;j<ref_lev;j++)
                    for(k=0;k<ref_lev;k++)
                    {
                        if(vartype == 0)  masstot = masstot + Gtot[l]*pow(grid_size*pcTOcm, 3);
                        if(vartype > 0)  //write out stuff that is NOT density!
                            Gtot[l] = (Gtot[l]/Grho_tot[l]);
                        fwrite(&Gtot[l], sizeof(double), 1, outfile);
                        if(l % 1000 == 0)
                            printf("vartype = %d, Gtot[l] = %lg\n", vartype, Gtot[l]);
                        l++;
                    }
            if(vartype == 0) printf("masstot_alt = %lg\n", masstot/1.989e33);
            fclose(outfile);
        }
        
        if(myrank==0) for(k=0;k<ref_lev*ref_lev*ref_lev;k++) printf("Printing %g\n",Gtot[k]);
        
    }
    
    
    free(Gdum);
    free(Grho);
    free(Gtot);
    //free(Gnpart);
    //free(Gnpart_tot);
    free(Grho_tot);
    
    return(Ngas);
}





/* this routine allocates the memory for the
 * particle data.
 */
int allocate_memory(void)
{
    printf("allocating memory...\n");
    
    if(!(P=(struct particle_data *) malloc(NumPart*sizeof(struct particle_data))))
    {
        fprintf(stderr,"failed to allocate memory.\n");
        exit(0);
    }
    
    //  P--;   /* start with offset 1 */
    
    
    if(!(Id=(int *) malloc(NumPart*sizeof(int))))
    {
        fprintf(stderr,"failed to allocate memory.\n");
        exit(0);
    }
    
    //  Id--;   /* start with offset 1 */
    
    printf("allocating particle memory...done\n");
    return(0);
}



/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
int reordering(void)
{
    int i,j;
    int idsource, idsave, dest;
    struct particle_data psave, psource;
    
    
    printf("reordering....\n");
    
    for(i=1; i<=NumPart; i++)
    {
        if(Id[i] != i)
        {
            psource= P[i];
            idsource=Id[i];
            dest=Id[i];
            
            do
            {
                psave= P[dest];
                idsave=Id[dest];
                
                P[dest]= psource;
                Id[dest]= idsource;
                
                if(dest == i)
                    break;
                
                psource= psave;
                idsource=idsave;
                
                dest=idsource;
            }
            while(1);
        }
    }
    
    printf("done.\n");
    
    Id++;
    free(Id);
    
    printf("space for particle ID freed\n");
    return(0);
}


double calc_kernel_tsc(int n, double xgrid, double ygrid, double zgrid, double hsm, double grid_size, double grid_size_half)
{
    double rad, kernel, ratio, xratio, yratio, zratio, Wx, Wy, Wz, grid_size_simu;
    
    rad = (P[n].disx - xgrid)*(P[n].disx - xgrid) + (P[n].disy - ygrid)*(P[n].disy - ygrid)+ (P[n].disz - zgrid)*(P[n].disz - zgrid);
    rad = pow(rad,0.5);
    
    ratio = rad/P[n].hsm_phys;
    
    //triangular-shaped cloud kernel
    grid_size_simu = grid_size/1.e3/Time*(hubble_param);
    
    xratio = fabs((P[n].disx-xgrid)/grid_size);
    yratio = fabs((P[n].disy-ygrid)/grid_size);
    zratio = fabs((P[n].disz-zgrid)/grid_size);
    
    Wx = Wy = Wz = 0;
    
    if(ratio >= 1. && rad < grid_size_half)
        xratio = yratio = zratio = rad/(grid_size_half);
    
    if(xratio < 0.5)
        Wx = 0.75 - pow(xratio, 2);
    if(xratio < 1.5 && xratio >= 0.5)
        Wx = 0.5*pow(1.5 - xratio,2);
    
    if(yratio < 0.5)
        Wy = 0.75 - pow(yratio, 2);
    if(yratio < 1.5 && yratio >= 0.5)
        Wy = 0.5*pow(1.5 - yratio,2);
    
    if(zratio < 0.5)
        Wz = 0.75 - pow(zratio, 2);
    if(zratio < 1.5 && zratio >= 0.5)
        Wz = 0.5*pow(1.5 - zratio,2);
    
    kernel = pow(1./grid_size_simu,3)*Wx*Wy*Wz;
    
    return(kernel);
}


double calc_kernel_spline(int n, double xgrid, double ygrid, double zgrid, double hsm, double grid_size, double grid_size_half)
{
    double rad, kernel, ratio, xratio, yratio, zratio, Wx, Wy, Wz, grid_size_simu;
    double grid_size_max, fac=1.0;
    double radx, rady, radz;
    
    grid_size_max = pow(grid_size_half*grid_size_half + grid_size_half*grid_size_half + grid_size_half*grid_size_half,0.5);
    
    rad = (P[n].disx - xgrid)*(P[n].disx - xgrid) + (P[n].disy - ygrid)*(P[n].disy - ygrid)+ (P[n].disz - zgrid)*(P[n].disz - zgrid);
    rad = pow(rad,0.5);
    
    radx = fabs(P[n].disx - xgrid);
    rady = fabs(P[n].disy - ygrid);
    radz = fabs(P[n].disz - zgrid);
    
    ratio = rad/P[n].hsm_phys;
    
    
    if(ratio >= 1. && radx < grid_size_half && rady < grid_size_half && radz < grid_size_half)
    {
        ratio = rad/(grid_size_half);
        if(hsm < grid_size_half/1.e3/Time*(hubble_param))
            hsm = grid_size_half/1.e3/Time*(hubble_param);
        //printf("Dense particle1! nh = %lg, hsm = %lg\n", P[n].nh, P[n].hsm_phys);
    }
    
    
    //Gadget kernel
    
    /*
     if(ratio <= 0.5)
     kernel = fac*(8./3.14159/pow(hsm,3)) * (1. - 6.*pow(ratio,2) + 6.*pow(ratio,3));
     if(ratio > 0.5 && ratio <= 1.)
     kernel = (8./3.14159/pow(hsm,3)) * 2.*pow(1. - ratio, 3);
     if(ratio > 1.)
     kernel = 0.;
     */
    
    if(ratio < 1)
        kernel = (1./3.14159/pow(hsm,3)) * (1. - 1.5*pow(ratio,2) + 0.75*pow(ratio,3));
    if(ratio >= 1 && ratio <= 2)
        kernel = (1./3.14159/pow(hsm,3)) * (0.25*pow(2. - ratio, 3));
    if(ratio > 2)
        kernel = 0;
    
    
    //if(radx > grid_size_half || rady > grid_size_half || radz > grid_size_half)
    //   kernel = 0;
    
    
    return(kernel);
}




/* this template shows how one may convert from Gadget's units
 * to cgs units.
 * In this example, the temperate of the gas is computed.
 * (assuming that the electron density in units of the hydrogen density
 * was computed by the code. This is done if cooling is enabled.)
 */
int unit_conversion(void)
{
    double GRAVITY, BOLTZMANN, PROTONMASS;
    double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
    double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;
    double G, Xh, HubbleParam;
    
    int i;
    double MeanWeight, u, gamma;
    double h2frac, muh2, muh2in, temp;
    
    /* physical constants in cgs units */
    GRAVITY   = 6.672e-8;
    BOLTZMANN = 1.3806e-16;
    PROTONMASS = 1.6726e-24;
    
    /* internal unit system of the code */
    UnitLength_in_cm= 3.085678e21;   /*  code length unit in cm/h */
    UnitMass_in_g= 1.989e43;         /*  code mass unit in g/h */
    UnitVelocity_in_cm_per_s= 1.0e5;
    
    UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
    UnitDensity_in_cgs= UnitMass_in_g/ pow(UnitLength_in_cm,3);
    UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ pow(UnitTime_in_s,2);
    UnitEnergy_in_cgs= UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);
    
    G=GRAVITY/ pow(UnitLength_in_cm,3) * UnitMass_in_g * pow(UnitTime_in_s,2);
    
    
    Xh= 0.76e0;  /* mass fraction of hydrogen */
    HubbleParam= 0.7e0;
    
    //for NON-comoving sims
    //HubbleParam = 1.0;
    
    for(i=0; i<=NumPart; i++)
    {
        if(P[i].Type==0)  /* gas particle */
        {
            /*	  MeanWeight= 4.0/(3*Xh+1+4*Xh*P[i].elec) * PROTONMASS;  */
            
            MeanWeight=1.2195;
            h2frac=2.0*P[i].H2I;
            
            muh2in=(0.24/4.0) + ((1.0-h2frac)*0.76) + (h2frac*.76/2.0);
            muh2=pow(muh2in, -1.0);
            
            if(muh2 >= 1.22)
            {
                MeanWeight=muh2;
            }
            
            MeanWeight=MeanWeight*PROTONMASS;
            
            //MeanWeight= 1.22e0 * PROTONMASS;
            
            /* convert internal energy to cgs units */
            
            u  = P[i].U * UnitEnergy_in_cgs/ UnitMass_in_g;
            
            //gamma= 5.0/3.0;
            gamma = P[i].gam;
            
            /* get temperature in Kelvin */
            
            temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
            P[i].Rho= P[i].Density * UnitDensity_in_cgs * HubbleParam * HubbleParam * pow(1.e0+zred,3.e0);
            P[i].nh= P[i].Rho / MeanWeight;
            P[i].Pres = (gamma-1)*u*P[i].Rho;
            P[i].U = P[i].U * UnitEnergy_in_cgs/ UnitMass_in_g;
            
            // Just to be safe
            P[i].mapped_mass = 0.0;
            P[i].edge_flag = 0;
            
            /*  printf("zred = %g", zred);*/
        }
    }
    return(0);
}


/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int read_snapshot(char *fname, int files, char *outname, double delx, double dely, double delz, double  boxsize)
{
    
    FILE *fd;
    FILE *outfile;
    char   buf[200];
    int    i,j,k,l,dummy,ntot_withmasses;
    int    t,n,off,pc,pc_new=0,pc_sph;
    int NumPart_new = 0, Ngas_new = 0;
    int Idnew, nnew=1;
    double *pos, massnew, hsmnew, x, y, z, nnew_doub=1.0;
    double randomx, randomy, randomz;
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
    
    pos= (double*)malloc(sizeof(double) * 3);
    
    
    for(i=0, pc=0; i<files; i++, pc=pc_new)
    {
        if(files>1)
            sprintf(buf,"%s.%d",fname,i);
        else
            sprintf(buf,"%s",fname);
        
        if(!(fd=fopen(buf,"r")))
        {
            printf("can't open file `%s`\n",buf);
            exit(0);
        }
        
        printf("reading `%s' ...\n",buf); fflush(stdout);
        
        fread(&dummy, sizeof(dummy), 1, fd);
        fread(&header1, sizeof(header1), 1, fd);
        fread(&dummy, sizeof(dummy), 1, fd);
        
        if(files==1)
        {
            for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
            {
                NumPart+= header1.npart[k];
                printf("NumPart[%d] = %d\n", k, header1.npart[k]);
            }
            Ngas= header1.npart[0];
        }
        else
        {
            for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
                NumPart+= header1.npartTotal[k];
            Ngas= header1.npartTotal[0];
        }
        
        for(k=0, ntot_withmasses=0; k<5; k++)
        {
            if(header1.mass[k]==0)
                ntot_withmasses+= header1.npart[k];
        }
        
        if(i==0)
            allocate_memory();
        
        SKIP;
        for(k=0,pc_new=pc;k<6;k++)
        {
            for(n=0;n<header1.npart[k];n++)
            {
                fread(&P[pc_new].Pos[0], sizeof(double), 3, fd);
                pc_new++;
            }
        }
        SKIP;
        
        SKIP;
        for(k=0,pc_new=pc;k<6;k++)
        {
            for(n=0;n<header1.npart[k];n++)
            {
                fread(&P[pc_new].Vel[0], sizeof(double), 3, fd);
                pc_new++;
            }
        }
        SKIP;
        
        
        SKIP;
        for(k=0,pc_new=pc;k<6;k++)
        {
            for(n=0;n<header1.npart[k];n++)
            {
                fread(&Id[pc_new], sizeof(int), 1, fd);
                pc_new++;
            }
        }
        SKIP;
        
        if(ntot_withmasses>0)
        {
            SKIP;
        }
        
        for(k=0, pc_new=pc; k<6; k++)
        {
            for(n=0;n<header1.npart[k];n++)
            {
                P[pc_new].Type=k;
                if(header1.mass[k]==0)
                {
                    fread(&P[pc_new].Mass, sizeof(double), 1, fd);
                }
                else
                    P[pc_new].Mass= header1.mass[k];
                pc_new++;
            }
        }
        if(ntot_withmasses>0)
        {
            SKIP;
        }
        
        
        if(header1.npart[0]>0)
        {
            SKIP;
            for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
                fread(&P[pc_sph].U, sizeof(double), 1, fd);
                pc_sph++;
            }
            SKIP;
            
            SKIP;
            for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
                fread(&P[pc_sph].Density, sizeof(double), 1, fd);
                pc_sph++;
            }
            SKIP;
            
            
            SKIP;
            for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
                fread(&P[pc_sph].hsm, sizeof(double), 1, fd);
                pc_sph++;
            }
            SKIP;
            
            SKIP;
            for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
                fread(&P[pc_sph].H2I, sizeof(double), 1, fd);
                fread(&P[pc_sph].HII, sizeof(double), 1, fd);
                fread(&P[pc_sph].DII, sizeof(double), 1, fd);
                fread(&P[pc_sph].HDI, sizeof(double), 1, fd);
                fread(&P[pc_sph].HeII, sizeof(double), 1, fd);
                fread(&P[pc_sph].HeIII, sizeof(double), 1, fd);
                pc_sph++;
            }
            SKIP;
            
            SKIP;
            for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
                fread(&P[pc_sph].gam, sizeof(double), 1, fd);
                pc_sph++;
            }
            SKIP;
            
            printf("gam = %lg\n",P[100].gam);
            
            SKIP;
            for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
                fread(&P[pc_sph].sink, sizeof(double), 1, fd);
                pc_sph++;
            }
            SKIP;
            
        }
        fclose(fd);
    }
    
    Time= header1.time;
    zred= header1.redshift;
    printf("z= %6.2f \n",zred);
    printf("Time= %12.7e \n",Time);
    
    //For NON-cosmological runs
    //Time = 1.0;
    printf("Time= %12.7e \n",Time);
    
    printf("L= %6.2f \n",header1.BoxSize);
    fflush(stdout);
    return(Ngas);
}



