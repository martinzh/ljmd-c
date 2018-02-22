/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

/* generic file- or pathname buffer length */
#define BLEN 200

/* a few physical constants */
const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */

/* structure to hold the complete information
 * about the MD system */
struct _mdsys {
    int natoms,nfi,nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
};

typedef struct _mdsys mdsys_t;

struct _params{
    int p_ints[3];
    double p_doubles[6];
};

typedef struct _params params_t;

/* helper function: read a line and then return
   the first string with whitespace stripped off */
static int get_a_line(FILE *fp, char *buf)
{
    char tmp[BLEN], *ptr;

    /* read a line and cut of comments and blanks */
    if (fgets(tmp,BLEN,fp)) {
        int i;

        ptr=strchr(tmp,'#');
        if (ptr) *ptr= '\0';
        i=strlen(tmp); --i;
        while(isspace(tmp[i])) {
            tmp[i]='\0';
            --i;
        }
        ptr=tmp;
        while(isspace(*ptr)) {++ptr;}
        i=strlen(ptr);
        strcpy(buf,tmp);
        return 0;
    } else {
        perror("problem reading input");
        return -1;
    }
    return 0;
}

/* helper function: zero out an array */
static void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i]=0.0;
    }
}

/* helper function: apply minimum image convention */
static double pbc(double x, const double boxby2)
{
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}

/* compute kinetic energy */
static void ekin(mdsys_t *sys)
{
    int i;

    sys->ekin=0.0;
    for (i=0; i<sys->natoms; ++i) {
        sys->ekin += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}

/* compute forces */
static void force(mdsys_t *sys)
{
    double r,ffac;
    double rx,ry,rz;
    int i,j;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    for(i=0; i < (sys->natoms); ++i) {
        for(j=0; j < (sys->natoms); ++j) {

            /* particles have no interactions with themselves */
            if (i==j) continue;

            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            r = sqrt(rx*rx + ry*ry + rz*rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
                ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
                                         +6*pow(sys->sigma/r,6.0)/r);

                sys->epot += 0.5*4.0*sys->epsilon*(pow(sys->sigma/r,12.0)
                                               -pow(sys->sigma/r,6.0));

                sys->fx[i] += rx/r*ffac;
                sys->fy[i] += ry/r*ffac;
                sys->fz[i] += rz/r*ffac;
            }
        }
    }
}

/* velocity verlet */
static void velverlet(mdsys_t *sys)
{
    int i;

    /* first part: propagate velocities by half and positions by full step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        sys->rx[i] += sys->dt*sys->vx[i];
        sys->ry[i] += sys->dt*sys->vy[i];
        sys->rz[i] += sys->dt*sys->vz[i];
    }

    /* compute forces and potential energy */
    force(sys);

    /* second part: propagate velocities by another half step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
    }
}

/* append data to output. */
static void output(mdsys_t *sys, FILE *erg, FILE *traj)
{
    int i;

    printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
    fprintf(erg,"% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
    fprintf(traj,"%d\n nfi=%d etot=%20.8f\n", sys->natoms, sys->nfi, sys->ekin+sys->epot);
    for (i=0; i<sys->natoms; ++i) {
        fprintf(traj, "Ar  %20.8f %20.8f %20.8f\n", sys->rx[i], sys->ry[i], sys->rz[i]);
    }
}


/* main */
int main(int argc, char **argv)
{
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;

    mdsys_t sys;
    params_t par_sys; // to send parameters to other pr

    int rank;
    int num_t;

    int chunk_size;

    int start_id;
    int temp_start_id;

    int end_id;

    MPI_Init(NULL, NULL);

    MPI_Request request;
    MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_t);

    if(rank == 0){
        /* read input file */
        if(get_a_line(stdin,line)) return 1;
        sys.natoms=atoi(line);
        if(get_a_line(stdin,line)) return 1;
        sys.mass=atof(line);
        if(get_a_line(stdin,line)) return 1;
        sys.epsilon=atof(line);
        if(get_a_line(stdin,line)) return 1;
        sys.sigma=atof(line);
        if(get_a_line(stdin,line)) return 1;
        sys.rcut=atof(line);
        if(get_a_line(stdin,line)) return 1;
        sys.box=atof(line);
        if(get_a_line(stdin,restfile)) return 1;
        if(get_a_line(stdin,trajfile)) return 1;
        if(get_a_line(stdin,ergfile)) return 1;
        if(get_a_line(stdin,line)) return 1;
        sys.nsteps=atoi(line);
        if(get_a_line(stdin,line)) return 1;
        sys.dt=atof(line);
        if(get_a_line(stdin,line)) return 1;
        nprint=atoi(line);

        par_sys.p_ints[0] = sys.natoms;
        par_sys.p_ints[1] = sys.nfi;
        par_sys.p_ints[2] = sys.nsteps;
        par_sys.p_doubles[0] = sys.dt;
        par_sys.p_doubles[1] = sys.mass;
        par_sys.p_doubles[2] = sys.epsilon;
        par_sys.p_doubles[3] = sys.sigma;
        par_sys.p_doubles[4] = sys.box;
        par_sys.p_doubles[5] = sys.rcut;

    }

    MPI_Datatype ParametersType;
    MPI_Datatype type[2] = {MPI_INT, MPI_DOUBLE};
    int          blocklen[2] = {3,6};
    MPI_Aint     displacements[2];

    MPI_Address(&par_sys.p_ints, &displacements[0]);
    MPI_Address(&par_sys.p_doubles, &displacements[1]);

    displacements[1] = displacements[1] - displacements[0];
    displacements[0] = 0;

    MPI_Type_struct( 2, blocklen, displacements, type, &ParametersType);
    MPI_Type_commit( &ParametersType);

    MPI_Bcast(&par_sys, 1, ParametersType, 0, MPI_COMM_WORLD);

    if( rank != 0){
        sys.natoms = par_sys.p_ints[0];
        sys.nfi = par_sys.p_ints[1];
        sys.nsteps =  par_sys.p_ints[2];
        sys.dt = par_sys.p_doubles[0];
        sys.mass =  par_sys.p_doubles[1];
        sys.epsilon = par_sys.p_doubles[2];
        sys.sigma = par_sys.p_doubles[3];
        sys.box = par_sys.p_doubles[4];
        sys.rcut = par_sys.p_doubles[5];
    }

    /* allocate memory */
   sys.rx=(double *)malloc(sys.natoms*sizeof(double));
   sys.ry=(double *)malloc(sys.natoms*sizeof(double));
   sys.rz=(double *)malloc(sys.natoms*sizeof(double));
   sys.vx=(double *)malloc(sys.natoms*sizeof(double));
   sys.vy=(double *)malloc(sys.natoms*sizeof(double));
   sys.vz=(double *)malloc(sys.natoms*sizeof(double));
   sys.fx=(double *)malloc(sys.natoms*sizeof(double));
   sys.fy=(double *)malloc(sys.natoms*sizeof(double));
   sys.fz=(double *)malloc(sys.natoms*sizeof(double));

    /* read restart */
    if(rank == 0){
        fp=fopen(restfile,"r");
        if(fp) {
           for (i=0; i<sys.natoms; ++i) {
               fscanf(fp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
           }
           for (i=0; i<sys.natoms; ++i) {
               fscanf(fp,"%lf%lf%lf",sys.vx+i, sys.vy+i, sys.vz+i);
           }
           fclose(fp);
        } else {
           perror("cannot read restart file");
           return 3;
        }

        erg=fopen(ergfile,"w");
        traj=fopen(trajfile,"w");

        printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
        printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
        output(&sys, erg, traj);
    }

    azzero(sys.fx, sys.natoms);
    azzero(sys.fy, sys.natoms);
    azzero(sys.fz, sys.natoms);

    MPI_Bcast(&sys.rx, 1, ParametersType, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.ry, 1, ParametersType, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.rz, 1, ParametersType, 0, MPI_COMM_WORLD);

    /* initialize forces and energies.*/
    sys.nfi=0;
    //force(&sys);
    //ekin(&sys);

    /**************************************************/

    chunk_size = sys.natoms / num_t;

    if(sys.natoms % num_t > rank) chunk_size++;

    if(rank == 0){
        start_id = 0;
        end_id = chunk_size;
        MPI_Isend(&end_id, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &request);
    }else{
        MPI_Irecv(&temp_start_id, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);

        start_id = temp_start_id + 1;

        end_id = temp_start_id + chunk_size;

        if(rank < num_t -1) MPI_Isend(&end_id, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &request);
    }

    // printf("rank %2d chunk_size = %4d total atoms = %6d\n", rank, chunk_size, sys.natoms);
    printf("rank %2d start_id = %4d end_id = %4d chunk = %4d\n", rank, start_id, end_id, end_id - start_id);


    /**************************************************/
    /* main MD loop */
    // for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {
    //
    //     /* write output, if requested */
    //     if ((sys.nfi % nprint) == 0)
    //         output(&sys, erg, traj);
    //
    //     /* propagate system and recompute energies */
    //     velverlet(&sys);
    //     ekin(&sys);
    // }
    /**************************************************/

    /* clean up: close files, free memory */
   if(rank == 0){
       printf("Simulation Done.\n");
       fclose(erg);
       fclose(traj);
   }

   free(sys.rx);
   free(sys.ry);
   free(sys.rz);
   free(sys.vx);
   free(sys.vy);
   free(sys.vz);
   free(sys.fx);
   free(sys.fy);
   free(sys.fz);

    MPI_Type_free(&ParametersType);
    MPI_Finalize();

    return 0;
}
