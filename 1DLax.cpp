// 1D Lax scheme code

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <time.h>
#include "mpi.h"
#include "input1Dhydro.h"
using namespace std;

// In preperation for setting error of TIMESTEP
#define UPPER_LIMIT_OF_TIMESTEP 1000

double u(double, double);         // velocity : u = u(rho, m)
double p(double, double, double); // pressure : p = p(rho, m, e)

void   outputData(int, double**, int);
void   FluxCalculation(double**, double**, int, int);
void   DataCommunication(double**, double**, int, int, int, int, int, int);
void   OutFlowBoundaryCondition(double**, double**, int);
void   /*Lax*/Scheme(double**, double**, double**);

const int    firstj        = 1;
const int    lastj         = (PART);
const int    SpaceArea     = (PART)+2;
const int    LeftBoundary  = 0;
const int    RightBoundary = (PART)+1;
const int    mesh          = (PART)*(CORE_NUM);
const double dx            = ((MAX)-(MIN)) / mesh;
const double dt            = (CFL_NUMBER)*dx / (ISO_SOUND_SPEED);
const double dt_dx         = dt/dx;

int main(int argc, char* argv[]) {

  int size, rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (size != CORE_NUM) {
    if (rank == 0) cerr << "====CORE_NUM_ERROR====" << endl;
    MPI_Finalize();
    return 1;
  }

  clock_t t_start, t_end;
  if (rank == 0) t_start = clock();

  double **Q1 = new double*[SpaceArea];
  double **Q2 = new double*[SpaceArea];
  double **E  = new double*[SpaceArea];

  for (int j=0; j<SpaceArea; j++) {
    Q1[j] = new double[3]; // [0]:rho, [1]:m, [2]:e
    Q2[j] = new double[3]; // [0]:rho, [1]:m, [2]:e
    E[j]  = new double[3]; // Each flux of element of Q
  }

  int tagsl = (10*rank) + 10;     // MPI tag for sending data to left
  int tagsr = (10*rank) + 10 + 1; // MPI tag for sending data to right
  int tagrl = (10*rank) + 1;      // MPI tag for receiving data from left
  int tagrr = (10*rank) + 20;     // MPI tag for receiving data from right
  int left  = rank - 1;
  int right = rank + 1;

  if (rank == 0)          left  = MPI_PROC_NULL;
  if (rank == CORE_NUM-1) right = MPI_PROC_NULL;

  // Initialization ( 0 time step )
  int timeStep = 0;
  if (rank == 0) cout << "Time step = " << timeStep << endl;
  InitialConditionInput(Q1, rank);
  outputData(timeStep, Q1, rank);

  // Calculation with a scheme ( 1 time step - )
  while (timeStep < UPPER_LIMIT_OF_TIMESTEP) {
    // Odd time step calculation with Q1, E ---> Q2
    timeStep++;
    if (rank == 0) cout << "Time step = " << timeStep << endl;

    FluxCalculation(Q1, E, firstj, lastj);
    DataCommunication(Q1, E, left, right, tagsl, tagsr, tagrl, tagrr);
    OutFlowBoundaryCondition(Q1, E, rank);
    Scheme(Q1, E, Q2);

    if (timeStep%OUTPUT_FREQ == 0) outputData(timeStep, Q2, rank);
    if (timeStep == TIMESTEP) break;

    // Even time step calculation with Q2, E ---> Q1
    timeStep++;
    if (rank == 0) cout << "Time step = " << timeStep << endl;

    FluxCalculation(Q2, E, firstj, lastj);
    DataCommunication(Q2, E, left, right, tagsl, tagsr, tagrl, tagrr);
    OutFlowBoundaryCondition(Q2, E, rank);
    Scheme(Q2, E, Q1);

    if (timeStep%OUTPUT_FREQ == 0) outputData(timeStep, Q1, rank);
    if (timeStep == TIMESTEP) break;
  }

  // Output result
  if (rank == 0) {
    t_end = clock();
    double duration = (double)(t_end - t_start) / CLOCKS_PER_SEC;
    cout << endl
         << " -------------------------------------------------" << endl
         << "   MESH            = " << mesh                      << endl
         << "   dx              = " << dx                        << endl
         << "   x range is [ " << MIN << " : " << MAX << " ]."   << endl
         << " -------------------------------------------------" << endl
         << "   TIMESTEP        = " << TIMESTEP                  << endl
         << "   dt              = " << dt                        << endl
         << "   t range is [ 0 : " << TIMESTEP*dt << " ]."       << endl
         << " -------------------------------------------------" << endl
         << "   GAMMA           = " << GAMMA                     << endl
         << "   ISO SOUND SPEED = " << ISO_SOUND_SPEED           << endl
         << "   CFL NUMBER      = " << CFL_NUMBER                << endl
         << " -------------------------------------------------" << endl
         << "   FILE NAME       = " << FILE_NAME                 << endl
         << " -------------------------------------------------" << endl
         << "   Duration        = " << duration << " s"          << endl
         << " -------------------------------------------------" << endl
         << endl;
  }

  for (int j=0; j<SpaceArea; j++) {
    delete[] Q1[j];
    delete[] Q2[j];
    delete[] E[j];
  }
  delete[] Q1;
  delete[] Q2;
  delete[] E;

  MPI_Finalize();
  return 0;

}

// Definition of velocity : u = u(rho, m)
double u(double rho, double m) { return m/rho; }

// Definition of pressure : p = p(rho, m, e)
double p(double rho, double m, double e) { return ( (GAMMA)-1.0 )*( e - 0.5*m*m/rho ); }

void outputData(int step, double** Q, int rank) {
  stringstream filename;
  filename << FILE_NAME << "_rank" << to_string(rank) << "_"
           << setw(4) << setfill('0') << step << ".tab";

  ofstream outputfile;
  outputfile.open( filename.str().c_str() );
  outputfile << "# rank = " << rank << endl
             << fixed << setprecision(7)
             << "# [t = " << step*dt << "] "
             << "(j x density pressure velocity) =" << endl;

  double minx = (MIN) + rank*(PART)*dx;
  for (int j=firstj; j<=lastj; j++) {
    outputfile << setw(4) << setfill('0') << j << "\t"
               << minx + (j-1)*dx              << "\t"
               << Q[j][0]                      << "\t"
               << p(Q[j][0], Q[j][1], Q[j][2]) << "\t"
               << u(Q[j][0], Q[j][1])          << endl;
  }

  outputfile.close();
}

void FluxCalculation(double** Q, double** E, int first, int last) {
  for (int j=first; j<=last; j++) {
    double frac_rho = 1.0 / Q[j][0];
    E[j][0]  = Q[j][1];
    E[j][1]  = ( (GAMMA)-1.0 )*Q[j][2];
    E[j][1] += 0.5*( 3.0-(GAMMA) )*Q[j][1]*Q[j][1]*frac_rho;
    E[j][2]  = (GAMMA)*Q[j][1]*Q[j][2]*frac_rho;
    E[j][2] -= 0.5*( (GAMMA)-1.0 )*Q[j][1]*Q[j][1]*Q[j][1]*frac_rho*frac_rho;
  }
}

void DataCommunication(double** Q, double** E, int left, int right,
                       int tagsl, int tagsr, int tagrl, int tagrr) {

  MPI_Bsend(Q[firstj], 3, MPI_DOUBLE, left,  tagsl, MPI_COMM_WORLD);
  MPI_Bsend(Q[lastj],  3, MPI_DOUBLE, right, tagsr, MPI_COMM_WORLD);
  MPI_Recv(Q[LeftBoundary],  3, MPI_DOUBLE, left,  tagrl, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(Q[RightBoundary], 3, MPI_DOUBLE, right, tagrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  MPI_Bsend(E[firstj], 3, MPI_DOUBLE, left,  tagsl, MPI_COMM_WORLD);
  MPI_Bsend(E[lastj],  3, MPI_DOUBLE, right, tagsr, MPI_COMM_WORLD);
  MPI_Recv(E[LeftBoundary],  3, MPI_DOUBLE, left,  tagrl, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(E[RightBoundary], 3, MPI_DOUBLE, right, tagrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void OutFlowBoundaryCondition(double** Q, double** E, int rank) {
  if (rank == 0) {
    for (int element=0; element<3; element++) {
      Q[LeftBoundary][element] = Q[firstj][element];
      E[LeftBoundary][element] = E[firstj][element];
    }
  }
  if (rank == (CORE_NUM)-1) {
    for (int element=0; element<3; element++) {
      Q[RightBoundary][element] = Q[lastj][element];
      E[RightBoundary][element] = E[lastj][element];
    }
  }
}

void /*Lax*/Scheme(double** Q, double** E, double** Qnext) {
  for (int j=firstj; j<=lastj; j++) {
  for (int element=0; element<3; element++) {
    Qnext[j][element]  = 0.5*( Q[j+1][element] + Q[j-1][element] );
    Qnext[j][element] -= 0.5*( E[j+1][element] - E[j-1][element] )*dt_dx;
  }}
}
