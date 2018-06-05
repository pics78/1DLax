/* Input data ============================================================ */
#define TIMESTEP             50
#define CORE_NUM              4
#define PART                 64
#define MIN                  -0.5
#define MAX                   0.5
#define GAMMA                 1.4
#define ISO_SOUND_SPEED       1.0
#define CFL_NUMBER            0.45
#define OUTPUT_FREQ           2
#define FILE_NAME          "ShockTube"
/* ======================================================================= */

void InitialConditionInput(double** Q0, int rank) {

  const int TotalMeshNum = (PART)*(CORE_NUM);
  double *input_den = new double[TotalMeshNum];
  double *input_pre = new double[TotalMeshNum];
  double *input_vel = new double[TotalMeshNum];

  /* Definition of initial condition */

  // shock tube condition
  for (int j=0; j<(int)(TotalMeshNum/2); j++) {
    input_den[j] = 1.0;
    input_pre[j] = 1.0;
    input_vel[j] = 0.0;
  }
  for (int j=(int)(TotalMeshNum/2); j<TotalMeshNum; j++) {
    input_den[j] = 0.125;
    input_pre[j] = 0.1;
    input_vel[j] = 0.0;
  }

  /* =============================== */

  int input_first = rank*(PART);
  int input_last  = input_first+(PART)-1;

  double frac_g = 1.0 / ( (GAMMA)-1.0 );

  int j = 1;
  for (int x=input_first; x<=input_last; x++) {
    Q0[j][0]  = input_den[x];
    Q0[j][1]  = input_den[x]*input_vel[x];
    Q0[j][2]  = frac_g*input_pre[x];
    Q0[j][2] += 0.5*input_den[x]*input_vel[x]*input_vel[x];
    j++;
  }

  delete[] input_den;
  delete[] input_pre;
  delete[] input_vel;
}
