#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>

#define  PROTONMASS  1.6726e-24
#define  BOLTZMANN   1.3806e-16

#define MAX_STRING_LEN 1024
#define TABSIZE 91

double get_metaldependent_cooling_rate(double logTemp, double logZ);  // pass: log10(temperatue/Kelvin), log10(metallicity) 

static double CoolRate[8][TABSIZE];

static char *name[] = {
	"stripped_mzero.cie",
	"stripped_m-30.cie",
	"stripped_m-20.cie",
	"stripped_m-15.cie",
	"stripped_m-10.cie",
	"stripped_m-05.cie",
	"stripped_m-00.cie",
	"stripped_m+05.cie"
};


// Metallicies with repect to solar. Will be converted to absolut metallicities by adding log10(Z_sun), Zsun=0.02 
static double metallicities[8] = {
	-5.0,   // actually primordial -> -infinity 
	-3.0,
	-2.0,
	-1.5,
	-1.0,
	-0.5,
	+0.0,
	+0.5
};

static double CoolRate[8][TABSIZE];

void read_cooling_functions(void)
{
  FILE *fd;
  char buf[MAX_STRING_LEN];
  int i, n;
  float sd_logT, sd_ne, sd_nh, sd_nt, sd_logLnet,
    sd_logLnorm, sd_logU, sd_logTau, sd_logP12, sd_logRho24, sd_ci, sd_mubar;

  for(i = 0; i < 8; i++)
    metallicities[i] += log10(0.02);     // add solar metallicity 

  for(i = 0; i < 8; i++)
  {
    snprintf(buf, MAX_STRING_LEN - 1, "../../sage/extra/CoolFunctions/%s", name[i]);

    if(!(fd = fopen(buf, "r")))
    {
      printf("file `%s' not found\n", buf);
      exit(EXIT_FAILURE); 
    }

    for(n = 0; n <= 90; n++)
    {
      fscanf(fd, " %f %f %f %f %f %f %f %f %f %f %f %f ",
        &sd_logT, &sd_ne, &sd_nh, &sd_nt,
        &sd_logLnet, &sd_logLnorm, &sd_logU,
        &sd_logTau, &sd_logP12, &sd_logRho24, &sd_ci, &sd_mubar);

      CoolRate[i][n] = sd_logLnorm;
    }

    fclose(fd);
  }

  printf("cooling functions read\n\n");

}

int32_t main(int argc, char **argv)
{

#define UnitMass_in_g 1.989e+43
#define UnitLength_in_cm 3.08568e+24 
#define UnitVelocity_in_cm_per_s 100000

  double UnitDensity_in_cgs = UnitMass_in_g / pow(UnitLength_in_cm, 3);
  double UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;


  read_cooling_functions();
  double temp, logZ, lambda, x;

  double Tlow = 3.0;
  double Thigh = 7.0;
  double dT = 0.01;

  double logZlow = -10.0;
  double logZhigh = 1.0;
  double dlogZ = 0.01;

  double *results;

  int32_t N_T = (Thigh - Tlow) / dT;
  int32_t N_Z = (logZhigh - logZlow) / dlogZ;

  int32_t count = 0;

  printf("%d %d %d\n", N_T, N_Z, N_T*N_Z);
  //results = malloc(N_T * N_Z * 10 * sizeof(*(results)));

  FILE *filetowrite;
  char fname[MAX_STRING_LEN];

  snprintf(fname, MAX_STRING_LEN, "../../sage/extra/CoolFunctions/lookuptable");

  filetowrite = fopen(fname, "wb");
  if (filetowrite == NULL)
  {
    fprintf(stderr, "Could not open file %s\n", fname);
    return EXIT_FAILURE;
  }
 
  fwrite(&Tlow, sizeof(double), 1, filetowrite);
  fwrite(&Thigh, sizeof(double), 1, filetowrite);
  fwrite(&dT, sizeof(double), 1, filetowrite);

  fwrite(&logZlow, sizeof(double), 1, filetowrite);
  fwrite(&logZhigh, sizeof(double), 1, filetowrite);
  fwrite(&dlogZ, sizeof(double), 1, filetowrite);
 
  for (temp = Tlow; temp < Thigh; temp += dT)
  {
    for (logZ = logZlow; logZ < logZhigh; logZ += dlogZ)
    {    
      lambda = get_metaldependent_cooling_rate(temp, logZ);
      x = PROTONMASS * BOLTZMANN * pow(10, temp) / lambda; // now this has units sec g/cm^3  
      x /= (UnitDensity_in_cgs * UnitTime_in_s); // now in internal units 

      fwrite(&x, sizeof(double), 1, filetowrite);

      ++count;

    }

  }
 
  printf("Wrote %d elements\n", count); 
  fclose(filetowrite);


  return EXIT_SUCCESS;
}


double get_metaldependent_cooling_rate(double logTemp, double logZ)  // pass: log10(temperatue/Kelvin), log10(metallicity) 
{
  int i;
  double get_rate(int tab, double logTemp);
  double rate1, rate2, rate;


  if(logZ < metallicities[0])
    logZ = metallicities[0];

  if(logZ > metallicities[7])
    logZ = metallicities[7];

  i = 0;
  while(logZ > metallicities[i + 1])
    i++;

  // look up at i and i+1 
  rate1 = get_rate(i, logTemp);
  rate2 = get_rate(i + 1, logTemp);

  rate = rate1 + (rate2 - rate1) / (metallicities[i + 1] - metallicities[i]) * (logZ - metallicities[i]);

  return pow(10, rate);
}

double get_rate(int tab, double logTemp)
{
  int index;
  double rate1, rate2, rate, logTindex;

  if(logTemp < 4.0)
    logTemp = 4.0;

  index = (logTemp - 4.0) / 0.05;
  if(index >= 90)
    index = 89;

  logTindex = 4.0 + 0.05 * index;

  rate1 = CoolRate[tab][index];
  rate2 = CoolRate[tab][index + 1];

  rate = rate1 + (rate2 - rate1) / (0.05) * (logTemp - logTindex);

  return rate;
}

