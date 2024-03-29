#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <gsl/gsl_spline.h>
#include <math.h>

#include "LALSimIMREOBNRv2.h"

#ifndef _LALSIMIMRSPINEOB_H
#define _LALSIMIMRSPINEOB_H

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif



/**
 * Set the total number of multipoles
 * */
#define MAX_NUM_MODES 7


struct
SpinEOBModes
{
  INT4 lmModes[MAX_NUM_MODES][2];
  REAL8TimeSeries *hlposm[MAX_NUM_MODES];
  REAL8TimeSeries *hlnegm[MAX_NUM_MODES];

  struct SpinEOBModes *next;
};

/**
 * Parameters for the spinning EOB model, used in calculating the Hamiltonian.
 * The Hamiltonian is given in Barausse and Buonanno (http://arxiv.org/pdf/0912.3517)
 * The parameters correspond to the following as found in the paper:
 * KK - K found in the equations for \f$\Delta_i\f$ (Eqn 5.77-5.81)
 * k0 - \f$\Delta_0\f$ given in Eqn 5.77
 * k1 - \f$\Delta_1\f$ given in Eqn 5.78
 * k2 - \f$\Delta_2\f$ given in Eqn 5.79
 * k3 - \f$\Delta_3\f$ given in Eqn 5.80
 * k4 - \f$\Delta_4\f$ given in Eqn 5.81
 * b3 - \f$\omega^{fd}_2\f$ given in Eqn 5.40
 * bb3 - \f$\omega^{fd}_1\f$ given in Eqn 5.40
 * d1 - SO calibration parameter of SEOBNRv1
 * d1v2 - SO calibration parameter of SEOBNRv2
 * dheffSS - SS calibration parameter of SEOBNRv1
 * dheffSSv2 - SS calibration parameter of SEOBNRv2
 * tidal1 - tidal params of body 1
 * tidal2 - tidal params of body 2
 */

typedef struct
tagSpinEOBHCoeffs
{
  double KK;
  double k0;
  double k1;
  double k2;
  double k3;
  double k4;
  double k5;
  double k5l;
  double b3;
  double bb3;
  double d1;
  double d1v2;
  double dheffSS;
  double dheffSSv2;
  UINT4    SpinAlignedEOBversion;
  int      updateHCoeffs;
  TidalEOBParams *tidal1;
  TidalEOBParams *tidal2;
}
SpinEOBHCoeffs;

typedef struct
tagSEOBHCoeffConstants
{

  double a0k2; //Coefficient of a^0 in k2
  double a1k2; //Coefficient of a^1 in k2

  double a0k3; //Coefficient of a^0 in k3
  double a1k3; //Coefficient of a^1 in k3

  double a0k4; //Coefficient of a^0 in k4
  double a1k4; //Coefficient of a^1 in k4
  double a2k4; //Coefficient of a^2 in k4

  double a0k5; //Coefficient of a^0 in k5
  double a1k5; //Coefficient of a^1 in k5
  double a2k5; //Coefficient of a^2 in k5

}
SEOBHCoeffConstants;

SEOBHCoeffConstants XLALEOBSpinPrecCalcSEOBHCoeffConstants(REAL8 eta);

/**
 * Parameters for the spinning EOB model.
 * 1) eobParams contains parameters common to nonspin and spin EOBNR models,
 * including mass ratio, masses, pre-computed coefficients for potential, flux and waveforms,
 * NQC coefficients and Newtonian multiple prefixes.
 * 2) seobCoeffs contans parameters for calculating the spin EOB Hamiltonian.
 * 3) s1Vec and s2Vec are individual spins, in unit of total mass.
 * 4) sigmaStar and sigmaKerr are effective spins of the test-particle and background.
 * 5) a is the spin value being used for test-particle limit spin terms.
 * 6) alignedSpins and tortoise are controling flags.
 */

typedef struct
tagSpinEOBParams
{
  EOBParams               *eobParams;
  SpinEOBHCoeffs          *seobCoeffs;
  SEOBHCoeffConstants     *seobCoeffConsts;
  EOBNonQCCoeffs          *nqcCoeffs;
  REAL8Vector             *s1Vec;
  REAL8Vector             *s2Vec;
  REAL8Vector             *sigmaStar;
  REAL8Vector             *sigmaKerr;
  REAL8                   a;
  REAL8                   chi1;
  REAL8                   chi2;
  REAL8                   prev_dr;
  int                     alignedSpins;
  Approximant             seobApproximant; /*OPTV3*/
  int                     tortoise;
  int                     ignoreflux;
  REAL8 deltaT;
  int                     use_hm;
}
SpinEOBParams;

/* We need to encapsulate the data for the GSL derivative function */
typedef
struct tagHcapDerivParams
{
   const REAL8   *values;
   SpinEOBParams *params;
   UINT4         varyParam;
}
HcapDerivParams;

/* We need to encapsulate the data for calculating spherical 2nd derivatives */
typedef
struct tagHcapSphDeriv2Params
{
  const REAL8     *sphValues;
  SpinEOBParams   *params;
  UINT4           varyParam1;
  UINT4           varyParam2;
  INT4            use_optimized;
}
HcapSphDeriv2Params;

/* We need to encapsulate the data for the GSL derivative function */
typedef
struct tagPrecEulerAnglesIntegration
{
   gsl_spline *alpha_spline;
   gsl_spline *beta_spline;
   gsl_interp_accel *alpha_acc;
   gsl_interp_accel *beta_acc;
}
PrecEulerAnglesIntegration;












#endif /* _LALSIMIMRSPINEOB_H */
