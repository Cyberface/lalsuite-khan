/*
 * Copyright (C) 2015 Michael Puerrer, Sebastian Khan, Frank Ohme, Ofek Birnholtz, Lionel London, Francesco Pannarale
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

/**
 * \author Michael Puerrer, Sebastian Khan, Frank Ohme, Ofek Birnholtz, Lionel London, Francesco Pannarale
 *
 * \file
 *
 * \brief Internal function for IMRPhenomD phenomenological waveform model.
 * See \ref LALSimIMRPhenom_c for more details.
 *
 * Note that the functions defined here in LALSimIMRPhenomD_internals.c / .h
 * must be used by BOTH LALSimIMRPhenomD.c AND LALSimIMRPhenomP.c.
 */

/*
This waveform uses the TaylorF2 coefficients for it's inspiral phase augmented
by higher order phenomenological terms tuned to SEOBv2-Hybrid waveforms.
Below are lines copied from LALSimInspiralPNCoefficients.c which are the TaylorF2
phase coefficients we have used.
We document them here in case changes to that file changes the behaviour
of this waveform.

    const REAL8 mtot = m1 + m2;
    const REAL8 d = (m1 - m2) / (m1 + m2);
    const REAL8 eta = m1*m2/mtot/mtot;
    const REAL8 m1M = m1/mtot;
    const REAL8 m2M = m2/mtot;
    // Use the spin-orbit variables from arXiv:1303.7412, Eq. 3.9
    // We write dSigmaL for their (\delta m/m) * \Sigma_\ell
    // There's a division by mtotal^2 in both the energy and flux terms
    // We just absorb the division by mtotal^2 into SL and dSigmaL

    const REAL8 SL = m1M*m1M*chi1L + m2M*m2M*chi2L;
    const REAL8 dSigmaL = d*(m2M*chi2L - m1M*chi1L);

    const REAL8 pfaN = 3.L/(128.L * eta);
    //Non-spin phasing terms - see arXiv:0907.0700, Eq. 3.18
    pfa->v[0] = 1.L;
    pfa->v[2] = 5.L*(743.L/84.L + 11.L * eta)/9.L;
    pfa->v[3] = -16.L*LAL_PI;
    pfa->v[4] = 5.L*(3058.673L/7.056L + 5429.L/7.L * eta
                     + 617.L * eta*eta)/72.L;
    pfa->v[5] = 5.L/9.L * (7729.L/84.L - 13.L * eta) * LAL_PI;
    pfa->vlogv[5] = 5.L/3.L * (7729.L/84.L - 13.L * eta) * LAL_PI;
    pfa->v[6] = (11583.231236531L/4.694215680L
                     - 640.L/3.L * LAL_PI * LAL_PI - 6848.L/21.L*LAL_GAMMA)
                 + eta * (-15737.765635L/3.048192L
                     + 2255./12. * LAL_PI * LAL_PI)
                 + eta*eta * 76055.L/1728.L
                 - eta*eta*eta * 127825.L/1296.L;
    pfa->v[6] += (-6848.L/21.L)*log(4.);
    pfa->vlogv[6] = -6848.L/21.L;
    pfa->v[7] = LAL_PI * ( 77096675.L/254016.L
                     + 378515.L/1512.L * eta - 74045.L/756.L * eta*eta);

    // Spin-orbit terms - can be derived from arXiv:1303.7412, Eq. 3.15-16
    const REAL8 pn_gamma = (554345.L/1134.L + 110.L*eta/9.L)*SL + (13915.L/84.L - 10.L*eta/3.)*dSigmaL;
    switch( spinO )
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
            pfa->v[7] += (-8980424995.L/762048.L + 6586595.L*eta/756.L - 305.L*eta*eta/36.L)*SL - (170978035.L/48384.L - 2876425.L*eta/672.L - 4735.L*eta*eta/144.L) * dSigmaL;
*/

#include "LALSimIMRPhenomInternalUtils.h"
#include "LALSimIMRPhenomD_internals.h"

/*
 *
 * Internal function implementations
 *
 */

////////////////////////////// Miscellaneous functions //////////////////////////////

/**
 * PN reduced spin parameter
 * See Eq 5.9 in http://arxiv.org/pdf/1107.1267v2.pdf
 */
static double chiPN(double Seta, double eta, double chi1, double chi2)
{
  // Convention m1 >= m2 and chi1 is the spin on m1
  // The 0.5 factor missing in the definitions of chi_s and chi_a is
  // recovered in the return expresion
  double chi_s = (chi1 + chi2);
  double chi_a = (chi1 - chi2);

  return 0.5 * (chi_s * (1.0 - eta * 76.0 / 113.0) + Seta * chi_a);
}

/**
 * Return the closest higher power of 2
 */
static size_t NextPow2(const size_t n)
{
  // use pow here, not bit-wise shift, as the latter seems to run against an upper cutoff long before SIZE_MAX, at least on some platforms
  return (size_t) pow(2,ceil(log2(n)));
}

///**
// * Step function
// */
//static double StepFunc(const double t, const double t1) {
//  if (t < t1)
//    return 0.0;
//  else
//    return 1.0;
//}

/**
 * Step function in boolean version
 */
static bool StepFunc_boolean(const double t, const double t1) {
  return (!(t < t1));
}

//////////////////////// Final spin, final mass, fring, fdamp ////////////////////////

// Final Spin and Radiated Energy formulas described in 1508.07250

/**
 * Formula to predict the final spin. Equation 3.6 arXiv:1508.07250
 * s defined around Equation 3.6.
 */
static double FinalSpin0815_s(double eta, double s) {
  double eta2 = eta*eta;
  double eta3 = eta2*eta;
  double s2 = s*s;
  double s3 = s2*s;

/* FIXME: there are quite a few int's withouth a . in this file */
//FP: eta2, eta3 can be avoided
return eta*(3.4641016151377544 - 4.399247300629289*eta +
      9.397292189321194*eta2 - 13.180949901606242*eta3 +
      s*((1.0/eta - 0.0850917821418767 - 5.837029316602263*eta) +
      (0.1014665242971878 - 2.0967746996832157*eta)*s +
      (-1.3546806617824356 + 4.108962025369336*eta)*s2 +
      (-0.8676969352555539 + 2.064046835273906*eta)*s3));
}

/**
 * Wrapper function for FinalSpin0815_s.
 */
static double FinalSpin0815(double eta, double chi1, double chi2) {
  // Convention m1 >= m2
  double Seta = sqrt(1.0 - 4.0*eta);
  double m1 = 0.5 * (1.0 + Seta);
  double m2 = 0.5 * (1.0 - Seta);
  double m1s = m1*m1;
  double m2s = m2*m2;
  // s defined around Equation 3.6 arXiv:1508.07250
  double s = (m1s * chi1 + m2s * chi2);
  return FinalSpin0815_s(eta, s);
}

/**
 * Formula to predict the total radiated energy. Equation 3.7 and 3.8 arXiv:1508.07250
 * Input parameter s defined around Equation 3.7 and 3.8.
 */
static double EradRational0815_s(double eta, double s) {
  double eta2 = eta*eta;
  double eta3 = eta2*eta;

  return (eta*(0.055974469826360077 + 0.5809510763115132*eta - 0.9606726679372312*eta2 + 3.352411249771192*eta3)*
    (1. + (-0.0030302335878845507 - 2.0066110851351073*eta + 7.7050567802399215*eta2)*s))/(1. + (-0.6714403054720589 - 1.4756929437702908*eta + 7.304676214885011*eta2)*s);
}

/**
 * Wrapper function for EradRational0815_s.
 */
static double EradRational0815(double eta, double chi1, double chi2) {
  // Convention m1 >= m2
  double Seta = sqrt(1.0 - 4.0*eta);
  double m1 = 0.5 * (1.0 + Seta);
  double m2 = 0.5 * (1.0 - Seta);
  double m1s = m1*m1;
  double m2s = m2*m2;
  // arXiv:1508.07250
  double s = (m1s * chi1 + m2s * chi2) / (m1s + m2s);

  return EradRational0815_s(eta, s);
}

/**
 * fring is the real part of the ringdown frequency
 * 1508.07250 figure 9
 */
static double fring(double eta, double chi1, double chi2, double finspin) {
  double return_val;

  if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *iFring = gsl_spline_alloc(gsl_interp_cspline, QNMData_length);
  gsl_spline_init(iFring, QNMData_a, QNMData_fring, QNMData_length);

  return_val = gsl_spline_eval(iFring, finspin, acc) / (1.0 - EradRational0815(eta, chi1, chi2));

  gsl_spline_free(iFring);
  gsl_interp_accel_free(acc);
  return return_val;
}

/**
 * fdamp is the complex part of the ringdown frequency
 * 1508.07250 figure 9
 */
static double fdamp(double eta, double chi1, double chi2, double finspin) {
  double return_val;

  if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fdamp function: final spin > 1.0 not supported\n");

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *iFdamp = gsl_spline_alloc(gsl_interp_cspline, QNMData_length);
  gsl_spline_init(iFdamp, QNMData_a, QNMData_fdamp, QNMData_length);

  return_val = gsl_spline_eval(iFdamp, finspin, acc) / (1.0 - EradRational0815(eta, chi1, chi2));

  gsl_spline_free(iFdamp);
  gsl_interp_accel_free(acc);
  return return_val;
}

static int init_useful_powers(UsefulPowers *p, REAL8 number)
{
  XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");
  XLAL_CHECK(number >= 0, XLAL_EDOM, "number must be non-negative");

  // consider changing pow(x,1/6.0) to cbrt(x) and sqrt(x) - might be faster
  double sixth = pow(number, 1.0 / 6.0);
  p->third = sixth * sixth;
  //p->third = cbrt(number);
  p->two_thirds = p->third * p->third;
  p->four_thirds = number * p->third;
  p->five_thirds = p->four_thirds * p->third;
  p->two = number * number;
  p->seven_thirds = p->third * p->two;
  p->eight_thirds = p->two_thirds * p->two;
  p->inv = 1. / number;
  double m_sixth = 1.0 / sixth;
  p->m_seven_sixths = p->inv * m_sixth;
  p->m_third = m_sixth * m_sixth;
  p->m_two_thirds = p->m_third * p->m_third;
  p->m_five_thirds = p->inv * p->m_two_thirds;

  return XLAL_SUCCESS;
}

/******************************* Amplitude functions *******************************/

/**
 * amplitude scaling factor defined by eq. 17 in 1508.07253
 */
static double amp0Func(double eta) {
  return sqrt(2.0/3.0*eta)*PI_M_SIXTH;
}

///////////////////////////// Amplitude: Inspiral functions /////////////////////////

// Phenom coefficients rho1, ..., rho3 from direct fit
// AmpInsDFFitCoeffChiPNFunc[eta, chiPN]

/**
 * rho_1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double rho1_fun(double eta, double eta2, double xi) {

  return 3931.8979897196696 - 17395.758706812805*eta
  + (3132.375545898835 + 343965.86092361377*eta - 1.2162565819981997e6*eta2
  + (-70698.00600428853 + 1.383907177859705e6*eta - 3.9662761890979446e6*eta2)*xi
  + (-60017.52423652596 + 803515.1181825735*eta - 2.091710365941658e6*eta2)*xi*xi)*xi;
}


/**
 * rho_2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double rho2_fun(double eta, double eta2, double xi) {

  return -40105.47653771657 + 112253.0169706701*eta
  + (23561.696065836168 - 3.476180699403351e6*eta + 1.137593670849482e7*eta2
  + (754313.1127166454 - 1.308476044625268e7*eta + 3.6444584853928134e7*eta2)*xi
  + (596226.612472288 - 7.4277901143564405e6*eta + 1.8928977514040343e7*eta2)*xi*xi)*xi;
}

/**
 * rho_3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double rho3_fun(double eta, double eta2, double xi) {

  return 83208.35471266537 - 191237.7264145924*eta
  + (-210916.2454782992 + 8.71797508352568e6*eta - 2.6914942420669552e7*eta2
  + (-1.9889806527362722e6 + 3.0888029960154563e7*eta - 8.390870279256162e7*eta2)*xi
  + (-1.4535031953446497e6 + 1.7063528990822166e7*eta - 4.2748659731120914e7*eta2)*xi*xi)*xi;
}

// The Newtonian term in LAL is fine and we should use exactly the same (either hardcoded or call).
// We just use the Mathematica expression for convenience.
/**
 * Inspiral amplitude plus rho phenom coefficents. rho coefficients computed
 * in rho1_fun, rho2_fun, rho3_fun functions.
 * Amplitude is a re-expansion. See 1508.07253 and Equation 29, 30 and Appendix B arXiv:1508.07253 for details
 */
static double AmpInsAnsatz(double Mf, UsefulPowers * powers_of_Mf, AmpInsPrefactors * prefactors) {

  return 1 + powers_of_Mf->two_thirds * prefactors->two_thirds
      + powers_of_Mf->four_thirds * prefactors->four_thirds
			+ powers_of_Mf->five_thirds * prefactors->five_thirds
			+ powers_of_Mf->seven_thirds * prefactors->seven_thirds + powers_of_Mf->eight_thirds * prefactors->eight_thirds
			+ Mf * (prefactors->one + Mf * prefactors->two + powers_of_Mf->two * prefactors->three);
}

static int init_amp_ins_prefactors(AmpInsPrefactors * prefactors, IMRPhenomDAmplitudeCoefficients* p)
{
	XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");
	XLAL_CHECK(0 != prefactors, XLAL_EFAULT, "prefactors is NULL");

	double eta = p->eta;

	prefactors->amp0 = amp0Func(eta);

	double chi1 = p->chi1;
	double chi2 = p->chi2;

	double chi12 = p->chi12;
	double chi22 = p->chi22;
	double eta2 = p->eta2;
	double eta3 = p->eta3;

	double Pi = LAL_PI;
	double Pi2 = powers_of_pi.two;
	double Seta = p->Seta;
  double SetaPlus1 = p->SetaPlus1;

	prefactors->two_thirds = ((-969 + 1804*eta)*powers_of_pi.two_thirds)/672.;
	prefactors->one = ((chi1*(81*SetaPlus1 - 44*eta) + chi2*(81 - 81*Seta - 44*eta))*Pi)/48.;
	prefactors->four_thirds = (	(-27312085.0 - 10287648*chi22 - 10287648*chi12*SetaPlus1 + 10287648*chi22*Seta
								 + 24*(-1975055 + 857304*chi12 - 994896*chi1*chi2 + 857304*chi22)*eta
								 + 35371056*eta2
								 )
							* powers_of_pi.four_thirds) / 8.128512e6;
	prefactors->five_thirds = (powers_of_pi.five_thirds * (chi2*(-285197*(-1 + Seta) + 4*(-91902 + 1579*Seta)*eta - 35632*eta2)
															+ chi1*(285197*SetaPlus1 - 4*(91902 + 1579*Seta)*eta - 35632*eta2)
															+ 42840*(-1.0 + 4*eta)*Pi
															)
								) / 32256.;
	prefactors->two = - (Pi2*(-336*(-3248849057.0 + 2943675504*chi12 - 3339284256*chi1*chi2 + 2943675504*chi22)*eta2
							  - 324322727232*eta3
							  - 7*(-177520268561 + 107414046432*chi22 + 107414046432*chi12*SetaPlus1
									- 107414046432*chi22*Seta + 11087290368*(chi1 + chi2 + chi1*Seta - chi2*Seta)*Pi
									)
							  + 12*eta*(-545384828789 - 176491177632*chi1*chi2 + 202603761360*chi22
										+ 77616*chi12*(2610335 + 995766*Seta) - 77287373856*chi22*Seta
										+ 5841690624*(chi1 + chi2)*Pi + 21384760320*Pi2
										)
								)
						)/6.0085960704e10;
	prefactors->seven_thirds = p->rho1;
	prefactors->eight_thirds = p->rho2;
	prefactors->three = p->rho3;

	return XLAL_SUCCESS;
}

/**
 * Take the AmpInsAnsatz expression and compute the first derivative
 * with respect to frequency to get the expression below.
 */
static double DAmpInsAnsatz(double Mf, UsefulPowers *powers_of_Mf, IMRPhenomDAmplitudeCoefficients* p) {
  double eta = p->eta;
  double chi1 = p->chi1;
  double chi2 = p->chi2;
  //double rho1 = p->rho1;
  //double rho2 = p->rho2;
  //double rho3 = p->rho3;

  double chi12 = p->chi12;
  double chi22 = p->chi22;
  double eta2 = p->eta2;
  double eta3 = p->eta3;
  double Pi = LAL_PI;
  double Pi2 = powers_of_pi.two;
  double Seta = p->Seta;
  double SetaPlus1 = p->SetaPlus1;

   return ((-969 + 1804*eta)*powers_of_pi.two_thirds)/(1008.*powers_of_Mf->third)
   + ((chi1*(81*SetaPlus1 - 44*eta) + chi2*(81 - 81*Seta - 44*eta))*Pi)/48.
   + ((-27312085 - 10287648*chi22 - 10287648*chi12*SetaPlus1
   + 10287648*chi22*Seta + 24*(-1975055 + 857304*chi12 - 994896*chi1*chi2 + 857304*chi22)*eta
   + 35371056*eta2)*powers_of_Mf->third*powers_of_pi.four_thirds)/6.096384e6
   + (5*powers_of_Mf->two_thirds*powers_of_pi.five_thirds*(chi2*(-285197*(-1 + Seta)
   + 4*(-91902 + 1579*Seta)*eta - 35632*eta2) + chi1*(285197*SetaPlus1
   - 4*(91902 + 1579*Seta)*eta - 35632*eta2) + 42840*(-1 + 4*eta)*Pi))/96768.
   - (Mf*Pi2*(-336*(-3248849057.0 + 2943675504*chi12 - 3339284256*chi1*chi2 + 2943675504*chi22)*eta2 - 324322727232*eta3
   - 7*(-177520268561 + 107414046432*chi22 + 107414046432*chi12*SetaPlus1 - 107414046432*chi22*Seta
   + 11087290368*(chi1 + chi2 + chi1*Seta - chi2*Seta)*Pi)
   + 12*eta*(-545384828789.0 - 176491177632*chi1*chi2 + 202603761360*chi22 + 77616*chi12*(2610335 + 995766*Seta)
   - 77287373856*chi22*Seta + 5841690624*(chi1 + chi2)*Pi + 21384760320*Pi2)))/3.0042980352e10
   + (7.0/3.0)*powers_of_Mf->four_thirds*p->rho1 + (8.0/3.0)*powers_of_Mf->five_thirds*p->rho2 + 3.*powers_of_Mf->two*p->rho3;
}

/////////////////////////// Amplitude: Merger-Ringdown functions ///////////////////////

// Phenom coefficients gamma1, ..., gamma3
// AmpMRDAnsatzFunc[]

/**
 * gamma 1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double gamma1_fun(double eta, double eta2, double xi) {

  return 0.006927402739328343 + 0.03020474290328911*eta
  + (0.006308024337706171 - 0.12074130661131138*eta + 0.26271598905781324*eta2
  + (0.0034151773647198794 - 0.10779338611188374*eta + 0.27098966966891747*eta2)*xi
  + (0.0007374185938559283 - 0.02749621038376281*eta + 0.0733150789135702*eta2)*xi*xi)*xi;
}

/**
 * gamma 2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double gamma2_fun(double eta, double eta2, double xi) {

  return 1.010344404799477 + 0.0008993122007234548*eta
  + (0.283949116804459 - 4.049752962958005*eta + 13.207828172665366*eta2
  + (0.10396278486805426 - 7.025059158961947*eta + 24.784892370130475*eta2)*xi
  + (0.03093202475605892 - 2.6924023896851663*eta + 9.609374464684983*eta2)*xi*xi)*xi;
}

/**
 * gamma 3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double gamma3_fun(double eta, double eta2, double xi) {

  return 1.3081615607036106 - 0.005537729694807678*eta
  + (-0.06782917938621007 - 0.6689834970767117*eta + 3.403147966134083*eta2
  + (-0.05296577374411866 - 0.9923793203111362*eta + 4.820681208409587*eta2)*xi
  + (-0.006134139870393713 - 0.38429253308696365*eta + 1.7561754421985984*eta2)*xi*xi)*xi;
}

/**
 * Ansatz for the merger-ringdown amplitude. Equation 19 arXiv:1508.07253
 */
static double AmpMRDAnsatz(double f, IMRPhenomDAmplitudeCoefficients* p) {
  double fRD = p->fRD;
  double fDM = p->fDM;
  double gamma1 = p->gamma1;
  double gamma2 = p->gamma2;
  double gamma3 = p->gamma3;
  double fDMgamma3 = fDM*gamma3;
  double fminfRD = f - fRD;
  return exp( -(fminfRD)*gamma2 / (fDMgamma3) )
    * (fDMgamma3*gamma1) / (pow_2_of(fminfRD) + pow_2_of(fDMgamma3));
}

/**
 * first frequency derivative of AmpMRDAnsatz
 */
static double DAmpMRDAnsatz(double f, IMRPhenomDAmplitudeCoefficients* p) {
  double fRD = p->fRD;
  double fDM = p->fDM;
  double gamma1 = p->gamma1;
  double gamma2 = p->gamma2;
  double gamma3 = p->gamma3;

  double fDMgamma3 = fDM * gamma3;
  double pow2_fDMgamma3 = pow_2_of(fDMgamma3);
  double fminfRD = f - fRD;
  double expfactor = exp(((fminfRD)*gamma2)/(fDMgamma3));
  double pow2pluspow2 = pow_2_of(fminfRD) + pow2_fDMgamma3;

  return ((-2*fDM*(fminfRD)*gamma3*gamma1) / pow2pluspow2 -
    (gamma2*gamma1)) / ( expfactor * (pow2pluspow2)) ;
}

/**
 * Equation 20 arXiv:1508.07253 (called f_peak in paper)
 * analytic location of maximum of AmpMRDAnsatz
 */
static double fmaxCalc(IMRPhenomDAmplitudeCoefficients* p) {
  double fRD = p->fRD;
  double fDM = p->fDM;
  double gamma2 = p->gamma2;
  double gamma3 = p->gamma3;

  // NOTE: There's a problem with this expression from the paper becoming imaginary if gamma2>=1
  // Fix: if gamma2 >= 1 then set the square root term to zero.
  if (!(gamma2 > 1))
    return fabs(fRD + (fDM*(-1 + sqrt(1 - pow_2_of(gamma2)))*gamma3)/gamma2);
  else
    return fabs(fRD + (-fDM*gamma3)/gamma2);
}

///////////////////////////// Amplitude: Intermediate functions ////////////////////////

// Phenom coefficients delta0, ..., delta4 determined from collocation method
// (constraining 3 values and 2 derivatives)
// AmpIntAnsatzFunc[]

/**
 * Ansatz for the intermediate amplitude. Equation 21 arXiv:1508.07253
 */
static double AmpIntAnsatz(double Mf, IMRPhenomDAmplitudeCoefficients* p) {
  double Mf2 = Mf*Mf;
  return p->delta0 + Mf*p->delta1 + Mf2*(p->delta2 + Mf*p->delta3 + p->delta4*Mf2);
}

/**
 * The function name stands for 'Amplitude Intermediate Collocation Fit Coefficient'
 * This is the 'v2' value in Table 5 of arXiv:1508.07253
 */
static double AmpIntColFitCoeff(double eta, double eta2, double chi) {
  double xi = -1.0 + chi;

  return 0.8149838730507785 + 2.5747553517454658*eta
  + (1.1610198035496786 - 2.3627771785551537*eta + 6.771038707057573*eta2
  + (0.7570782938606834 - 2.7256896890432474*eta + 7.1140380397149965*eta2)*xi
  + (0.1766934149293479 - 0.7978690983168183*eta + 2.1162391502005153*eta2)*xi*xi)*xi;
}

  /**
  * The following functions (delta{0,1,2,3,4}_fun) were derived
  * in mathematica according to
  * the constraints detailed in arXiv:1508.07253,
  * section 'Region IIa - intermediate'.
  * These are not given in the paper.
  * Can be rederived by solving Equation 21 for the constraints
  * given in Equations 22-26 in arXiv:1508.07253
  */
static double delta0_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d) {
  double f1 = p->f1;
  double f2 = p->f2;
  double f3 = p->f3;
  double v1 = p->v1;
  double v2 = p->v2;
  double v3 = p->v3;
  double d1 = p->d1;
  double d2 = p->d2;

  double f12 = d->f12;
  double f13 = d->f13;
  double f14 = d->f14;
  double f15 = d->f15;
  double f22 = d->f22;
  double f23 = d->f23;
  double f24 = d->f24;
  double f32 = d->f32;
  double f33 = d->f33;
  double f34 = d->f34;
  double f35 = d->f35;

  return -((d2*f15*f22*f3 - 2*d2*f14*f23*f3 + d2*f13*f24*f3 - d2*f15*f2*f32 + d2*f14*f22*f32
  - d1*f13*f23*f32 + d2*f13*f23*f32 + d1*f12*f24*f32 - d2*f12*f24*f32 + d2*f14*f2*f33
  + 2*d1*f13*f22*f33 - 2*d2*f13*f22*f33 - d1*f12*f23*f33 + d2*f12*f23*f33 - d1*f1*f24*f33
  - d1*f13*f2*f34 - d1*f12*f22*f34 + 2*d1*f1*f23*f34 + d1*f12*f2*f35 - d1*f1*f22*f35
  + 4*f12*f23*f32*v1 - 3*f1*f24*f32*v1 - 8*f12*f22*f33*v1 + 4*f1*f23*f33*v1 + f24*f33*v1
  + 4*f12*f2*f34*v1 + f1*f22*f34*v1 - 2*f23*f34*v1 - 2*f1*f2*f35*v1 + f22*f35*v1 - f15*f32*v2
  + 3*f14*f33*v2 - 3*f13*f34*v2 + f12*f35*v2 - f15*f22*v3 + 2*f14*f23*v3 - f13*f24*v3
  + 2*f15*f2*f3*v3 - f14*f22*f3*v3 - 4*f13*f23*f3*v3 + 3*f12*f24*f3*v3 - 4*f14*f2*f32*v3
  + 8*f13*f22*f32*v3 - 4*f12*f23*f32*v3) / (pow_2_of(f1 - f2)*pow_3_of(f1 - f3)*pow_2_of(f3-f2)));
}

static double delta1_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d) {
  double f1 = p->f1;
  double f2 = p->f2;
  double f3 = p->f3;
  double v1 = p->v1;
  double v2 = p->v2;
  double v3 = p->v3;
  double d1 = p->d1;
  double d2 = p->d2;

  double f12 = d->f12;
  double f13 = d->f13;
  double f14 = d->f14;
  double f15 = d->f15;
  double f22 = d->f22;
  double f23 = d->f23;
  double f24 = d->f24;
  double f32 = d->f32;
  double f33 = d->f33;
  double f34 = d->f34;
  double f35 = d->f35;

  return -((-(d2*f15*f22) + 2*d2*f14*f23 - d2*f13*f24 - d2*f14*f22*f3 + 2*d1*f13*f23*f3
  + 2*d2*f13*f23*f3 - 2*d1*f12*f24*f3 - d2*f12*f24*f3 + d2*f15*f32 - 3*d1*f13*f22*f32
  - d2*f13*f22*f32 + 2*d1*f12*f23*f32 - 2*d2*f12*f23*f32 + d1*f1*f24*f32 + 2*d2*f1*f24*f32
  - d2*f14*f33 + d1*f12*f22*f33 + 3*d2*f12*f22*f33 - 2*d1*f1*f23*f33 - 2*d2*f1*f23*f33
  + d1*f24*f33 + d1*f13*f34 + d1*f1*f22*f34 - 2*d1*f23*f34 - d1*f12*f35 + d1*f22*f35
  - 8*f12*f23*f3*v1 + 6*f1*f24*f3*v1 + 12*f12*f22*f32*v1 - 8*f1*f23*f32*v1 - 4*f12*f34*v1
  + 2*f1*f35*v1 + 2*f15*f3*v2 - 4*f14*f32*v2 + 4*f12*f34*v2 - 2*f1*f35*v2 - 2*f15*f3*v3
  + 8*f12*f23*f3*v3 - 6*f1*f24*f3*v3 + 4*f14*f32*v3 - 12*f12*f22*f32*v3 + 8*f1*f23*f32*v3)
  / (pow_2_of(f1 - f2)*pow_3_of(f1 - f3)*pow_2_of(-f2 + f3)));
}

static double delta2_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d) {
  double f1 = p->f1;
  double f2 = p->f2;
  double f3 = p->f3;
  double v1 = p->v1;
  double v2 = p->v2;
  double v3 = p->v3;
  double d1 = p->d1;
  double d2 = p->d2;

  double f12 = d->f12;
  double f13 = d->f13;
  double f14 = d->f14;
  double f15 = d->f15;
  double f23 = d->f23;
  double f24 = d->f24;
  double f32 = d->f32;
  double f33 = d->f33;
  double f34 = d->f34;
  double f35 = d->f35;

  return -((d2*f15*f2 - d1*f13*f23 - 3*d2*f13*f23 + d1*f12*f24 + 2*d2*f12*f24 - d2*f15*f3
  + d2*f14*f2*f3 - d1*f12*f23*f3 + d2*f12*f23*f3 + d1*f1*f24*f3 - d2*f1*f24*f3 - d2*f14*f32
  + 3*d1*f13*f2*f32 + d2*f13*f2*f32 - d1*f1*f23*f32 + d2*f1*f23*f32 - 2*d1*f24*f32 - d2*f24*f32
  - 2*d1*f13*f33 + 2*d2*f13*f33 - d1*f12*f2*f33 - 3*d2*f12*f2*f33 + 3*d1*f23*f33 + d2*f23*f33
  + d1*f12*f34 - d1*f1*f2*f34 + d1*f1*f35 - d1*f2*f35 + 4*f12*f23*v1 - 3*f1*f24*v1 + 4*f1*f23*f3*v1
  - 3*f24*f3*v1 - 12*f12*f2*f32*v1 + 4*f23*f32*v1 + 8*f12*f33*v1 - f1*f34*v1 - f35*v1 - f15*v2
  - f14*f3*v2 + 8*f13*f32*v2 - 8*f12*f33*v2 + f1*f34*v2 + f35*v2 + f15*v3 - 4*f12*f23*v3 + 3*f1*f24*v3
  + f14*f3*v3 - 4*f1*f23*f3*v3 + 3*f24*f3*v3 - 8*f13*f32*v3 + 12*f12*f2*f32*v3 - 4*f23*f32*v3)
  / (pow_2_of(f1 - f2)*pow_3_of(f1 - f3)*pow_2_of(-f2 + f3)));
}

static double delta3_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d) {
  double f1 = p->f1;
  double f2 = p->f2;
  double f3 = p->f3;
  double v1 = p->v1;
  double v2 = p->v2;
  double v3 = p->v3;
  double d1 = p->d1;
  double d2 = p->d2;

  double f12 = d->f12;
  double f13 = d->f13;
  double f14 = d->f14;
  double f22 = d->f22;
  double f24 = d->f24;
  double f32 = d->f32;
  double f33 = d->f33;
  double f34 = d->f34;

  return -((-2*d2*f14*f2 + d1*f13*f22 + 3*d2*f13*f22 - d1*f1*f24 - d2*f1*f24 + 2*d2*f14*f3
  - 2*d1*f13*f2*f3 - 2*d2*f13*f2*f3 + d1*f12*f22*f3 - d2*f12*f22*f3 + d1*f24*f3 + d2*f24*f3
  + d1*f13*f32 - d2*f13*f32 - 2*d1*f12*f2*f32 + 2*d2*f12*f2*f32 + d1*f1*f22*f32 - d2*f1*f22*f32
  + d1*f12*f33 - d2*f12*f33 + 2*d1*f1*f2*f33 + 2*d2*f1*f2*f33 - 3*d1*f22*f33 - d2*f22*f33
  - 2*d1*f1*f34 + 2*d1*f2*f34 - 4*f12*f22*v1 + 2*f24*v1 + 8*f12*f2*f3*v1 - 4*f1*f22*f3*v1
  - 4*f12*f32*v1 + 8*f1*f2*f32*v1 - 4*f22*f32*v1 - 4*f1*f33*v1 + 2*f34*v1 + 2*f14*v2
  - 4*f13*f3*v2 + 4*f1*f33*v2 - 2*f34*v2 - 2*f14*v3 + 4*f12*f22*v3 - 2*f24*v3 + 4*f13*f3*v3
  - 8*f12*f2*f3*v3 + 4*f1*f22*f3*v3 + 4*f12*f32*v3 - 8*f1*f2*f32*v3 + 4*f22*f32*v3)
  / (pow_2_of(f1 - f2)*pow_3_of(f1 - f3)*pow_2_of(-f2 + f3)));
}

static double delta4_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d) {
  double f1 = p->f1;
  double f2 = p->f2;
  double f3 = p->f3;
  double v1 = p->v1;
  double v2 = p->v2;
  double v3 = p->v3;
  double d1 = p->d1;
  double d2 = p->d2;

  double f12 = d->f12;
  double f13 = d->f13;
  double f22 = d->f22;
  double f23 = d->f23;
  double f32 = d->f32;
  double f33 = d->f33;

  return -((d2*f13*f2 - d1*f12*f22 - 2*d2*f12*f22 + d1*f1*f23 + d2*f1*f23 - d2*f13*f3 + 2*d1*f12*f2*f3
  + d2*f12*f2*f3 - d1*f1*f22*f3 + d2*f1*f22*f3 - d1*f23*f3 - d2*f23*f3 - d1*f12*f32 + d2*f12*f32
  - d1*f1*f2*f32 - 2*d2*f1*f2*f32 + 2*d1*f22*f32 + d2*f22*f32 + d1*f1*f33 - d1*f2*f33 + 3*f1*f22*v1
  - 2*f23*v1 - 6*f1*f2*f3*v1 + 3*f22*f3*v1 + 3*f1*f32*v1 - f33*v1 - f13*v2 + 3*f12*f3*v2 - 3*f1*f32*v2
  + f33*v2 + f13*v3 - 3*f1*f22*v3 + 2*f23*v3 - 3*f12*f3*v3 + 6*f1*f2*f3*v3 - 3*f22*f3*v3)
  / (pow_2_of(f1 - f2)*pow_3_of(f1 - f3)*pow_2_of(-f2 + f3)));
}

/**
 * Calculates delta_i's
 * Method described in arXiv:1508.07253 section 'Region IIa - intermediate'
 */
static void ComputeDeltasFromCollocation(IMRPhenomDAmplitudeCoefficients* p) {
  // Three evenly spaced collocation points in the interval [f1,f3].
  double f1 = AMP_fJoin_INS;
  double f3 = p->fmaxCalc;
  double dfx = 0.5*(f3 - f1);
  double f2 = f1 + dfx;

  UsefulPowers powers_of_f1;
  int status = init_useful_powers(&powers_of_f1, f1);
  XLAL_CHECK_VOID ( status == XLAL_SUCCESS, XLAL_EFUNC, "Failed to initialize useful powers of f1.");

  AmpInsPrefactors prefactors;
  status = init_amp_ins_prefactors(&prefactors, p);
  XLAL_CHECK_VOID ( status == XLAL_SUCCESS, XLAL_EFUNC, "Failed to initialize amplitude prefactors for inspiral range.");


  // v1 is inspiral model evaluated at f1
  // d1 is derivative of inspiral model evaluated at f1
  double v1 = AmpInsAnsatz(f1, &powers_of_f1, &prefactors);
  double d1 = DAmpInsAnsatz(f1, &powers_of_f1, p);

  // v3 is merger-ringdown model evaluated at f3
  // d2 is derivative of merger-ringdown model evaluated at f3
  double v3 = AmpMRDAnsatz(f3, p);
  double d2 = DAmpMRDAnsatz(f3, p);

  // v2 is the value of the amplitude evaluated at f2
  // they come from the fit of the collocation points in the intermediate region
  double v2 = AmpIntColFitCoeff(p->eta, p->eta2, p->chi);

  p->f1 = f1;
  p->f2 = f2;
  p->f3 = f3;
  p->v1 = v1;
  p->v2 = v2;
  p->v3 = v3;
  p->d1 = d1;
  p->d2 = d2;

  // Now compute the delta_i's from the collocation coefficients
  // Precompute common quantities here and pass along to delta functions.
  DeltaUtility d;
  d.f12 = f1*f1;
  d.f13 = f1*d.f12;
  d.f14 = f1*d.f13;
  d.f15 = f1*d.f14;
  d.f22 = f2*f2;
  d.f23 = f2*d.f22;
  d.f24 = f2*d.f23;
  d.f32 = f3*f3;
  d.f33 = f3*d.f32;
  d.f34 = f3*d.f33;
  d.f35 = f3*d.f34;
  p->delta0 = delta0_fun(p, &d);
  p->delta1 = delta1_fun(p, &d);
  p->delta2 = delta2_fun(p, &d);
  p->delta3 = delta3_fun(p, &d);
  p->delta4 = delta4_fun(p, &d);
}


///////////////////////////// Amplitude: glueing function ////////////////////////////

/**
 * A struct containing all the parameters that need to be calculated
 * to compute the phenomenological amplitude
 */
static void ComputeIMRPhenomDAmplitudeCoefficients(IMRPhenomDAmplitudeCoefficients *p, double eta, double chi1, double chi2, double finspin) {

  p->eta = eta;
  p->etaInv = 1./eta;
  p->chi1 = chi1;
  p->chi2 = chi2;
  p->chi12 = chi1*chi1;
  p->chi22 = chi2*chi2;
  double eta2 = eta*eta;
  p->eta2 = eta2;
  p->eta3 = eta*eta2;
  double Seta = sqrt(1.0 - 4.0*eta);
  p->Seta = Seta;
  p->SetaPlus1 = 1.0 + Seta;

  p->q = 0.5 * (1.0 + Seta - 2.0*eta) * p->etaInv;
  p->chi = chiPN(Seta, eta, chi1, chi2);
  double xi = -1.0 + p->chi;

  p->fRD = fring(eta, chi1, chi2, finspin);
  p->fDM = fdamp(eta, chi1, chi2, finspin);

  // Compute gamma_i's, rho_i's first then delta_i's
  p->gamma1 = gamma1_fun(eta, eta2, xi);
  p->gamma2 = gamma2_fun(eta, eta2, xi);
  p->gamma3 = gamma3_fun(eta, eta2, xi);

  p->fmaxCalc = fmaxCalc(p);

  p->rho1 = rho1_fun(eta, eta2, xi);
  p->rho2 = rho2_fun(eta, eta2, xi);
  p->rho3 = rho3_fun(eta, eta2, xi);

  // compute delta_i's
  ComputeDeltasFromCollocation(p);

}

// Call ComputeIMRPhenomDAmplitudeCoefficients() first!
/**
 * This function computes the IMR amplitude given phenom coefficients.
 * Defined in VIII. Full IMR Waveforms arXiv:1508.07253
 */
static double IMRPhenDAmplitude(double f, IMRPhenomDAmplitudeCoefficients *p, UsefulPowers *powers_of_f, AmpInsPrefactors * prefactors) {
  // Defined in VIII. Full IMR Waveforms arXiv:1508.07253
  // The inspiral, intermediate and merger-ringdown amplitude parts

  // Transition frequencies
  p->fInsJoin = AMP_fJoin_INS;
  p->fMRDJoin = p->fmaxCalc;

  double AmpPreFac = prefactors->amp0 * powers_of_f->m_seven_sixths;

  // split the calculation to just 1 of 3 possible mutually exclusive ranges

  if (!StepFunc_boolean(f, p->fInsJoin))	// Inspiral range
  {
	  double AmpIns = AmpPreFac * AmpInsAnsatz(f, powers_of_f, prefactors);
	  return AmpIns;
  }

  if (StepFunc_boolean(f, p->fMRDJoin))	// MRD range
  {
	  double AmpMRD = AmpPreFac * AmpMRDAnsatz(f, p);
	  return AmpMRD;
  }

  //	Intermediate range
  double AmpInt = AmpPreFac * AmpIntAnsatz(f, p);
  return AmpInt;
}

/********************************* Phase functions *********************************/

////////////////////////////// Phase: Ringdown functions ///////////////////////////

// alpha_i i=1,2,3,4,5 are the phenomenological intermediate coefficients depending on eta and chiPN
// PhiRingdownAnsatz is the ringdown phasing in terms of the alpha_i coefficients

/**
 * alpha 1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double alpha1Fit(double eta, double eta2, double xi) {

  return 43.31514709695348 + 638.6332679188081*eta
    + (-32.85768747216059 + 2415.8938269370315*eta - 5766.875169379177*eta2
    + (-61.85459307173841 + 2953.967762459948*eta - 8986.29057591497*eta2)*xi
    + (-21.571435779762044 + 981.2158224673428*eta - 3239.5664895930286*eta2)*xi*xi)*xi;
}

/**
 * alpha 2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double alpha2Fit(double eta, double eta2, double xi) {

  return -0.07020209449091723 - 0.16269798450687084*eta
  + (-0.1872514685185499 + 1.138313650449945*eta - 2.8334196304430046*eta2
  + (-0.17137955686840617 + 1.7197549338119527*eta - 4.539717148261272*eta2)*xi
  + (-0.049983437357548705 + 0.6062072055948309*eta - 1.682769616644546*eta2)*xi*xi)*xi;
}

/**
 * alpha 3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double alpha3Fit(double eta, double eta2, double xi) {

  return 9.5988072383479 - 397.05438595557433*eta
  + (16.202126189517813 - 1574.8286986717037*eta + 3600.3410843831093*eta2
  + (27.092429659075467 - 1786.482357315139*eta + 5152.919378666511*eta2)*xi
  + (11.175710130033895 - 577.7999423177481*eta + 1808.730762932043*eta2)*xi*xi)*xi;
}

/**
 * alpha 4 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double alpha4Fit(double eta, double eta2, double xi) {

  return -0.02989487384493607 + 1.4022106448583738*eta
  + (-0.07356049468633846 + 0.8337006542278661*eta + 0.2240008282397391*eta2
  + (-0.055202870001177226 + 0.5667186343606578*eta + 0.7186931973380503*eta2)*xi
  + (-0.015507437354325743 + 0.15750322779277187*eta + 0.21076815715176228*eta2)*xi*xi)*xi;
}

/**
 * alpha 5 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double alpha5Fit(double eta, double eta2, double xi) {

  return 0.9974408278363099 - 0.007884449714907203*eta
  + (-0.059046901195591035 + 1.3958712396764088*eta - 4.516631601676276*eta2
  + (-0.05585343136869692 + 1.7516580039343603*eta - 5.990208965347804*eta2)*xi
  + (-0.017945336522161195 + 0.5965097794825992*eta - 2.0608879367971804*eta2)*xi*xi)*xi;
}

/**
 * Ansatz for the merger-ringdown phase Equation 14 arXiv:1508.07253
 * Rholm was added when IMRPhenomHM (high mode) was added.
 * Rholm = fRD22/fRDlm. For PhenomD (only (l,m)=(2,2)) this is just equal
 * to 1. and PhenomD is recovered.
 * Taulm = fDMlm/fDM22. Ratio of ringdown damping times.
 * Again, when Taulm = 1.0 then PhenomD is recovered.
 */
static double PhiMRDAnsatzInt(double f, IMRPhenomDPhaseCoefficients *p, double Rholm, double Taulm)
{
  double sqrootf = sqrt(f);
  double fpow1_5 = f * sqrootf;
  // check if this is any faster: 2 sqrts instead of one pow(x,0.75)
  double fpow0_75 = sqrt(fpow1_5); // pow(f,0.75);

  return -(p->alpha2/f)
		 + (4.0/3.0) * (p->alpha3 * fpow0_75)
		 + p->alpha1 * f
		 + p->alpha4 * Rholm * atan((f - p->alpha5 * p->fRD) / (Rholm * p->fDM * Taulm));
}

/**
 * First frequency derivative of PhiMRDAnsatzInt
 * Rholm was added when IMRPhenomHM (high mode) was added.
 * Rholm = fRD22/fRDlm. For PhenomD (only (l,m)=(2,2)) this is just equal
 * to 1. and PhenomD is recovered.
 * Taulm = fDMlm/fDM22. Ratio of ringdown damping times.
 * Again, when Taulm = 1.0 then PhenomD is recovered.
 */
static double DPhiMRD(double f, IMRPhenomDPhaseCoefficients *p, double Rholm, double Taulm) {
  return ( p->alpha1 + p->alpha2/pow_2_of(f) + p->alpha3/pow(f,0.25)+ p->alpha4/(p->fDM * Taulm * (1 + pow_2_of(f - p->alpha5 * p->fRD)/(pow_2_of(p->fDM * Taulm * Rholm)))) ) * p->etaInv;
}

///////////////////////////// Phase: Intermediate functions /////////////////////////////

// beta_i i=1,2,3 are the phenomenological intermediate coefficients depending on eta and chiPN
// PhiIntAnsatz is the intermediate phasing in terms of the beta_i coefficients


// \[Beta]1Fit = PhiIntFitCoeff\[Chi]PNFunc[\[Eta], \[Chi]PN][[1]]

/**
 * beta 1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double beta1Fit(double eta, double eta2, double xi) {

  return 97.89747327985583 - 42.659730877489224*eta
  + (153.48421037904913 - 1417.0620760768954*eta + 2752.8614143665027*eta2
  + (138.7406469558649 - 1433.6585075135881*eta + 2857.7418952430758*eta2)*xi
  + (41.025109467376126 - 423.680737974639*eta + 850.3594335657173*eta2)*xi*xi)*xi;
}

/**
 * beta 2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double beta2Fit(double eta, double eta2, double xi) {

  return -3.282701958759534 - 9.051384468245866*eta
  + (-12.415449742258042 + 55.4716447709787*eta - 106.05109938966335*eta2
  + (-11.953044553690658 + 76.80704618365418*eta - 155.33172948098394*eta2)*xi
  + (-3.4129261592393263 + 25.572377569952536*eta - 54.408036707740465*eta2)*xi*xi)*xi;
}

/**
 * beta 3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double beta3Fit(double eta, double eta2, double xi) {

  return -0.000025156429818799565 + 0.000019750256942201327*eta
  + (-0.000018370671469295915 + 0.000021886317041311973*eta + 0.00008250240316860033*eta2
  + (7.157371250566708e-6 - 0.000055780000112270685*eta + 0.00019142082884072178*eta2)*xi
  + (5.447166261464217e-6 - 0.00003220610095021982*eta + 0.00007974016714984341*eta2)*xi*xi)*xi;
}


/**
 * ansatz for the intermediate phase defined by Equation 16 arXiv:1508.07253
 */
static double PhiIntAnsatz(double Mf, IMRPhenomDPhaseCoefficients *p) {
  // 1./eta in paper omitted and put in when need in the functions:
  // ComputeIMRPhenDPhaseConnectionCoefficients
  // IMRPhenDPhase
  return  p->beta1*Mf - p->beta3/(3.*pow_3_of(Mf)) + p->beta2*log(Mf);
}

/**
 * First frequency derivative of PhiIntAnsatz
 * (this time with 1./eta explicitly factored in)
 */
static double DPhiIntAnsatz(double Mf, IMRPhenomDPhaseCoefficients *p) {
  return (p->beta1 + p->beta3/pow_4_of(Mf) + p->beta2/Mf) * p->etaInv;
}

/**
 * temporary instance of DPhiIntAnsatz used when computing
 * coefficients to make the phase C(1) continuous between regions.
 */
static double DPhiIntTemp(double ff, IMRPhenomDPhaseCoefficients *p) {
  double etaInv = p->etaInv;
  double beta1 = p->beta1;
  double beta2 = p->beta2;
  double beta3 = p->beta3;
  double C2Int = p->C2Int;

  return C2Int + (beta1 + beta3/pow_4_of(ff) + beta2/ff)*etaInv;
}

///////////////////////////// Phase: Inspiral functions /////////////////////////////

// sigma_i i=1,2,3,4 are the phenomenological inspiral coefficients depending on eta and chiPN
// PhiInsAnsatzInt is a souped up TF2 phasing which depends on the sigma_i coefficients

/**
 * sigma 1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double sigma1Fit(double eta, double eta2, double xi) {

  return 2096.551999295543 + 1463.7493168261553*eta
  + (1312.5493286098522 + 18307.330017082117*eta - 43534.1440746107*eta2
  + (-833.2889543511114 + 32047.31997183187*eta - 108609.45037520859*eta2)*xi
  + (452.25136398112204 + 8353.439546391714*eta - 44531.3250037322*eta2)*xi*xi)*xi;
}

/**
 * sigma 2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double sigma2Fit(double eta, double eta2, double xi) {

  return -10114.056472621156 - 44631.01109458185*eta
  + (-6541.308761668722 - 266959.23419307504*eta + 686328.3229317984*eta2
  + (3405.6372187679685 - 437507.7208209015*eta + 1.6318171307344697e6*eta2)*xi
  + (-7462.648563007646 - 114585.25177153319*eta + 674402.4689098676*eta2)*xi*xi)*xi;
}

/**
 * sigma 3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double sigma3Fit(double eta, double eta2, double xi) {

  return 22933.658273436497 + 230960.00814979506*eta
  + (14961.083974183695 + 1.1940181342318142e6*eta - 3.1042239693052764e6*eta2
  + (-3038.166617199259 + 1.8720322849093592e6*eta - 7.309145012085539e6*eta2)*xi
  + (42738.22871475411 + 467502.018616601*eta - 3.064853498512499e6*eta2)*xi*xi)*xi;
}

/**
 * sigma 4 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double sigma4Fit(double eta, double eta2, double xi) {

  return -14621.71522218357 - 377812.8579387104*eta
  + (-9608.682631509726 - 1.7108925257214056e6*eta + 4.332924601416521e6*eta2
  + (-22366.683262266528 - 2.5019716386377467e6*eta + 1.0274495902259542e7*eta2)*xi
  + (-85360.30079034246 - 570025.3441737515*eta + 4.396844346849777e6*eta2)*xi*xi)*xi;
}

/**
 * Ansatz for the inspiral phase.
 * We call the LAL TF2 coefficients here.
 * The exact values of the coefficients used are given
 * as comments in the top of this file
 * Defined by Equation 27 and 28 arXiv:1508.07253
 */
static double PhiInsAnsatzInt(double Mf, UsefulPowers *powers_of_Mf, PhiInsPrefactors *prefactors, IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn)
{
  XLAL_CHECK(0 != pn, XLAL_EFAULT, "pn is NULL");

  // Assemble PN phasing series
  const double v = powers_of_Mf->third * powers_of_pi.third;
  const double logv = log(v);

  double phasing = prefactors->initial_phasing;

  phasing += prefactors->two_thirds * powers_of_Mf->two_thirds;
  phasing += prefactors->third * powers_of_Mf->third;
  phasing += prefactors->third_with_logv * logv * powers_of_Mf->third;
  phasing += prefactors->logv * logv;
  phasing += prefactors->minus_third * powers_of_Mf->m_third;
  phasing += prefactors->minus_two_thirds * powers_of_Mf->m_two_thirds;
  phasing += prefactors->minus_one * powers_of_Mf->inv;
  phasing += prefactors->minus_four_thirds / powers_of_Mf->four_thirds;
  phasing += prefactors->minus_five_thirds * powers_of_Mf->m_five_thirds; // * v^0

  // Now add higher order terms that were calibrated for PhenomD
  phasing += ( prefactors->one * Mf + prefactors->four_thirds * powers_of_Mf->four_thirds
			   + prefactors->five_thirds * powers_of_Mf->five_thirds
			   + prefactors->two * powers_of_Mf->two
        ) * p->etaInv;

  return phasing;
}

static int init_phi_ins_prefactors(PhiInsPrefactors * prefactors, IMRPhenomDPhaseCoefficients* p, PNPhasingSeries *pn)
{
	XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");
	XLAL_CHECK(0 != prefactors, XLAL_EFAULT, "prefactors is NULL");

	double sigma1 = p->sigma1;
	double sigma2 = p->sigma2;
	double sigma3 = p->sigma3;
	double sigma4 = p->sigma4;

  // PN phasing series
	prefactors->initial_phasing = pn->v[5] - LAL_PI_4;
	prefactors->two_thirds = pn->v[7] * powers_of_pi.two_thirds;
	prefactors->third = pn->v[6] * powers_of_pi.third;
	prefactors->third_with_logv = pn->vlogv[6] * powers_of_pi.third;
	prefactors->logv = pn->vlogv[5];
	prefactors->minus_third = pn->v[4] * powers_of_pi.m_third;
	prefactors->minus_two_thirds = pn->v[3] * powers_of_pi.m_two_thirds;
	prefactors->minus_one = pn->v[2] * powers_of_pi.inv;
  prefactors->minus_four_thirds = pn->v[1] / powers_of_pi.four_thirds;
  prefactors->minus_five_thirds = pn->v[0] * powers_of_pi.m_five_thirds; // * v^0

  // higher order terms that were calibrated for PhenomD
	prefactors->one = sigma1;
	prefactors->four_thirds = sigma2 * 0.75;
	prefactors->five_thirds = sigma3 * 0.6;
	prefactors->two = sigma4 * 0.5;

	return XLAL_SUCCESS;
}

/**
 * First frequency derivative of PhiInsAnsatzInt
 */
static double DPhiInsAnsatzInt(double Mf, IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn) {
  double sigma1 = p->sigma1;
  double sigma2 = p->sigma2;
  double sigma3 = p->sigma3;
  double sigma4 = p->sigma4;
  double Pi = LAL_PI;

  // Assemble PN phasing series
  const double v = cbrt(Pi*Mf);
  const double logv = log(v);
  const double v2 = v * v;
  const double v3 = v * v2;
  const double v4 = v * v3;
  const double v5 = v * v4;
  const double v6 = v * v5;
  const double v7 = v * v6;
  const double v8 = v * v7;

  // Apply the correct prefactors to LAL phase coefficients to get the
  // phase derivative dphi / dMf = dphi/dv * dv/dMf
  double Dphasing = 0.0;
  Dphasing += 2.0 * pn->v[7] * v7;
  Dphasing += (pn->v[6] + pn->vlogv[6] * (1.0 + logv)) * v6;
  Dphasing += pn->vlogv[5] * v5;
  Dphasing += -1.0 * pn->v[4] * v4;
  Dphasing += -2.0 * pn->v[3] * v3;
  Dphasing += -3.0 * pn->v[2] * v2;
  Dphasing += -4.0 * pn->v[1] * v;
  Dphasing += -5.0 * pn->v[0];
  Dphasing /= v8 * 3.0;
  Dphasing *= Pi;

  // Now add higher order terms that were calibrated for PhenomD
  Dphasing += (
          sigma1
        + sigma2 * v * powers_of_pi.m_third
        + sigma3 * v2 * powers_of_pi.m_two_thirds
        + (sigma4*powers_of_pi.inv) * v3
        ) * p->etaInv;

  return Dphasing;
}

///////////////////////////////// Phase: glueing function ////////////////////////////////

/**
 * A struct containing all the parameters that need to be calculated
 * to compute the phenomenological phase
 */
static void ComputeIMRPhenomDPhaseCoefficients(IMRPhenomDPhaseCoefficients *p, double eta, double chi1, double chi2, double finspin, LALDict *extraParams) {

  // Convention m1 >= m2
  p->eta = eta;
  p->etaInv = 1./eta;
  p->chi1 = chi1;
  p->chi2 = chi2;
  double eta2 = eta*eta;
  p->eta2 = eta2;
  p->Seta = sqrt(1.0 - 4.0*eta);

  p->q = 0.5*(1.0 + p->Seta - 2.0*eta)*p->etaInv;
  p->chi = chiPN(p->Seta, eta, chi1, chi2);
  double xi = -1.0 + p->chi;

  p->sigma1 = sigma1Fit(eta, eta2, xi);
  p->sigma2 = sigma2Fit(eta, eta2, xi);
  p->sigma3 = sigma3Fit(eta, eta2, xi);
  p->sigma4 = sigma4Fit(eta, eta2, xi);

  p->beta1 = beta1Fit(eta, eta2, xi);
  p->beta2 = beta2Fit(eta, eta2, xi);
  p->beta3 = beta3Fit(eta, eta2, xi);

  p->alpha1 = alpha1Fit(eta, eta2, xi);
  p->alpha2 = alpha2Fit(eta, eta2, xi);
  p->alpha3 = alpha3Fit(eta, eta2, xi);
  p->alpha4 = alpha4Fit(eta, eta2, xi);
  p->alpha5 = alpha5Fit(eta, eta2, xi);

  p->fRD = fring(eta, chi1, chi2, finspin);
  p->fDM = fdamp(eta, chi1, chi2, finspin);

  p->sigma1*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDSigma1(extraParams));
  p->sigma2*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDSigma2(extraParams));
  p->sigma3*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDSigma3(extraParams));
  p->sigma4*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDSigma4(extraParams));
  p->beta1*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDBeta1(extraParams));
  p->beta2*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDBeta2(extraParams));
  p->beta3*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDBeta3(extraParams));
  p->alpha1*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDAlpha1(extraParams));
  p->alpha2*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDAlpha2(extraParams));
  p->alpha3*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDAlpha3(extraParams));
  p->alpha4*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDAlpha4(extraParams));
  p->alpha5*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDAlpha5(extraParams));

}

/**
 * This function aligns the three phase parts (inspiral, intermediate and merger-rindown)
 * such that they are c^1 continuous at the transition frequencies
 * Defined in VIII. Full IMR Waveforms arXiv:1508.07253
 * Rholm was added when IMRPhenomHM (high mode) was added.
 * Rholm = fRD22/fRDlm. For PhenomD (only (l,m)=(2,2)) this is just equal
 * to 1. and PhenomD is recovered.
 * Taulm = fDMlm/fDM22. Ratio of ringdown damping times.
 * Again, when Taulm = 1.0 then PhenomD is recovered.
 */
static void ComputeIMRPhenDPhaseConnectionCoefficients(IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn, PhiInsPrefactors *prefactors, double Rholm, double Taulm)
{
  double etaInv = p->etaInv;

  // Transition frequencies
  // Defined in VIII. Full IMR Waveforms arXiv:1508.07253
  p->fInsJoin=PHI_fJoin_INS;
  p->fMRDJoin=0.5*p->fRD;

  // Compute C1Int and C2Int coeffs
  // Equations to solve for to get C(1) continuous join
  // PhiIns (f)  =   PhiInt (f) + C1Int + C2Int f
  // Joining at fInsJoin
  // PhiIns (fInsJoin)  =   PhiInt (fInsJoin) + C1Int + C2Int fInsJoin
  // PhiIns'(fInsJoin)  =   PhiInt'(fInsJoin) + C2Int
  double DPhiIns = DPhiInsAnsatzInt(PHI_fJoin_INS, p, pn);
  double DPhiInt = DPhiIntAnsatz(PHI_fJoin_INS, p);
  p->C2Int = DPhiIns - DPhiInt;

  UsefulPowers powers_of_fInsJoin;
  init_useful_powers(&powers_of_fInsJoin, PHI_fJoin_INS);
  p->C1Int = PhiInsAnsatzInt(PHI_fJoin_INS, &powers_of_fInsJoin, prefactors, p, pn)
    - etaInv * PhiIntAnsatz(PHI_fJoin_INS, p) - p->C2Int * PHI_fJoin_INS;

  // Compute C1MRD and C2MRD coeffs
  // Equations to solve for to get C(1) continuous join
  // PhiInsInt (f)  =   PhiMRD (f) + C1MRD + C2MRD f
  // Joining at fMRDJoin
  // Where \[Phi]InsInt(f) is the \[Phi]Ins+\[Phi]Int joined function
  // PhiInsInt (fMRDJoin)  =   PhiMRD (fMRDJoin) + C1MRD + C2MRD fMRDJoin
  // PhiInsInt'(fMRDJoin)  =   PhiMRD'(fMRDJoin) + C2MRD
  // temporary Intermediate Phase function to Join up the Merger-Ringdown
  double PhiIntTempVal = etaInv * PhiIntAnsatz(p->fMRDJoin, p) + p->C1Int + p->C2Int*p->fMRDJoin;
  double DPhiIntTempVal = DPhiIntTemp(p->fMRDJoin, p);
  double DPhiMRDVal = DPhiMRD(p->fMRDJoin, p, Rholm, Taulm);
  p->C2MRD = DPhiIntTempVal - DPhiMRDVal;
  p->C1MRD = PhiIntTempVal - etaInv * PhiMRDAnsatzInt(p->fMRDJoin, p, Rholm, Taulm) - p->C2MRD*p->fMRDJoin;
}

/**
 * This function computes the IMR phase given phenom coefficients.
 * Defined in VIII. Full IMR Waveforms arXiv:1508.07253
 * Rholm was added when IMRPhenomHM (high mode) was added.
 * Rholm = fRD22/fRDlm. For PhenomD (only (l,m)=(2,2)) this is just equal
 * to 1. and PhenomD is recovered.
 * Taulm = fDMlm/fDM22. Ratio of ringdown damping times.
 * Again, when Taulm = 1.0 then PhenomD is recovered.
 */
static double IMRPhenDPhase(double f, IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn, UsefulPowers *powers_of_f, PhiInsPrefactors *prefactors, double Rholm, double Taulm)
{
  // Defined in VIII. Full IMR Waveforms arXiv:1508.07253
  // The inspiral, intermendiate and merger-ringdown phase parts

  // split the calculation to just 1 of 3 possible mutually exclusive ranges
  if (!StepFunc_boolean(f, p->fInsJoin))	// Inspiral range
  {
	  double PhiIns = PhiInsAnsatzInt(f, powers_of_f, prefactors, p, pn);
	  return PhiIns;
  }

  if (StepFunc_boolean(f, p->fMRDJoin))	// MRD range
  {
	  double PhiMRD = p->etaInv * PhiMRDAnsatzInt(f, p, Rholm, Taulm) + p->C1MRD + p->C2MRD * f;
	  return PhiMRD;
  }

  //	Intermediate range
  double PhiInt = p->etaInv * PhiIntAnsatz(f, p) + p->C1Int + p->C2Int * f;
  return PhiInt;
}

/**
 * Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
 * (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
 * was not available when PhenomD was tuned.
 */
static double Subtract3PNSS(double m1, double m2, double M, double eta, double chi1, double chi2){
  REAL8 m1M = m1 / M;
  REAL8 m2M = m2 / M;
  REAL8 pn_ss3 =  (326.75L/1.12L + 557.5L/1.8L*eta)*eta*chi1*chi2;
  pn_ss3 += ((4703.5L/8.4L+2935.L/6.L*m1M-120.L*m1M*m1M) + (-4108.25L/6.72L-108.5L/1.2L*m1M+125.5L/3.6L*m1M*m1M)) *m1M*m1M * chi1*chi1;
  pn_ss3 += ((4703.5L/8.4L+2935.L/6.L*m2M-120.L*m2M*m2M) + (-4108.25L/6.72L-108.5L/1.2L*m2M+125.5L/3.6L*m2M*m2M)) *m2M*m2M * chi2*chi2;
  return pn_ss3;
}
