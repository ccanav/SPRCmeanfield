#define FILELN 38
#define TAU_TWO 0.3
#define TAU_FALL 2.0
#define THRESHOLD -30.0
#define INCREMENT 100
#define LIGHT 1
#define N 5 // number of state variables
#define C 6 // number of currents
#define V1 0  /* mV */
#define H1 1  
#define N1 2 
#define A1 3 
#define M1 4 
#define I_NA1 0  
#define I_K1 1 
#define I_L1 2
#define I_S1 3 
#define I_K3 4
#define I_S2 5 
#define CM 1.00 /*uF/cm2*/
#define I_APP 0.00 /*uA/cm2*/
#define E_NA  50.0
#define E_K  -90.0
#define E_L  -71.975 //-70.4
#define E_SYN  -75.0 /*reversal potential of GABA_Asynapse*/
#define G_NA 218.8 //239.13  mS/cm2 all conductances were divided by 0.0000768 cm2
#define G_K3   8.225 // 10.22469
#define G_K1   0.769 // 1.1058
#define G_L   0.1915 // 0.2224
#define G_SYN   35*0.0214844 // 36 single synapse value is 0.02148 mS/cm2*/ 
#define G_SYN2  0.09115 // equal to 7 nS ChR*/ 
#define THETA_M -52.995 // -56.42
#define SIGMA1_M 4.0
#define SIGMA2_M -13.0
#define K1_M 0.25
#define K2_M 0.1
#define THETA_H -55.711 //-56.59
#define SIGMA1_H -20.0
#define SIGMA2_H 3.5
#define K1_H 0.012
#define K2_H 0.2
#define THETA_N 5.9 // -4.23
#define SIGMA1_N 12.0
#define SIGMA2_N -8.5
#define K1_N 1.0 
#define K2_N 0.001
#define THETA_A 51.355 // 51.4
#define SIGMA1_A 12.0
#define SIGMA2_A -80.0
#define K1_A 1.0 
#define K2_A 0.02
#define ALPHA 6.25 /* ms */
#define TAUSYN 1.0 /*  1.0 3.0 2.0   */
#define SS1 0.00
#define SS2 0.00
#define I_STIM 0.0
#define FORCEPERIOD 125.0


/* INTEGRATION PARAMETERS */

/* these only apply to C code */
#define LWORK 327 /* N(4*N + 8)+7*/
#define LIWORK 31 /*3*N+7*/
#define LRCONT 36 /*4*N+4*/
#define DURATION 1.0  
#define TIMEON 0
#define TIMEOFF 0
#define PERIOD  250
#define ENDTIME 100
#define START_TIME  0

#define PRINT  0
#define PI 3.1614
#define DEBUG 0
#define PULSEON 300000
#define PULSEOFF 300000
#define T 13.068  /*3.99 for 0.8 ms delay  forcing period for tonic/phasic PRC */
#define PF_MIN 2.0 // upper frequency range 250 Hz, 4 ms
#define PF_MAX 20.0 // upper frequency range 50 Hz, 20 ms
#define DELAY 3.0 //for dviprc loop through forcing frequencies at constant DELAY


