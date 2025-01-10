#include <stdio.h>
#include <math.h>
#include "bvia.h"
#include "dprcvia.c"
#include "cblock.h"

// biexponential prc

/* global variable declarations */

/* double F[N];
double Y[N];
long int *np;
double *xp; */
 
double gsyn1=0.0,gsyn2=G_SYN2,current[C],vset1,initial[N],period[2];
int flag1=0,flag2=0,flag3=0,flag4=0,pflag=0;
double f1,f2,f3,ts,forcingt;

void main()
{
double state[N];
double fstate[N];
double rtol[N],atol[N];

int itol,ijac,mljac,mujac;
int imas,mlmas,mumas;
int iout,lwork,liwork,lrcont;
int iwork[LIWORK];
double work[LWORK];
extern int out_(),altout_();
extern int flag1,flag2,flag3,flag4,pflag;
extern double f1,f2,f3;
int n,idid,counter=0;

double tau_2,arg,f,oldpf=0.0,oldpn=0.0,newpn,m,b,intersection,delay;
double x,xend,h;
FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*pi,*ip;
int i,j,k;
fp1 = fopen("prediction.data","w");
//fp2 = fopen("f21.data","w");
//fp3 = fopen("f31.data","w");
pi = fopen("period.data","w");
ip = fopen("initial.data","w");
n=N;
/* housekeeping for numerical integration */
for(i=0;i<7;i++) {iwork[i]=0;
                  work[i]=0.0;}
for(i=1;i<N;i++) {rtol[i]=1.0e-8; 
                  atol[i]=1.0e-14;}
rtol[0]=1.0e-7;
atol[0]=1.0e-13;

iwork[2] = 10000000;
itol=1;
ijac=0;
mljac=n;
mujac=0;
imas=0;
mlmas=n;
mumas=0;
lwork=LWORK;
liwork=LIWORK;
lrcont=LRCONT;
iout=1;
// loop thorugh delay values
for(k=0;k<60;k++)
{
	delay= k*0.1;
	printf("%f\n",delay);
	fflush(stdout);
//Mean field loop through forcing frequencies
if (!PRINT) {for(counter=0;counter<INCREMENT+1;counter++)
{
  forcingt =  PF_MIN + counter*(PF_MAX-PF_MIN)/(1.0*INCREMENT) ;
h=1e-5;
x=START_TIME;
/* initialize state variables from "state.data"*/
/* gsyn2 has tonic component from stimulus train*/
arg=(TAU_TWO*TAU_FALL)*log(TAU_FALL/TAU_TWO)/(TAU_FALL-TAU_TWO);
f = 1.0/(exp(-arg/TAU_FALL)-exp(-arg/TAU_TWO));
gsyn2=(1.0/(1-exp(-forcingt/TAU_FALL))) -(1.0/(1.0-exp(-forcingt/TAU_TWO)));
gsyn2=f*gsyn2*G_SYN;
scan_(state);
deriv_(&n,&x,state,fstate);
xend=ENDTIME;
tau_2=TAU_TWO;
// tau_2=(TAU_RISE*TAU_FALL)/(TAU_FALL - TAU_RISE); /* from Abbott and Dayan page 182*/
// printf("%f\n",tau_2);
// run the code long enough that the limit cycle is at steady state to determine period
// and state variables at spike threshold stored in "initial.data"
radau5_(&n,deriv_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        out_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);
/* loop through forcing frequencies at constant delay */
pflag=1;
   x = 0.0;
  gsyn1 = 0.0;
  h=1e-5;
  xend = delay;
  ts=xend;
  for(n=0;n<N;n++) state[n]=initial[n];
  flag1=0;
  if(xend>0) radau5_(&n,deriv_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        altout_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);
//  h=1e-4;
  flag2=0;
  flag3=0;
  flag4=0;
ts=x;
// apply the perturbation (note switch to aderiv() subroutine) see "bprcvia.c"
// then integrate until 3 spikes have been emitted
// in order to capture 1st, 2nd and 3rd order PRCs in "f1.data", "f2.data" and "f3.data"
// if too much time elapses between spikes, the loop is exited without printing see "pviaprcout.c"
// need to make sure flag1 is not set, ts would be greater the entrained period 
// if period[0] then the forcing silenced the neuron
if((period[0]>0.0)&&(!flag1)) {
 do{
  xend = xend + 0.1;
  //printf("%f %f\n",x-ts,gsyn1);
  radau5_(&n,aderiv_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        altout_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);
        } while (!flag3);
 if(x<ENDTIME)
  {
// fprintf(fp1,"%f %f \n",forcingt, period[0]*(1.0+f1));
 gsyn2= gsyn2/(1.0*G_SYN);
 newpn= period[0]*(1.0+f1 +f2);
 m=(newpn-oldpn)/(forcingt-oldpf);
 b=oldpn-m*oldpf;
 intersection= b/(1-m);
 //if((oldpf<=oldpn)&& (forcingt>=period[0]*(1.0+f1))) fprintf(fp1,"%f %f\n",delay,1000.0/intersection);

 if((oldpf<=oldpn)&& (forcingt>=period[0]*(1.0+f1+f2))) fprintf(fp1,"%f %f\n",delay,1000.0/intersection);
 oldpf=forcingt;
 oldpn=period[0]*(1.0+f1+f2);
 //fprintf(fp1,"%f %f %f\n",forcingt, period[0]*(1.0+f1),gsyn2);
 fflush(fp1);
 // fprintf(fp2,"%f %f\n",(1.0*counter)/(1.0*INCREMENT),f2);
 // fprintf(fp3,"%f %f\n",(1.0*counter)/(1.0*INCREMENT),f3);
   }
 }
pflag=0;
 }
}
// end of forcing frequency loop for each delay
}

  fprintf(pi,"%f \n",period[0]);
  for(n=0;n<N;n++) fprintf(ip,"%f\n",initial[n]);
fflush(stdout);
fclose(fp1);
//fclose(fp2);
//fclose(fp3);
fclose(pi);
fclose(ip);
}

