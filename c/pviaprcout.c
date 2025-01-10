/* out.c */
#include <math.h>
#include <stdio.h>
#include "bvia.h"

#define T_RES 1
#define V_LOW_RES 2.0
#define V_HIGH_RES 0.05

int out_(nr,xold,x,state,n,irtrn)
        int *nr;
        double *xold;
        double *x;
        double *state;
        int *n;
        int *irtrn;
{
static double t_old,yold,yolder,yprev,tprev;
static double tcross,yprev2,yold2,yolder2,last_t,last_t2;
extern double current[C],initial[N],period[2];
static double oldcurrent[C];
static double oldstate[N];
extern double phi();
int i;
if(*nr==1) 
{
yolder=yold;
yold=state[0];
yprev=yold;
t_old=*x;}
if( (yold < THRESHOLD) && (state[0] > THRESHOLD) ) 
   { tcross = t_old + (*x - t_old)*(THRESHOLD - yold)/(state[0] - yold);
   period[0] = tcross - last_t;
   initial[0] = THRESHOLD;
   for(i=1; i<N;i++) initial[i] = oldstate[i] + (state[i] - oldstate[i])*(THRESHOLD - yold)/(state[0] - yold);
    last_t = tcross;}
yolder=yold;
yold=state[0];
t_old=*x;
for(i=0;i<N;i++) oldstate[i]=state[i];
 if (PRINT) printf("%f %f\n",*x,state[0]);
//f (PRINT) printf("%f %f %f %f %f %f\n",*x,state[0],state[1],state[2],state[3],state[4]);
}
int altout_(nr,xold,x,state,n,irtrn)
        int *nr;
        double *xold;
        double *x;
        double *state;
        int *n;
        int *irtrn;
{
static double t_old,yold,yolder,yprev,tprev;
static double tcross,yprev2,yold2,yolder2,last_t,last_t2;
extern double current[C],initial[N],period[2];
extern double f1,f2,f3;
extern int flag1,flag2,flag3,flag4;
static double oldcurrent[C];
static double oldstate[N];
extern double gsyn1,gsyn2,phi();
double scale;
int i;
scale=1.0/(1.0*G_SYN);
if (DEBUG) printf("%f %f\n",*x,state[0]);
//if (DEBUG) printf("%f %f\n",*x,scale*(gsyn1));
//printf("%f ",*x);
//for(i=0;i<N;i++) printf("%f ",state[i]);
//printf("\n");
//fflush(stdout);
if(*nr==1) 
{
yolder=yold;
yold=state[0];
yprev=yold;
t_old=*x;}
if( (yold < THRESHOLD) && (state[0] > THRESHOLD) ) 
   { tcross = t_old + (*x - t_old)*(THRESHOLD - yold)/(state[0] - yold);
   if(flag3) flag4=1;
   if(flag2) flag3=1;
   if(flag1) flag2=1;
   flag1 = 1;
   if((flag1==1)&& (flag2==0)) 
   {f1 = (tcross - period[0])/period[0];
//   if(DEBUG) printf("flag1 %f %f %f\n", tcross, *x,state[0]);
   } 
   if((flag2==1)&& (flag3==0)) 
   {f2 = (tcross - last_t - period[0])/period[0];
//   if(DEBUG) printf("flag2 %f %f %f\n", tcross, *x,state[0]);
   } 
   if((flag3==1)&& (!flag4))
   {f3 = (tcross -last_t - period[0])/period[0];
//   if(DEBUG) printf("flag3 %f %f %f\n", tcross, *x,state[0]);
   } 
    last_t=tcross;}
yolder=yold;
yold=state[0];
t_old=*x;
if(*x>ENDTIME) flag3=1;
for(i=0;i<N;i++) oldstate[i]=state[i];
for(i=0;i<C;i++) oldcurrent[i]=current[i];
}
