/* out.c */
#include <math.h>
#include <stdio.h>
#include "avia.h"

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
static double tcross,yprev2,yold2,yolder2,last_t=0.0,last_t2;
extern double gsyn1,current[C],initial[N],period[2];
extern int gflag;
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
 //  printf("1 last_t %f %f\n",last_t,tcross);
   if(tcross>= last_t+DELAY) last_t=tcross;
 //  printf("1 last_t %f %f\n",last_t,tcross);
    }
   //if(t_old<last_t+DELAY&&*x>last_t+DELAY&&last_t>0.0) {gflag=1;
   if(t_old<last_t+DELAY&&*x>last_t+DELAY) {gflag=1;
 //  printf("2 last_t %f %f\n",last_t,tcross);
	   last_t=tcross;
//   printf("3 last_t %f %f\n",last_t,tcross);
    }
yolder=yold;
yold=state[0];
t_old=*x;
for(i=0;i<N;i++) oldstate[i]=state[i];
if (FREERUN) printf("%f %f\n",*x,state[0]);
//if (PRINT) printf("%f %f\n",*x,state[0]);
// (PRINT) printf("%f %f %f %f %f %f\n",*x,state[0],state[1],state[2],state[3],state[4]);
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
extern int gflag,flag1,flag2,flag3,flag4,clearflag;
static double oldcurrent[C];
static double oldstate[N];
extern double gsyn1,phi();
int i;
//if (PRINT) printf("%f %f\n",*x,gsyn1/(1.0*G_SYN));
if (PRINT) printf("%f %f\n",*x,state[0]);
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
if(clearflag){last_t=0;
	      tcross=0;
	      clearflag=0;}
if( (yold < THRESHOLD) && (state[0] >= THRESHOLD) ) 
   { tcross = t_old + (*x - t_old)*(THRESHOLD - yold)/(state[0] - yold);
   if(flag3) flag4=1;
   if(flag2) flag3=1;
   if(flag1) flag2=1;
   flag1 = 1;
   //printf("flag %d %d %d\n",flag1,flag2,flag3);
   if((flag1==1)&& (flag2==0)) f1 = (tcross - period[0])/period[0];
   if((flag2==1)&& (flag3==0)) f2 = (tcross - last_t - period[0])/period[0];
   if((flag3==1)&& (!flag4)) f3 = (tcross -last_t - period[0])/period[0];
   if(tcross>= last_t+DELAY) last_t=tcross;
    //printf("1 last_t %f %f\n",last_t,tcross);
   }
   if(t_old<last_t+DELAY&&*x>=last_t+DELAY) {gflag=1;
    //printf("2 last_t %f %f\n",last_t,tcross);
	   last_t=tcross;
    //printf("3 last_t %f %f\n",last_t,tcross);
	                                    }
//                                    printf("gflag %d\n",gflag);
  //printf("%f %f %f\n", t_old,*x,last_t+DELAY);
yolder=yold;
yold=state[0];
t_old=*x;
if(*x>ENDTIME) flag3=1;
for(i=0;i<N;i++) oldstate[i]=state[i];
for(i=0;i<C;i++) oldcurrent[i]=current[i];
}

