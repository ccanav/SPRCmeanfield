#include <math.h>
#include <stdio.h>


double efun(z) 
double z;
{
	if (fabs(z) < 1e-4) {
	return( 1 - z/2);
	}else{
		return( z/(exp(z) - 1));
	}
}
 
 
double  gaussian(v,a,b,c,d) 
double  v,a,b,c,d ;
{
double arg;
arg = pow(((v-c)/b),2.0);
if ((arg>50.0) || (arg<-50.0))
{if (arg>50.0) return(d);
else  return(d+a);}
else return  (d + a*exp(-arg));
}


double boltz(v,half,slope)
double v,half,slope;
{
double arg;
arg = -(v-half)/slope;
if ((arg>50.0) || (arg<-50.0))
{if (arg>50.0) return(0.0);
else  return(1.0);}
else return(1.0/(1.0 + exp(arg)));}


double alpha(v,k,theta,sigma)
double v,k,theta,sigma;
{
double arg;
arg = (theta-v)/sigma;
if(v==theta) return(k*sigma);// L'Hopital's rule
else {
if ((arg>50.0) || (arg<-50.0))
{if (arg>50.0) return(0.0);
else  return(k*(v-theta));} 
else return(k*(theta-v)/(exp(arg) - 1.0));
}
}

double beta(v,k,sigma)
double v,k,sigma;
{ return(k*exp(v/sigma));}
 
int deriv_(np,xp,Y,F)
double *F,*Y;
int *np;
double *xp;
{
extern double  ts,gsyn1,gsyn2,current[C]; 
extern int pflag;
int el;
double time_;
double tau_2,f,arg,theta,k1,sigma1,k2,sigma2;
time_ = *xp;

el = *np;
tau_2=TAU_TWO;
arg=(TAU_TWO*TAU_FALL)*log(TAU_FALL/TAU_TWO)/(TAU_FALL-TAU_TWO);
f = 1.0/(exp(-arg/TAU_FALL)-exp(-arg/TAU_TWO));
if(pflag&&*xp>ts-T) gsyn1=f*G_SYN*((exp((ts-T-*xp)/TAU_FALL)-1.0)/(1.0 - exp(-T/TAU_FALL))-(exp((ts-T-*xp)/tau_2)-1.0)/(1.0-exp(-T/tau_2)));
else gsyn1=0.0;
current[I_NA1] = G_NA*pow(Y[M1],3.0)*Y[H1]*(Y[V1] - E_NA);
current[I_K1] = G_K1*pow(Y[A1],4.0)*(Y[V1] - E_K);
current[I_K3] = G_K3*pow(Y[N1],4.0)*(Y[V1] - E_K);
current[I_L1] =  G_L*(Y[V1] - E_L);
current[I_S1] =  gsyn1*(Y[V1] - E_SYN);
current[I_S2] =  gsyn2*(Y[V1]-E_SYN) + G_SYN2*Y[V1] ;
F[V1] = (I_STIM-current[I_S2] - current[I_NA1] - current[I_K1] - current[I_K3] - current[I_L1] - current[I_S1])/CM;
theta = THETA_H; k1=K1_H; k2=K2_H; sigma1=SIGMA1_H; sigma2=SIGMA2_H;
F[H1] = beta(Y[V1],k1,sigma1)*(1.0 - Y[H1]) - alpha(Y[V1],k2,theta,sigma2)*Y[H1];
theta = THETA_N; k1=K1_N; k2=K2_N; sigma1=SIGMA1_N; sigma2=SIGMA2_N;
F[N1] = alpha(Y[V1],k1,theta,sigma1)*(1.0 - Y[N1]) - beta(Y[V1],k2,sigma2)*Y[N1];
theta = THETA_A; k1=K1_A; k2=K2_A; sigma1=SIGMA1_A; sigma2=SIGMA2_A;
F[A1] = alpha(Y[V1],k1,theta,sigma1)*(1.0 - Y[A1]) - beta(Y[V1],k2,sigma2)*Y[A1];
theta = THETA_M; k1=K1_M; k2=K2_M; sigma1=SIGMA1_M; sigma2=SIGMA2_M;
F[M1] = alpha(Y[V1],k1,theta,sigma1)*(1.0 - Y[M1]) - beta(Y[V1],k2,sigma2)*Y[M1];

return 0;
   }


int aderiv_(np,xp,Y,F)
double *F,*Y;
int *np;
double *xp;
{
extern double ts, gsyn1,gsyn2,current[C]; 
int el;
double time_;
double theta,k1,sigma1,k2,sigma2;
double f,arg,tau_2;
time_ = *xp;
//tau_2=(TAU_RISE*TAU_FALL)/(TAU_FALL - TAU_RISE);
tau_2=TAU_TWO;
arg=(TAU_TWO*TAU_FALL)*log(TAU_FALL/TAU_TWO)/(TAU_FALL-TAU_TWO);
f = 1.0/(exp(-arg/TAU_FALL)-exp(-arg/TAU_TWO));
//gsyn1=f*G_SYN*(exp((ts-*xp)/TAU_FALL)-exp((ts-*xp)/tau_2));
if(*xp < T) gsyn1=f*G_SYN*((exp((ts-*xp)/TAU_FALL)-1.0)/(1.0 - exp(-T/TAU_FALL))-(exp((ts-*xp)/tau_2)-1.0)/(1.0-exp(-T/tau_2)));
else gsyn1=0.0;

el = *np;
current[I_NA1] = G_NA*pow(Y[M1],3.0)*Y[H1]*(Y[V1] - E_NA);
current[I_K1] = G_K1*pow(Y[A1],4.0)*(Y[V1] - E_K);
current[I_K3] = G_K3*pow(Y[N1],4.0)*(Y[V1] - E_K);
current[I_L1] =  G_L*(Y[V1] - E_L);
current[I_S1] =  gsyn1*(Y[V1] - E_SYN);
current[I_S2] =  gsyn2*(Y[V1]-E_SYN) + G_SYN2*Y[V1] ;
F[V1] = (I_STIM-current[I_S2] - current[I_NA1] - current[I_K1] - current[I_K3] - current[I_L1] - current[I_S1])/CM;
theta = THETA_H; k1=K1_H; k2=K2_H; sigma1=SIGMA1_H; sigma2=SIGMA2_H;
F[H1] = beta(Y[V1],k1,sigma1)*(1.0 - Y[H1]) - alpha(Y[V1],k2,theta,sigma2)*Y[H1];
theta = THETA_N; k1=K1_N; k2=K2_N; sigma1=SIGMA1_N; sigma2=SIGMA2_N;
F[N1] = alpha(Y[V1],k1,theta,sigma1)*(1.0 - Y[N1]) - beta(Y[V1],k2,sigma2)*Y[N1];
theta = THETA_A; k1=K1_A; k2=K2_A; sigma1=SIGMA1_A; sigma2=SIGMA2_A;
F[A1] = alpha(Y[V1],k1,theta,sigma1)*(1.0 - Y[A1]) - beta(Y[V1],k2,sigma2)*Y[A1];
theta = THETA_M; k1=K1_M; k2=K2_M; sigma1=SIGMA1_M; sigma2=SIGMA2_M;
F[M1] = alpha(Y[V1],k1,theta,sigma1)*(1.0 - Y[M1]) - beta(Y[V1],k2,sigma2)*Y[M1];

return 0;
   }

void scan_(Y) 
double Y[N];
{FILE *fopen(),*sp;
int i;
sp = fopen("state.data","r");
for(i=0;i<N;i++) fscanf(sp,"%lf\n",&Y[i]);
fclose(sp);}

void dump_(Y) 
double Y[N];
{FILE *fopen(),*sp;
int i;
sp = fopen("end.data","w");
for(i=0;i<N;i++) fprintf(sp,"%.16f\n",Y[i]);
fclose(sp);}

int mas(n,amas,l)
        int *n;
        double *amas;
        int *l;
{return 0;}

int dummy(n,t,y,ydot)
        int *n;
        double *t;
        double *y;
        double *ydot;
{return 0;}

