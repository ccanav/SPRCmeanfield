#include <stdio.h>
#include <math.h>
#define TAU1 2.0 //2.3
#define TAU2 0.3 //0.484
#define TAURISE
#define P 4.0   
main()
{
int i;
double gsyn, f,c1,c2,arg,ts=0,x=0;
c1= 1.0/(1.0*TAU1);
c2= 1.0/(1.0*TAU2);
arg= log(c2/c1)/(c1-c2);
arg=(TAU2*TAU1)*log(TAU1/TAU2)/(TAU1-TAU2);
f = 1.0/(exp(c1*arg)-exp(c2*arg));
f = 1.0/(exp(-arg/TAU1)-exp(-arg/TAU2));
for(i=0;i<1600;i++)
{ x=0.01*i;
	if(x > ts +P-0.01) ts=ts+P; 
  gsyn=exp((ts-x)/TAU1)*(1.0/(1-exp(-P/TAU1))) -exp((ts-x)/TAU2)*(1.0/(1.0-exp(-P/TAU2)));
  gsyn=f*gsyn;
 printf("%f %f\n",x,gsyn);
}
}
