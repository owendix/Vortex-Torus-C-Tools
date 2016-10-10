#include "ubiq.h"
#include <string>
#include <stdio.h>

double radcurv(point, point, point);
void readcores(point[], char[]);

vortex startcore[NVORTMAX];
point core[NMAX];
ptadmin periph[NMAX];
int nvort, totalpts;
double radius[NMAX], tinit, r0=0.1;
string filename;
char wfile[20];  

main ()
{
    int i, counter=0, filenum, lastfile=10, firstfile=1;
    double avgrad=0;

    filename = "may15b";
	for (filenum=firstfile; filenum<lastfile; filenum++) {
	  sprintf(wfile,"%s.%d",filename.c_str(), filenum);      
      readcores(core, wfile);   
      for (i=0; i < totalpts; i++)     
      {
         periph[i].ip = i + 1;  periph[i].im = i - 1;
         periph[i].recon = 1;     /* initialize reconnection parameter */
      }

      /* Now make rings close up. */
      for (i=0; i < nvort; i++) {
         periph[startcore[i].start].im = startcore[i].end;
         periph[startcore[i].end].ip = startcore[i].start;
      }

      counter = 0;
      for (i=0; i < totalpts; i++) {
         radius[i] = radcurv(core[periph[i].im], core[i], core[periph[i].ip]);
         if (radius[i] > -1) {
            counter += 1;
            avgrad += radius[i];
         }
      }
      avgrad = avgrad/counter;
      cout <<tinit<<"	"<< avgrad <<"  "<<counter<<"\n";
	}
}

double radcurv(point start, point mid, point end)
{
    point lp, lm, vel, curv, tangent;
    double lpdotp, lmdotm, lpdotm, x, norm, invrad, rad;
    
    lp = end - mid;
    lm = mid - start;
 
    lpdotm = dot(lp, lm);
    lpdotp = dot(lp, lp);
    lmdotm = dot(lm, lm);
    x = sqrt(lpdotp + lmdotm + 2*lpdotm);
    tangent = (lp + lm)/x;
    
    curv = lp - dot(lp, tangent) * tangent;
    norm = dot(curv, curv);
    if (norm > 1e-24) {
        curv = -1.*curv/sqrt(norm);
	invrad = 2/x * sqrt(1 - lpdotm*lpdotm/lpdotp/lmdotm);
	if (invrad > 1/r0)
		invrad = sqrt(invrad*invrad-1/r0/r0);
	else invrad = 0;        // indicates infinite radius, i.e. straight line
    } 
    else {
	invrad = 0;
    }
	
    return invrad;    
}


void readcores(point core[], char writefile[])
{
    int i;
    ifstream from(writefile);

    from >> totalpts;  //totalpts; must not exceed NMAX in ubiq.h
    from >> nvort;  //nvort; it must not exceed NVORTMAX in ubiq.h

    // restart file overrides tinit from params.dat 
    from >> tinit;

    for (i=0; i< nvort; i++) {
      from >> startcore[i].start;
      from >> startcore[i].end;
      from >> startcore[i].term;
    }
    for (i=0; i < totalpts; i++)
      from >> core[i];
    return;
}
