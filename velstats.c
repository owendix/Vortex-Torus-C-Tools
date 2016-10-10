#include "ubiq.h"
#include <string>
#include <sstream>
#include <cstdio>
#include <cerrno>
using namespace std;
/* 3-torus analysis of velocity stats. Operates on restart files. 
Outputs to stdout
*************************************************************** 
If  params files were pre aug0412, need to include
absmax_fac: absmax2=absmax_fac*absmax_fac*r0*r0;
absmin_fac: absmin2=absmin_fac*absmin_fac*emax*emax*r0*r0;
limbiot_fac: limbiot=limbiot_fac*limbiot_fac*r0*r0;
After aug0412, these are part of the params files
Compile:
	Requires velstats.c linalg.c ubiq.h, linalg.h
	g++ -g velstats.c linalg.c -o vels.out

Input: OR CAN ADD abscomps BEFORE nskip	
	Pre-aug0412 trials:
	./vels.out ../../code/pfiles/parjul3112a.dat absmax_fac absmin_fac limbiot_fac ../../code/jul3112a nskip

	Post-aug0412 trials:
	./vels.out ../../code/pfiles/paraug0512a.dat ../../code/aug0512a nskip

Output:
	//cout: 
	//standard deviation of preceding quantity is printed on even columns
	//Can select if tan,norm,bnorm are actually just x,y,z (z=tan to v_{ns})
	Column    Quantity (Averages of)
        [1]     [time]
        2       tangent.z
        4       radius      
        6       beta
        8       |vlocal|
        10      vlocal-sp(or x)
        12      vlocal-norm(or y)
        14      vlocal-bnorm(or z)
        16      |v_nonlocal|
        18      v_nl-sp(or x)
        20      v_nl-norm(or y)
        22      v_nl-bnorm(or z)
        24      |v_fric|
        26      v_fric-sp(or x)
        28      v_fric-norm(or y)
        30      v_fric-bnorm(or z)
        32      |v_total|
        34      v_tot-sp(or x)
        36      v_tot-norm(or y)
        38      v_tot-bnorm(or z)
        40      Fsn.x
        42      Fsn.y
        44      Fsn.z 
        46      spacing (no std dev)
***********************************************************************
\rho(total)=145kg/m^3 (0<T<T\lambda)
\rho=\rho_n+\rho_s;
Barenghi,Donnelly,Vinen. J. Low Temp.Phys., 52,189 (1983):
@T=1.6K:\alpha=0.098,\alpha'=0.012,\rho=0.1452g/cm^3, \rho_s=0.1213 g/cm^3
***********************************************************************
Reads if NONLOCAL was on from params file, and sets global vble accordingly 
Assumes SPACCUTOFF was on: variable pt spacing and reconnect dist
*************************************************************************/

//IO function prototypes**************************************************
void readparams(string, int);
void readcores(point[], ptadmin[],ifstream&);
//rkf stuff
void nonlocalvel(point*, point*, double);
void derivs(point[],point[]);
point biot(point, point, point);
point localrz(point, point, point, point&, double&, double&);
//************************************************************************

/* Global variables for the module, entered in inout.c from params.dat */
int comps=1;//=0: use tan,norm,bnorm;=1:use x,y,z
int abscomps=0;//=1, use abs(components) in derivs, or =0
vortex startcore[NVORTMAX];
int outwrite;
double r0, absmax2, absmin2;

double tinit;	//really t: for a particular file
double vn, vs, emax, skew;
double a0, alpha;
double absmin_fac, absmax_fac, limbiot_fac;
int ilocal, ifric, inotan, ncycles;
int nonlocalon;

/* Other global variables */
double const1 = KAPPA / 4. / M_PI, const2;
double rhos=0.1213;//g/cm^3: if \alpha\approx0.1 (T=1.6K)
/* Global variables */
int nvort, totalpts;
double radius[NMAX];
point core[NMAX];
ptadmin periph[NMAX];

#ifdef VELSTRUCT
char velstrfile[39];
#endif

int main(int argc, char* argv[]) 
{
	/*argv order
	IF params file is later than aug 04 2012 (argc=3)
	./velstats ../../code/pfiles/parfilename.dat ../../code/filename nskip
	IF params file is earlier than aug 04 2012 (argc=5)
	./velstats ../../code/pfiles/parfilename.dat absmax_fac absmin_fac limbiot_fac ../../code/filename nskip */
	string infilepre, infilemid, infilebase, infile, paramsfile, strnum;
	stringstream strhandle;
	ifstream infhandle;
	double limbiot;
	int loc, nskip, done, num, numlow;
	point  nonloc[NMAX];

	if (argc==4 || argc==7){	//argv[argc] is null pointer
		nskip=atoi(argv[argc-1]);
		infilepre=(string)argv[argc-2];//no . or final /
		paramsfile=(string)argv[1];	//full params file path in pfiles/
	}else if (argc==5 || argc==8){
		nskip=atoi(argv[argc-1]);
		abscomps=atoi(argv[argc-2]);//=1, use abs(components) in derivs, or =0
		infilepre=(string)argv[argc-3];//no . or final /
		paramsfile=(string)argv[1];	//full params file path in pfiles/
	}else{
		cerr << "If params file is pre-aug0412, Input:\n";
		cerr << "  paramsfile absmax_fac absmin_fac limbiot_fac restartfilepre";
		cerr << " nskip\n";
		cerr << "If params file is post-aug0412, Input:\n";
		cerr << "  paramsfile restartfilepre nskip\n";
		cerr << "argc=" << argc << endl;
		cerr << "argv[argc-1]=" << argv[argc-1] << endl;
		exit(1);
	}
	if (nskip<1)
		nskip=1;
	if (abscomps!=1)
		abscomps=0;
	cerr << "infilepre = " << infilepre.c_str() << endl;
	cerr << "paramsfile = " << paramsfile.c_str() << endl;
	cerr << "nskip = " << nskip << endl;
	cerr << "abscomps = " << abscomps << endl;

	//Concatenate filename to path
    if ((loc=infilepre.find_last_of("/"))==string::npos)
        perror("Error in find_last_of()\n");
    infilemid=infilepre;
	//loc is position (last / before filename), 9 is length of "/mon3112a"
	//assign replaces old content with this "/mon3112a"
    infilemid.assign(infilepre,loc+1,8);//filename with leading /
    infilebase=(infilepre+"/")+infilemid; 
		//concat with num to make infile, inside while loop

	/*I want readparams to read in r0, NONLOCAL, absmax_fac, absmin_fac,
	limbiot_fac. Some I may have to input myself with argv. If 
	argc certain size, override
	*/
	//Read params from a specific file, included in argv
	//This params file is common to all in the list of filenames (restarts)
	readparams(paramsfile, argc);
	const2=exp(0.25)*a0;//for derivs()
	
	//read absmax_fac, absmin_fac, limbiot_fac from argv, or from params 
	//if params file is later than aug 04 2012
	if (argc==7 || argc==8){//strtod more robust than atof (sets errno)
		//sets to 0.0 upon failure
		absmax_fac=strtod(argv[2],NULL);
		absmin_fac=strtod(argv[3],NULL);
		limbiot_fac=strtod(argv[4],NULL);
		if (absmax_fac==0 || absmin_fac==0 || limbiot_fac<0 || errno){
			cerr << "Problem assigning vals from argv[], ";
			cerr << strerror(errno) << endl;
			exit(1);
		}
		cerr << "absmax=" << absmax_fac << "*r0" << endl;
		cerr << "absmin=" << absmin_fac << "*emax*r0" << endl;
		cerr << "sqrt(limbiot)=" << limbiot_fac << "*r0" << endl;
		cerr << "r0=" << r0 << ", vn=" << vn << ", vs=";
		cerr << vs << ", emax=" << emax << endl;
		if (nonlocalon)
			cerr << "NONLOCAL on" << endl;
		else
			cerr << "NONLOCAL off" << endl;
	}
	absmax2=absmax_fac*absmax_fac*r0*r0;
	absmin2=absmin_fac*absmin_fac*emax*emax*r0*r0;
	limbiot=limbiot_fac*limbiot_fac*r0*r0;

	#ifdef VELSTRUCT
    sprintf(velstrfile,"%s.%s.%s","velstruct",infilemid.c_str(),"dat");
    #endif
	//while opening the filename file was successful
	num=1;
	numlow=num;
	cerr << "Column    Quantity (Averages of)\n";
    cerr << "[1]    [time]\n";
    cerr << "2      tangent.z\n";
    cerr << "4      radius\n";
    cerr << "6      beta\n";
    cerr << "8      |vlocal|\n";
    cerr << "10     vlocal-tan\n";
    cerr << "12     vlocal-norm\n";
    cerr << "14     vlocal-bnorm\n";
    cerr << "16     |v_nonlocal|\n";
    cerr << "18     v_nl-tan\n";
    cerr << "20     v_nl-norm\n";
    cerr << "22     v_nl-bnorm\n";
    cerr << "24     |v_fric|\n";
    cerr << "26     v_fric-tan\n";
    cerr << "28     v_fric-norm\n";
    cerr << "30     v_fric-bnorm\n";
    cerr << "32     |v_total|\n";
    cerr << "34     v_tot-tan\n";
    cerr << "36     v_tot-norm\n";
    cerr << "38     v_tot-bnorm\n";
    cerr << "40     Fsn.x\n";
    cerr << "42     Fsn.y\n";
    cerr << "44     Fsn.z\n";
    cerr << "46     spacing (no std dev)\n";
	while (1){
		strhandle << num;
		strnum=strhandle.str();
		infile=infilebase+"."+strnum;
		strhandle.str("");//clears it for next use
		infhandle.open(infile.c_str());
		if (!infhandle.good()){//would exit upon eofbit being set
			//Not necessarily an error, but I'll be > stdout to file
			cerr << "Unable to open file " << infile.c_str() << endl;
			exit(0);
		}
	
		//Read file from "restart" file: one of a list included in argv
		//may not be able to pass ifstream handle
		readcores(core, periph, infhandle);
	
		//Calculate the nonlocal velocity field
		nonlocalvel(core, nonloc,limbiot);
		//Call derivs which prints either velstats (stdout) 
		//or velstats (stdout) AND velstruct (velstruct.filename.dat)
		derivs(core,nonloc);

		infhandle.close();
		num+=nskip;
	}	
	
	return (0);
}

//how does nonloc[i] get passed/defined/in derivs()...?
void nonlocalvel(point* y, point* nonloc, const double limbiot){
	static int i, j, ivort;
	static double dist2, skiplength, jnkdb;//full limbiot=3*r0*r0/4
	//double limbiot=3*r0*r0/4; Set in argv
	static point errsum, tmppt, ypt;
	static const double maxlfactor=0.2, minlfactor=0.0833333333333;//only maxlfactor needed
	//the values are taken from reconnect.c
	//y is the core value

	for (i=0; i<totalpts; i++){//i is point calculating nonlocal ON, by pts j
		nonloc[i].x = 0; nonloc[i].y = 0; nonloc[i].z = 0;
		//Error for compensated summation: Higham (2002) pg 84
		errsum.x=0.0; errsum.y=0.0; errsum.z=0.0;
		for (ivort=0; ivort<nvort; ivort++){
			j=startcore[ivort].start;
			do{
				if (j!=i && periph[j].ip !=i){//j is source of nonlocal on i
					if ((dist2=dot(y[j]-y[i],y[j]-y[i]))<limbiot){
						tmppt = nonloc[i];
						ypt = errsum + biot(y[j],y[periph[j].ip],y[i]);
						nonloc[i] = tmppt + ypt;
						errsum = (tmppt + (-1.0*nonloc[i])) + ypt;	//P.B.C.s
						//w/ compensated sum method
					}else{
						skiplength=sqrt(dist2)-sqrt(limbiot);
						while (j!=startcore[ivort].end && skiplength > 0){
							j=periph[j].ip;
							jnkdb=maxlfactor*radius[j];
							if (jnkdb*jnkdb>absmax2)
								jnkdb=sqrt(absmax2);
							else if (jnkdb*jnkdb<absmin2)
								jnkdb=sqrt(absmin2);
							skiplength -= jnkdb;
						}
					}
				}
				j=periph[j].ip;
			}while (j!=startcore[ivort].start);
		}//end for (ivort)
		nonloc[i] = nonloc[i] + errsum;	//Final correction: comp sum
	}
}

point biot(point stseg, point endseg, point testpt)
{
      static point l1, l2, l3, deldot, dir;
      static double dot11, dot12, dot22, dot33, dot23, factor1, factor2;

      l1 = testpt - stseg;
      l2 = endseg - stseg;
      l3 = testpt - endseg;      
      
      dot11 = dot(l1, l1);
      dot12 = dot(l1, l2);
      dot22 = dot(l2, l2);
      dot23 = dot12 - dot22;
      dot33 = dot(l3, l3);

      factor1 = dot11*dot22-dot12*dot12;

      if (fabs(factor1/dot11/dot22) > 1e-14){
        factor2 = dot23/sqrt(dot33) - dot12/sqrt(dot11);
        deldot = (factor2 / factor1) * cross(l1, l2);
      }else{
        factor1 = dot22*dot33-dot23*dot23;
        if  (fabs(factor1/dot22/dot33) > 1e-14){
            factor2 = dot23/sqrt(dot33) - dot12/sqrt(dot11);
            deldot = (factor2 / factor1) * cross(l3, l2);
        }else{
            if (dot11*dot33 > 1e-12 && dot12*dot23 > 0.0) {
                factor2 = fabs(dot11-dot33)/dot11/dot33/sqrt(dot22);
                deldot = (factor2 / 2) * cross(l1, l2);
            }else if (dot12*dot23 <= 0.0){
				factor2 = dot23/sqrt(dot33) - dot12/sqrt(dot11);
				dir = cross(l1, l2);
				factor1 = dot(dir, dir);
				deldot = (factor2 / factor1) * dir;
			}else{
                deldot.x = 0; deldot.y = 0; deldot.z = 0;
            }
        }
      }
      return deldot;
}

/* Evaluates the derivatives (RHS) of the system of ODE's x' = f(t,x) */
void derivs(point coretemp[], point nonloc[]) {
	static point f[NMAX];
	static int i, j;
	static point sp; /* unit vector approximating tangent to vortex line */
	static point tmp;
	static double beta, space, spaceprev, weight, sumweight, tempwt, Rwt;
	static point fD;
	static const double cunity=1.1;
	//const2 set in main()
	//Print, for each vortex, the avg and std dev in several
    //characteristics of the velocity terms
    static struct statsdb{
        double m;   //Like avg
        double q;   //Like variance
    }stats[22];
	static double num[22], dif;
	static point temppt, norm, bnorm;
	static int k, ivort;
	//initialize stats
	for (k=0;k<22;k++){
		stats[k].m=0.0;
		stats[k].q=0.0;
	}

	sumweight=0.0;//line length
#ifdef VELSTRUCT
    //prints, for each point on a vortex, radius, spacing relative
    //to absmax2, self-induced velocity, nonlocal velocity, 
    //and the tangential component of nonlocal velocity values around
    //the vortex. No averaging, lots of printing: expensive
    //point binorm;
    ofstream fvstruct(velstrfile, ios_base::out | ios_base::app);
#endif
	//loops do not need to be done this way. Could be done, in i=0->totalpts
	//order. This is how it was done in rkf.c, so I don't have to change
	//i->j, etc.
	for (ivort=0;ivort<nvort;ivort++){
		i=0;
		j=startcore[ivort].start;
		do {//using this so it ends when j=startcore[ivort].start again
			/* Initialize the derivative vector */
			f[i].x = 0.;
			f[i].y = 0.;
			f[i].z = 0.;
	
			/* Find the local self-induced field */
			if (ilocal+ifric+inotan !=0){
				f[i] = localrz(coretemp[periph[j].im], coretemp[j],
					coretemp[periph[j].ip], sp, radius[j], beta);
				if (ilocal !=0){
					if (!nonlocalon){
						//beta = const1*log(r0/a0);
						beta = cunity*const1*log(radius[j]/a0);
					}
					f[i] = beta*f[i];
				}
			}
			
			if (!comps){//don't use R^3 components
				bnorm=f[i]/num[3];
				norm=cross(bnorm,sp);
			}
		
    	    temppt=coretemp[periph[j].ip]-coretemp[j];
			space=sqrt(dot(temppt,temppt));//spacing: needed for weighted avg
			//establish weighting factor
			if (i==0){
				temppt=coretemp[j]-coretemp[periph[j].im];
				spaceprev=sqrt(dot(temppt,temppt));
			}
			weight=0.5*(space+spaceprev);
			num[0] = sp.z;//velocity direction
			num[1] = radius[j];	//radius
			num[2] = beta;		//beta
			num[3]=sqrt(dot(f[i],f[i]));//local velocity
			num[4]=dot(f[i],sp);
			num[5]=dot(f[i],norm);
			num[6]=dot(f[i],bnorm);
    		if (nonlocalon){
    			num[7]=const1*sqrt(dot(nonloc[j],nonloc[j]));//v_nl
				if (!comps){//don't use R^3 components
    				num[8]=const1*dot(nonloc[j],sp);//v_nl||tangent
					num[9]=const1*dot(nonloc[j],norm);//v_nl||normal
					num[10]=const1*dot(nonloc[j],bnorm);//v_nl||binorm
				}else{
					num[8]=const1*nonloc[j].x;
					num[9]=const1*nonloc[j].y;
					num[10]=const1*nonloc[j].z;//tangent to v_{ns}
				}
    		}else{
				num[7]=0.0;
    			num[8]=0.0;
    			num[9]=0.0;
    			num[10]=0.0;
    		}
			
			tempwt=weight+sumweight;//do this once, now 
    		for (k=0;k<11;k++){
				if (k>3 && k!=7 && abscomps)
    		    	dif = abs(num[k]) - stats[k].m;
				else
    		    	dif = num[k] - stats[k].m;
				Rwt=dif*weight/tempwt;
    		    stats[k].m += Rwt;
    		   	stats[k].q += sumweight*dif*Rwt;
    		}
			//sumweight=tempwt;do this at the end
			if (nonlocalon){
				//Add nonlocal field
				f[i] = f[i] + const1*nonloc[j];
			}
			/* Add imposed superfluid flow velocity */
		#ifdef VSHELIX
			//Add helical velocity field
			tmp.x = -sin(2*M_PI*ncycles*coretemp[i].z/r0);
			tmp.y = cos(2*M_PI*ncycles*coretemp[i].z/r0);
			tmp.z = skew;
			f[i] = f[i] + vs/sqrt(1+skew*skew)*tmp;
		#else
			#ifndef VNHELIX
			//Constant velocity field
			tmp.x =0; tmp.y =0; tmp.z =1;
			f[i] = f[i] + vs*tmp;
			#endif
		#endif
			
			// Add normal component for friction term
			temppt=f[i];	//Save prior to adding friction
		#ifdef VNHELIX
			//Add helical velocity field
			tmp.x = -sin(2*M_PI*ncycles*coretemp[i].z/r0);
			tmp.y = cos(2*M_PI*ncycles*coretemp[i].z/r0);
			tmp.z = skew;
			/* Find the friction-induced field */
			if (ifric != 0)
				f[i]=f[i]+alpha*cross(sp,(vn*tmp/sqrt(1+skew*skew)+(-1)*f[i]));
		#else	
			#ifdef VSHELIX
			tmp.x = 0; tmp.y = 0; tmp.z = 0;
			/* Find the friction-induced field */
			if (ifric != 0)
				f[i] = f[i] + alpha*cross(sp, (-1)*f[i]);
			#else
			tmp.x = 0; tmp.y = 0; tmp.z = 1;
			/* Find the friction-induced field */
			temppt=f[i];	//Save prior to adding friction
			if (ifric != 0)
				f[i] = f[i] + alpha*cross(sp,(vn*tmp+(-1)*f[i]));
			#endif
		#endif
			temppt=f[i]+(-1*temppt);//get friction term;
			//get drag force on each vortex seg
			fD=-rhos*KAPPA*cross(sp,temppt);
			num[11]=sqrt(dot(temppt,temppt));//friction magnitude
			if (!comps){//if not R^3 components
				num[12]=dot(temppt,sp);
				num[13]=dot(temppt,norm);
				num[14]=dot(temppt,bnorm);
			}else{
				//friction
				num[12]=temppt.x;
				num[13]=temppt.y;
				num[14]=temppt.z;
			}
			
			/* Reject tangential velocity */
			if (inotan != 0)
				f[i]=f[i]+(-1.0*dot(sp,f[i]))/dot(sp,sp)*sp;
			//total velocity and its components
        	num[15] = sqrt(dot(f[i],f[i]));
			if (!comps){//if not R^3 components
				num[16] = dot(f[i],sp);
				num[17] = dot(f[i],norm);
				num[18] = dot(f[i],bnorm);
			}else{
				num[16]=f[i].x;
				num[17]=f[i].y;
				num[18]=f[i].z;
			}
			//Fsn mutual friction force density
			num[19]=fD.x;
			num[20]=fD.y;
			num[21]=fD.z;
    	    for (k=11;k<22;k++){
				if (k>11 && k!=15 && k<19 && abscomps)//only components
					dif=abs(num[k])-stats[k].m;
				else
					dif=num[k]-stats[k].m;
				Rwt=dif*weight/tempwt;
    	        stats[k].m += Rwt;
    		   	stats[k].q += sumweight*dif*Rwt;
    	    }
			sumweight=tempwt;//wait until now, the end, to do this
		#ifdef VELSTRUCT//Doesn't use abs(vel components)
			//Printing to VELSTRUCT starts here
	        fvstruct << "(ivort,i,j)=(" << ivort << ",";
   	    	fvstruct << i << "," << j << ") ";
   	    	fvstruct << "(x,y,z)=(" << coretemp[j].x << ",";
			fvstruct << coretemp[j].y << "," << coretemp[j].z << ") ";
			for (k=0;k<22;k++){
					fvstruct << num[k];
				if (k!=21)
					fvstruct << " ";
				else
					fvstruct << endl;
			}
    	#endif
			i++;
			j=periph[j].ip;
			spaceprev=space;
		} while (j!=startcore[ivort].start);
		//end: for all points within each vortex
	}//end: for all vortices
	
	//for Fsn only, all others are averaged over line length
    //Checking velocity statistics
    //avg = stats[i].m
	//var = stats[i].q/sumweight;
	//cout: (fvstruct prints only values, one line each quantity)	
	//standard deviation of preceding quantity is printed on even columns
	/*Column    Quantity (Averages of)
        [1]     [time]
        2       tangent.z
        4       radius      
        6       beta
        8       |vlocal|
        10      vlocal-sp(or x)
        12      vlocal-norm(or y)
        14      vlocal-bnorm(or z)
        16      |v_nonlocal|
        18      v_nl-sp(or x)
        20      v_nl-norm(or y)
        22      v_nl-bnorm(or z)
        24      |v_fric|
        26      v_fric-sp(or x)
        28      v_fric-norm(or y)
        30      v_fric-bnorm(or z)
        32      |v_total|
        34      v_tot-sp(or x)
        36      v_tot-norm(or y)
        38      v_tot-bnorm(or z)
        40      Fsn.x
        42      Fsn.y
        44      Fsn.z 
        46      spacing (no std dev)
    */
	//Output time, from ifname (tinit)
	cout.precision(8);
	cout << tinit << " ";
	for (k=0;k<22;k++){
		if (k<19){
			cout << stats[k].m << " " << sqrt(stats[k].q/sumweight) << " ";
		}else{//var for Fsn has factor of sumweight/
			cout << (stats[k].m*sumweight/(r0*r0*r0)) << " ";
			cout << sqrt(stats[k].q/(r0*r0*r0)) << " ";
		}
	}
	cout << sumweight/(totalpts-1) << endl;
	
#ifdef VELSTRUCT
	fvstruct.close();
#endif
	return;
}

point localrz(point start, point mid, point end, point& tangent, double& rad,
		double& prefac)
{
	point lp, lm, vel, curv;
	double lpdotp, lmdotm, lpdotm, norm;

	point sumpt, combopt;
	double denom, dp, dm, tmpdb, tmpdb1, tmpdb2;
	
	lp = end - mid;		//Subtracting point structures contains PBC condition
	lm = mid - start;

	lpdotm = dot(lp, lm);
	lpdotp = dot(lp, lp);
	lmdotm = dot(lm, lm);
	//New definition of tangent vector with error 
	//2nd order in neighboring arclength
	//Previous definition was 1st order
	tangent = lmdotm*lp + lpdotp*lm;
	tangent = tangent/sqrt(dot(tangent,tangent));
	//same as newpt() and newfakept()
	if ((tmpdb=acos(lpdotm/sqrt(lpdotp*lmdotm)))>0.2){
        denom = lpdotp*lmdotm - lpdotm*lpdotm;
        sumpt=lp+lm;
        dp = lmdotm*dot(lp,sumpt);
        dm = lpdotp*dot(lm,sumpt);
        combopt=dp*lp+(-1)*dm*lm;
        curv = 2*denom/dot(combopt,combopt)*combopt;
    }else{
        if (tmpdb<1e-15){
            curv.x=0.0; curv.y=0.0; curv.z=0.0;
        }else{
            /*Linear in s_{+-}, derived by Aarts 3.9 & 3.10: 
            YIELDED BAD RESULTS AT HIGH ANGLES (truncation error), 
            EXCELLENT AT LOW ANGLES (low rounding error)*/
            curv=2*(lp/sqrt(lpdotp)+(-1)*lm/sqrt(lmdotm))
                /sqrt(lpdotp+lmdotm+2*lpdotm);
        }
    }

	norm = dot(curv, curv);
	if (norm > 1e-20){
		vel = cross(tangent, curv); 
		rad = 1.0/sqrt(norm);
	}else{
		vel.x = 0; vel.y = 0; vel.z = 0; 
		rad = 1e10; //radius set here too
	}

	prefac = const1*log(2.*sqrt(sqrt(lpdotp*lmdotm))/const2);

	if (isnan(vel.x) || isnan(tangent.x) || isnan(rad)){
		cerr.precision(17);
		cerr << "localrz contains nan" << endl;
        cerr << "totalpts=" << totalpts << " nvort=" << nvort << endl;
        cerr << "start=" << start << endl;
        cerr << "mid=" << mid << endl;
        cerr << "end=" << end << endl;
		cerr << "lpdotp=" << lpdotp << " lpdotm=" << lpdotm;
		cerr << " lmdotm=" << lmdotm << endl;
		cerr << "rad=" << rad << " norm=" << norm << endl;
		cerr << "tangent=" << tangent << endl;
		cerr << "curv=" << curv << endl;
		exit(1);
	}

	return vel;
}

/* The point of dummy is so the compiler correctly identifies the class T. 
 This routine is first because I haven't figured out how to prototype template
 classes properly. */
template<class T> T readvar(T dummy, ifstream& ffrom, string varname)
{	//passes ffrom by reference	
	/*params must be in same order as readvar(), but some in 
	params.dat can be ignored in readvar() list.
	When using old params.dat or not, must be able to skip useless
	quantities.
	When using old params.dat need to skip quantities to get to NONLOCAL
	When using new params.dat just need to skip intermediate quantities*/
	string phrase; T thenum;
	//Not used for setting nonlocalon when "NONLOCAL on" detected

	//ifstream operator>> is the extraction operator
	while(ffrom.good()){//allows for skipping of varnames I don't want
					//set if badbit, failbit or eofbit
		ffrom >> phrase;	//ends at valid whitespace
		while (phrase[0]=='#') {
			if (phrase[1]=='#'){//reached the macros added to params.dat
				//go back (-) what was just put down (gcount) from "cur"(rent) 
				//position
				ffrom.seekg((streamoff)(-(int)ffrom.gcount()),ios_base::cur);
				return -1;//no parameters have negative val, use as error
			}else if (!ffrom.good()){
				if (ffrom.eof()){
					cerr << "Reached EOF in params file, in readvar\n";
				}else{
					cerr << "Either badbit or failbit set, in readvar\n";
				}
				exit(1);
			}
			getline(ffrom, phrase);//gets rest of line until \n
			ffrom >> phrase;//get beginning of new line until whitespace
		}
		if (phrase==varname) {//already took varname (string)
			ffrom >> phrase;	//then takes '=' sign
			ffrom >> thenum;	//then takes thenum
			break;
		}
	}
	if (!ffrom.good()){
		cerr << "Problem in params readvar finding "<< varname << endl;
		exit(1);
	}

	return thenum;
}

void readparams(string params, int argc) {
	int i;
	size_t found;
	string line;
	ifstream from(params.c_str());
	double cutoff_fac;//not needed as global variable in this program

	//in this particular version, variables must be in order 
	//but there can be variables in the params file that aren't 
	//present here
	vn = readvar(vn, from, (string) "vn");
	vs = readvar(vs, from, (string) "vs");
	ilocal = readvar(ilocal, from, (string) "ilocal");
	ifric = readvar(ifric, from, (string) "ifric");
	inotan = readvar(inotan, from, (string) "inotan");
	a0 = readvar(a0, from, (string) "a0");
	alpha = readvar(alpha, from, (string) "alpha");
	emax = readvar(emax, from, (string) "emax");
	ncycles = readvar(ncycles, from, (string) "ncycles");
	skew = readvar(skew, from, (string) "skew");
	r0 = readvar(r0, from, (string) "r0");
	//If others =-1, r0==-1 too:
	if (r0==-1){//never set/read r0 from params file
		cerr << "params file in readvar bad, r0=-1" << endl;
		exit(1);
	}
	if (argc==4 || argc==5){//Didn't include these in argv[] in main()
		absmax_fac = readvar(absmax_fac, from, (string)"absmax_fac");
		absmin_fac = readvar(absmin_fac, from, (string)"absmin_fac");
		cutoff_fac = readvar(cutoff_fac, from, (string)"cutoff_fac");
		limbiot_fac = readvar(limbiot_fac, from, (string)"limbiot_fac");
	}
	//look for NONLOCAL on, if so, set nonlocalon=1
	//ifstream operator>> is the extraction operator
	i=0;//if last one should end long before ##, if 
	while(from.good()){//found r0 properly, then exited
		getline(from,line);//retrieves entire line, unlike 
		if (line[1]=='#')//operator>> which ends at valid whitespace
			if(i==0)//2nd occurrence of ## at beginning of line
				i++;
			else
				break;
	}//#defined macros at end of params file
	nonlocalon=0;
	if (line[1]=='#'){//second occurrence of ## at beginning of line
		found=line.find("NONLOCAL on");//finds first occurrence of string
		if (found!=string::npos)//returned if string not found
			nonlocalon=1;
	}

	return;
}

void readcores(point core[], ptadmin periph[], ifstream& from) {
	int i, j, jmin, ivort;
	double hinit;	//normally global but useless here

	from >> totalpts; //totalpts; must not exceed NMAX in ubiq.h
	from >> nvort; //nvort; it must not exceed NVORTMAX in ubiq.h

	vortex * vlims = new vortex [nvort];
	// restart file overrides tinit from params.dat
	from >> tinit;
	from >> hinit;

	for (i = 0; i < nvort; i++) {
		from >> vlims[i].start;
		from >> vlims[i].end;
		from >> startcore[i].term;
	}

	for (ivort = 0; ivort < nvort; ivort++) {
		for (i = vlims[ivort].start; i != vlims[ivort].end; i++) {
			from >> j;
			from >> periph[j].recon;
			from >> core[j];
			if (i == vlims[ivort].start) {
				startcore[ivort].start = j;
			} else {
				periph[j].im = jmin;
				periph[jmin].ip = j;
			}
			jmin = j;
		}
		from >> j;
		from >> periph[j].recon;
		from >> core[j];
		startcore[ivort].end = j;
		periph[j].im = jmin;
		periph[jmin].ip = j;
		periph[j].ip = startcore[ivort].start;
		periph[startcore[ivort].start].im = j;
	}
	delete [] vlims;
	return;
}

