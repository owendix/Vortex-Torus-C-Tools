#include "ubiq.h"
#include <sstream>
/**********************************************
Measures anisotropy of vortex tangles in trial. See 
bottom for what is printed. If isotropic, 
ipar = iperp1 = iperp2 = 2/3 and ibinpar=0. If lying 
within plane normal to velocity field, ipar=1, 
iperp1 =? iperp2 =? 1/2 and ibinpar is more complicated. 
The following relation should hold: ipar/2 + iperp1,2 = 1. 
***********************************************
Compile with linalg.c
Input: g++ -g anisotropy.c linalg.c -o isot.out
		./isot.out ../../code/mar1011a 10
-> use trial mar1011a only analyzing every 10th file
Output: (cout) 
				Column	Data
				1 		time
				2		ldens (optional)
				3		ipar
				4		iperp1
				5		iperp2
				6		ibin (forced parallel 2 v_{ns})
				7		ibinpar (not parallel)
				8		ibin/ibinpar
				9		ipar/2+iperp1
				10		ipar/2+iperp2
************************************************/
using namespace std;

void localrz(point, point, point, point&, point&);

double r0=0.01;

int main(int argc, char* argv[]){

	int printdist=1;	//=1, print dist (c.f. ldens), =0, don't print
	
	int num, numlow, nskip, done, totalpts, nvort;
	int line, i, j, ivort, tmpint1, tmpint2, vlims[NVORTMAX];
	double a0=1.3e-8, volinv;
	double dist, ipar, iperp1, iperp2;
	double time, tmpdb1, tmpdb2, tmpdist;
	point t, b, r1st, r2nd, r, rm, rmm;
	point ibin, ibinpar, v, vp1, vp2, tmppt;
	string ifnamemid, ifnamepre, ifnamebase, ifname;
	string strnum;
	stringstream strhandle;
	ifstream ifhandle;
	int loc;

	if (argc==4){
		nskip=atoi(argv[argc-1]);
		nskip=(nskip<1)?1:nskip;
		r0=atof(argv[argc-2]);
		if (r0<=0.0){
			cerr << "Input r0 must be positive\n";
			exit(1);
		}
		ifnamepre=(string)argv[1];//no . or last /
	}else if (argc==3){	//argv[argc] is null pointer
		nskip=atoi(argv[argc-1]);
		if (nskip<1)
			nskip=1;
		ifnamepre=(string)argv[1];//no . or last /
	}else if (argc==2){
		nskip=1;
		ifnamepre=(string)argv[1];//no . or last /
	}else{
		cerr << "Input requires 2 arguments: ./a.out filename nskip\n";
		exit(1);
	}

	//Tell me the r0 value being used: do it on standard error
	cerr << "The r0 value is: " << r0 << endl;
	
	volinv=1/(r0*r0*r0);
	//directory where files are located	
	//ifnamepre="../../code/filename";

	//Concatenate filename to path
	if ((loc=ifnamepre.find_last_of("/"))==string::npos)
		perror("Error in find_last_of()\n");
	ifnamemid=ifnamepre;
	ifnamemid.assign(ifnamepre,loc,9);//filename with leading /

	ifnamebase=ifnamepre+ifnamemid;	//concat with num to make ifname

	num=1;
	numlow=num;
	done=0;
	//calculate applied velocity field at rm (unit vec)
	v.x=0; v.y=0; v.z=1;
	//calculate perpendicular directions to applied field
	vp1.x=1; vp1.y=0; vp1.z=0;
	vp2.x=0; vp2.y=1; vp2.z=0;
	while (!done){
		dist=0.0;
		ipar=0.0;
		iperp1=0.0;
		iperp2=0.0;
		ibinpar.x=0.0;ibinpar.y=0.0;ibinpar.z=0.0;
		strhandle << num;
		strnum=strhandle.str();
		ifname=ifnamebase+"."+strnum;
		strhandle.str("");	//clears it for next use
		ifhandle.open(ifname.c_str());
		if (!ifhandle.good()){
			done=1;
		}else{
			ifhandle >> totalpts;
			ifhandle >> nvort;
			ifhandle >> time;
			ifhandle >> tmpdb1;
			//get vortex limits
			for (ivort=0;ivort<nvort;ivort++){
				ifhandle >> tmpint1;
				ifhandle >> tmpint2;
				vlims[ivort]=tmpint2-tmpint1+1;	//num of pts per vortex
				ifhandle >> tmpint1;
			}
			//start going through pts
			line=0;
			for (ivort=0; ivort<nvort; ivort++){
				//first point
				ifhandle >> tmpint1; ifhandle >> tmpint1;
				ifhandle >> r1st;	//1st pt within vortex
				rm=r1st;
				line++;
				//r1st is actually i=0 pt, r2nd is i=1 pt
				for (i=1;i<vlims[ivort];i++){	//i: pt number within vortex
					ifhandle >> tmpint1; ifhandle >> tmpint1;
					ifhandle >> r;
					if (i==1)
						r2nd=r;	//2nd pt within vortex
					tmppt=r-rm;	//diff btw pt and prev pt
					tmpdist=sqrt(dot(tmppt,tmppt));
					dist += tmpdist;
					
					//compute anisotropy measures
					if (i>1){ //r, rm, rmm all stored by this iteration
						//associate rm with distance rm->r, calculate
						//tangent and binormal vectors at rm
						localrz(rmm, rm, r, t, b);	//sets tangent, t, and 
													//binormal, b
						//ipar
						tmpdb1=dot(t,v);
						ipar += (1-tmpdb1)*(1+tmpdb1)*tmpdist;
						//iperp1
						tmpdb1=dot(t,vp1);
						iperp1 += (1-tmpdb1)*(1+tmpdb1)*tmpdist;	
						//iperp2
						tmpdb1=dot(t,vp2);
						iperp2 += (1-tmpdb1)*(1+tmpdb1)*tmpdist;	
						//ibinpar
						ibinpar =  ibinpar + b*tmpdist;	
					}
		
					//circulate old pts
					rmm=rm;
					rm=r;
					
					line++;	//cumulative pt number within tangle
				}
				//Complete ring: get last distance rlast->r1st
				//and calculate anisotropy at rm=rlast, and rm=r1st
				for (j=0;j<2;j++){
					if (j==0)
						r=r1st;
					else
						r=r2nd;
					tmppt=r-rm;	//diff btw pt and prev pt
					tmpdist=sqrt(dot(tmppt,tmppt));
					if (j==0)
						dist += tmpdist;	//close ring: with distance
					
					//compute anisotropy measures
					//tangent and binormal vectors at rm
					localrz(rmm, rm, r, t, b);	//sets tangent, t, and 
												//binormal, b
					//ipar
					tmpdb1=dot(t,v);
					ipar += (1-tmpdb1)*(1+tmpdb1)*tmpdist;
					//iperp1
					tmpdb1=dot(t,vp1);
					iperp1 += (1-tmpdb1)*(1+tmpdb1)*tmpdist;	
					//iperp2
					tmpdb1=dot(t,vp2);
					iperp2 += (1-tmpdb1)*(1+tmpdb1)*tmpdist;	
					//ibinpar
					ibinpar = ibinpar + b*tmpdist;	
					//circulate old pts
					rmm=rm;
					rm=r;
				}
			}
			//ibinpar is actually not parallel
			tmpdb1=sqrt(dot(ibinpar, ibinpar));
			//ibin is parallel
			ibin=dot(ibinpar,v)*v;		
			tmpdb2=sqrt(dot(ibin,ibin));

			ipar /= dist;
			iperp1 /= dist;
			iperp2 /= dist;
			dist *= volinv;
			tmpdb1 *= volinv/(dist*sqrt(dist)); //schwarz has /dist^(3/2)
			tmpdb2 *= volinv/(dist*sqrt(dist)); //schwarz has /dist^(3/2)
	
			cout.precision(15);

			/*	Column	Data
				1 		time
				2		ldens (optional)
				3		ipar
				4		iperp1
				5		iperp2
				6		ibin (forced parallel to v_{ns})
				7		ibinpar (not parallel, necessarily)
				8		ibin/ibinpar
				9		ipar/2+iperp1
				10		ipar/2+iperp2
			*/
			cout << time << " ";
			if (printdist)
				cout << dist << " ";
			cout.precision(7);
			cout << ipar << " " << iperp1;
			//print just component parallel to v
			cout << " " << iperp2 << " " << tmpdb2 << " ";
			//print fraction of total
			cout << tmpdb1 << " " << tmpdb2/tmpdb1 << " ";		

			tmpdb1 = ipar/2 + iperp1;
			cout << tmpdb1 << " ";
			tmpdb1 = ipar/2 + iperp2;
			cout << tmpdb1 << endl;
			
			ifhandle.close();
			num += nskip;
		}
	}
	
	return 0;
}

void localrz(point start, point mid, point end, point& tangent, point& vel){
	
	point lp, lm, curv, tmppt, sumpt, combopt;
    double lpdotp, lmdotm, lpdotm, norm;
    double tmpdb1, tmpdb2, tmpdb;
	double denom, dp, dm;
    
	lp = end - mid;
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
    }else{
        vel.x = 0; vel.y = 0; vel.z = 0;
    }

    return;
}

