#include <cstdlib>
#include "linalg.h"
/*********************************************
 * Tests the function biot() from rkf.c in torus
 * Compares ways to evaluate the integral for numerical stability
 * Does this by rotating testpt about origin, located at 
 * midpoint along source vortex segment
 * Rotates about vectors forming vortex segment and original
 * testpt location (orthogonal & symmetric position wrt vortex seg)
***********************************************
* Input: nothing
* Output:
* 	#ifdef POINTS: {x y z} for endpt, startpt, testpt(difft rotations)
*
* 	#ifdef ORTHOG: prints dot product of basis vectors and
* 					testpt wrt normal vector to plane of rotation
* 	#else:
* 			prints biot() value calculated at difft rotations
* 			of testpt about midpoint of vortex segment
************************************************
*/

//#define POINTS
#ifndef POINTS
	//#define ORTHOG
#endif
point biot(point, point, point);
double TINY1=1e-24, TINY2=1e-12;
int numer=0;	//=0: biot uses std numerator
				//=1: biot uses new numerator
int main(){
	int ptspread=0;		//=1: test point distributed from 0->pi/2
						//=0: distributed in small region near starting pt
	int testinit=1;		//=1: starts colinear and same dir as l2
						//=0; starts colinear and opposite dir as l2
	int num=15000;		//How many evaluation points
	double factr = 1e-11;	//Factor to multiply by if ptspread=0
	double outdb, cntang, testfactr;
	point test, end, st, v[3];
	double scale=RAND_MAX+1.0, max=.01, min=-.01, range;
	double tmpd, dbpt[3], dbpt2[3], va[3][3];	//Basis vecs
	double angle, dang, c, s, Rij;
	point out, tmppt;
	int i, j, k, l, m;

	//Pick random source vortex segment, with
	//midpoint at origin
	range=max-min;
	//srand(time(0));	//Comment to keep same seed
	while(1){
		end.x = range*(rand()/scale)+min;
		end.y = range*(rand()/scale)+min;
		end.z = range*(rand()/scale)+min;
		if ((tmpd=dot(end,end))>1.0e-12){	
			st.x=-end.x;	//origin btwn st and end
			st.y=-end.y;
			st.z=-end.z;
			v[0]=end/sqrt(tmpd);	//basis vec 1
			break;
		}
	}
	#ifdef POINTS
	cout << end << endl;
	cout << st << endl;
	#endif
	//Build other basis vectors for rotation of testpoint
	v[1].x=1.0; v[1].y=0.0; v[1].z=0.0;
	for (i=0;i<2;i++){
		tmppt=v[1]-dot(v[0],v[1])*v[0];
		if (dot(tmppt,tmppt)<1e-8){
			v[1].x=0.0;
			v[1].y=1.0;
		}else{
			tmppt=tmppt/sqrt(dot(tmppt,tmppt));
			v[1]=tmppt;
			break;
		}
	}
	v[2]=cross(v[0],v[1]);
	for (m=0;m<2;m++){
		if (m){
			TINY1=1e-14;
			TINY2=1e-12;
		}
		testfactr = 0.05;
		for (l=0;l<10;l++){
			if (l)
				testfactr += 0.1;
			cntang=0.0;
			//INITIAL TEST POINT LOCATION
			if (testinit==0){
				//test starts colinear (and opposite) of l2=vortex seg
				test = -testfactr*max*v[0];
			}else{
				//test starts colinear and same dir as l2=vortex seg
				test = testfactr*max*v[0];
			}
		
			#ifdef POINTS
			cout << test << endl;
			#elif defined(ORTHOG)
			cout << dot(test,v[2]) << endl;
			#else
			out=biot(st, end, test);
			outdb=sqrt(dot(out,out));
			cout << 0.0 << " " << outdb << endl;
			#endif
			//To iterate Rij, need to make point structures into arrays:
			//v[i].x,y,z = va[i][0,1,2]
			for (i=0;i<3;i++){
				va[i][0]=v[i].x;	//Notation: va = varray
				va[i][1]=v[i].y;
				va[i][2]=v[i].z;
			}
			if (ptspread)
				dang=M_PI/num;
			else
				dang=M_PI*factr;
			if (testinit==0){
				angle = -dang;	//Starts colinear and opposite of l2
			}else{ 
				angle = dang;	//Starts colinear and same dir as l2
			}
			//Rotate testpt about center of vortex segment (=origin)
			for (k=0;k<num;k++){
				//Rounding error might be creeping in
				c=cos(angle);
				s=sin(angle);
				dbpt[0]=test.x;
				dbpt[1]=test.y;
				dbpt[2]=test.z;
				for(i=0;i<3;i++){
					dbpt2[i]=0.0;
					for(j=0;j<3;j++){
						Rij = (c*va[0][i]+s*va[1][i])*va[0][j] + 
								(-s*va[0][i]+c*va[1][i])*va[1][j] +
								va[2][i]*va[2][j];	//3d rotation matrix
													//abt v[0],v[1] plane
						dbpt2[i] += Rij*dbpt[j];
						
					}
				}
				test.x=dbpt2[0];
				test.y=dbpt2[1];
				test.z=dbpt2[2];
			#ifdef POINTS	
				cout << test << endl;
			#elif defined(ORTHOG)
				cout << dot(test,v[2]) << endl;
			#else
				cntang+=dang;
				out=biot(st, end, test);
				outdb=sqrt(dot(out,out));
				cout << cntang << " " << outdb << endl;
			#endif
				if (testinit==0)
					angle = -dang;	//Starts colinear and opposite
				else
					angle = +dang;	//Starts colinear and same dir
			}
		}
	}
}

point biot(point stseg, point endseg, point testpt)
{
	point l1, l2, l3, deldot;
	double dot11, dot12, dot22, dot33, dot23, factor1, factor2;
	double tmp1, tmp2, tmp3, d1d2, dot13, d1d3;
	point dir;

	l1 = testpt - stseg;
	l2 = endseg - stseg;
	l3 = testpt - endseg;      
      
	dot11 = dot(l1, l1);
	dot12 = dot(l1, l2);
	dot22 = dot(l2, l2);
	dot23 = dot12 - dot22;
	dot33 = dot(l3, l3);
	
	factor1 = dot11*dot22 - dot12*dot12;
		
	if ((tmp1=fabs(factor1/dot11/dot22)) > TINY1){
		factor2 = dot23/sqrt(dot33) - dot12/sqrt(dot11);
		deldot = (factor2 / factor1) * cross(l1,l2);
	}else{
		factor1 = dot22*dot33-dot23*dot23;
		if  ((tmp2=fabs(factor1/dot22/dot33)) > TINY1){
			factor2 = dot23/sqrt(dot33) - dot12/sqrt(dot11);
			deldot = (factor2 / factor1) * cross(l3, l2);
		}else{
			if ((tmp3=dot11*dot33) > TINY2 && dot12*dot23 > 0.0) {
				factor2 = fabs(dot11-dot33)/dot11/dot33/sqrt(dot22);
				deldot = (factor2 / 2) * cross(l1,l2);
			}else if (dot12*dot23<=0.0){	//test point near vortex seg
				dot13=dot(l1,l3);
				d1d3=dot11*dot33;
				if (!numer){
					factor2 = dot23/sqrt(dot33) - dot12/sqrt(dot11);
				}else{
					factor2=-2.0*sqrt(dot22/d1d3)*dot13;
				}
				dir = cross(l1,l2);
				factor1 = dot(dir, dir);
				deldot = (factor2 / factor1) * dir;
			}else{
				deldot.x = 0; deldot.y = 0; deldot.z = 0;
			}
		}
	}
	return deldot;
}
