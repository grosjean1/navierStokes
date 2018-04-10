#include <vector>
#include <cassert>
#include <fstream>
#include <map>
#include "MatNS.hpp"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ctime>
#include "umfpack.h"
#include <string> 
//#define GNUPLOT_PATH "/usr/bin/gnuplot"


using namespace std;

int  main(int argc, const char** argv)
{
	MatMap M1,M2;
  double nu = 1;
	vector<double> xprec;
	vector<double> X;
	cout << " lecture de " << argv[1] << endl;
  Mesh2d Th(argv[1]);
	int n=Th.PointsMil();
/*	for(int i=0;i<n;i++){
			cout<<Th.v[i].getX()<<" "<<Th.v[i].getY()<<" "<<Th.v[i].getNum()<<" "<<Th.v[i].getLab().lab<<endl;
	}*/

	xprec=resolution_Stokes(Th,nu,M1,n);
	double dt=0.1;
	double alpha=1./dt;
	ofstream file("solution.txt");
	int i;
	
	for(int k=0;k<Th.nbt;k++){
		for(int il=0;il<15;il++){
			if(il<6){
				i=Th(k,il);
				file<<xprec[i]<<" ";
			}
			else if(il<12){
				i=Th(k,il-6);
				file<<xprec[i+n]<<" ";
			}
			else{
				i=Th(k,il-12);
				file<<xprec[i+2*n]<<" ";
			}
		}
		
		file<<"\n";
	}
	for(int t=0;t<5;t++){
cout << t << endl;
		string s = "sol_"+to_string(t)+".txt";
		ofstream file2(s.c_str());
		cout<<"pas de temps"<<t<<endl;
		X=resolution_Navier_Stokes(Th,alpha,nu,M2,n,xprec);
		xprec=X;
		for(int k=0;k<Th.nbt;k++){
			for(int il=0;il<15;il++){
				if(il<6){
					i=Th(k,il);
					file2<<xprec[i]<<" ";
				}
				else if(il<12){
					i=Th(k,il-6);
					file2<<xprec[i+n]<<" ";
				}
				else{
					i=Th(k,il-12);
					file2<<xprec[i+2*n]<<" ";
				}
			}
		
			file2<<"\n";
		}
		file2.close();
	}

	return 0;
}
