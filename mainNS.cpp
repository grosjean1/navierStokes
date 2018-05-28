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

using namespace std;

int  main(int argc, const char** argv)
{
	MatMap M1,M2;
  double nu = 0.0025;
	double dt=0.1;
	double alpha=1./dt;
	vector<double> xprec;
	vector<double> X;
	cout << " lecture de " << argv[1] << endl;
  Mesh2d Th(argv[1]);
	int n=Th.PointsMil();

	X=resolution_Stokes(Th,0,nu,M1,n,xprec,0,0); //RESOLUTION STOKES

	ofstream file("plot/solution.txt");
	int i;

	for(int k=0;k<Th.nbt;k++){
		for(int il=0;il<15;il++){
			if(il<6){
				i=Th(k,il);
				file<<X[i]<<" ";
			}
			else if(il<12){
				i=Th(k,il-6);
				file<<X[i+n]<<" ";
			}
			else{
				i=Th(k,il-12);
				file<<X[i+2*n]<<" ";
			}
		}
		
		file<<"\n";
	}
	file.close();
	xprec=X;

	ofstream file1("plot/sol_0.txt");
	cout<< "pas de temps 0"<<endl;
	M2.clear();
	X=resolution_Stokes(Th,alpha,nu,M2,n,xprec,1,0); //RESOLUTION NAVIER-STOKES
	xprec=X;
	for(int k=0;k<Th.nbt;k++){
		for(int il=0;il<15;il++){
			if(il<6){
				i=Th(k,il);
				file1<<xprec[i]<<" ";
			}
			else if(il<12){
				i=Th(k,il-6);
				file1<<xprec[i+n]<<" ";
			}
			else{
				i=Th(k,il-12);
				file1<<xprec[i+2*n]<<" ";
			}
		}
		file1<<"\n";
	}
	file1.close();

	for(int t=1;t<80;t++){
		string s = "plot/sol_"+to_string(t)+".txt";
		ofstream file2(s.c_str());
		cout<<"pas de temps"<<t<<endl;

		X=resolution_Stokes(Th,alpha,nu,M2,n,xprec,1,1); ////RESOLUTION NAVIER-STOKES EN REUTILISANT LA MAP
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
