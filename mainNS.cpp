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
	xprec=resolution_Stokes(Th,nu,M1,n);
//	xprec=X;

	double dt=0.1;
	double alpha=1./dt;
	//fostream file("solution.dat");
	//for(int i=0;i
	//cout<<"youlou"<<M.size()<<endl;
	for(int i=0;i<5;i++){
		cout<<"pas de temps"<<i<<endl;
		X=resolution_Navier_Stokes(Th,alpha,nu,M2,n,xprec);
		xprec=X;
		for(int i=0;i<X.size();i++){
			cout<<"X["<<i<<"]="<<X[i]<<endl;
		}
	}
	/*	if(i==4){
			for(int i=0;i<X.size();i++){
				cout<<"x["<<i<<"] = "<<X[i
	}*/
	
	return 0;
}

