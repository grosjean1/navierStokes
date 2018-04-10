#include <cassert>
#include "Fonctions_Utiles.hpp"
#include <fstream>
#include <iostream>
#include <math.h> 
# include <iomanip>
# include <ctime>
# include "umfpack.h"

using namespace std;


typedef map< pair<int,int>,double> MatMap;

# define TIME_SIZE 40
void timestamp(){
  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;
  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );
  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );
  std::cout << time_buffer << "\n";
}


//////////////////////////////////////// Equation de Stokes stationnaire /////////////////////////
vector<double> resolution_Stokes(Mesh2d Th,double nu, MatMap & M,int n){
	int nt=Th.nbt;
	//ofstream StokesMatElement("MaMat.txt");
	vector<double> solution;
	double b[2*n+Th.nv]; //2nd membre
	for(int i=0;i<2*n+Th.nv;i++){
		b[i]=0; //second membre
	}
	for(int k=0;k<nt;k++){
		double A[15][15];
		BuildMatNS(Th, 0,nu, A,k,0);
		int i,j;
		for(int il=0;il<15;il++){//mat glob assemblage: Th(k,il) renvoie le num glob d'un sommet si il<3, et sinon d'un milieu si il<6
			for(int jl=0;jl<15;jl++){
				if(il<6){
					i=Th(k,il);
				}
				else if(il<12){
					i=Th(k,il-6)+n;
				}
				else{
					i=Th(k,il-12)+2*n;
				}
				if(jl<6){
					j=Th(k,jl);
				}
				else if(jl<12){
					j=Th(k,jl-6)+n;
				}
				else{
					j=Th(k,jl-12)+2*n;
				}
				if(fabs(A[il][jl])>1e-15){
					M[make_pair(i,j)]+=A[il][jl];
				}
			}
		}
	}

	cout<<endl;
	/*for (std::map< pair<int,int>,double>::iterator it=M.begin(); it!=M.end(); ++it)
  	 std::cout << it->first.first<<" "<< it->first.second<<" "<<it->second << endl;*/
	int lab[6];//Condition aux limites
	for(int k=0;k<nt;k++){
		for(int il=0;il<6;il++){
			if(il<3)
				lab[il]=Th.t[k].v[il].getLab().OnGamma();
			else
				lab[il]=Th.t[k].mil[il-3].getLab().OnGamma();
			if(lab[il]==10||lab[il]==20||lab[il]==40){//BORDS (30 = sortie) 
				int i1,i2,i3;
				i1=Th(k,il);
				i2=Th(k,il)+n;
				//i3=Th(k,il)+2*n;
				pair <int,int> key1=make_pair(i1,i1);
				pair <int,int> key2=make_pair(i2,i2);
				//pair <int,int> key3=make_pair(i3,i3);

				if(M.find(key1)!= M.end())//si on trouve dans la map
					M[key1]=tgv;
				if(M.find(key2)!= M.end())//si on trouve dans la map
					M[key2]=tgv;
				if(il<3){
					b[i1]=g(Th.t[k].v[il],lab[il])*tgv;  ///Conditions aux bords
					b[i2]=0;
				}
				else{
					b[i1]=g(Th.t[k].mil[il-3],lab[il])*tgv;
					b[i2]=0;
				}	
			}
		}
	}

/*	ofstream file2("mat.txt");
  for ( int i = 0; i < 2*n+Th.nv; i++ ) {
 		for ( int j = 0; j < 2*n+Th.nv; j++ ){
			pair<int, int> key=make_pair(i,j);
			if(M.find(key)!= M.end()){//si on trouve dans la map
    		file2 << i<<" "<<j <<" "<< M[key]<<endl; 
  		}
		}
	}*/
	//SparseMatrix
	int VNN=M.size();
	int AI[VNN]; //indices des lignes
	double Ax[VNN]; //valeurs
	int Ap[2*n+Th.nv+1]; //nombre d'elmts non nuls par col
	Ap[0]=0;
	int cpt=0;
	for(int j=0;j<2*n+Th.nv;j++){
		int cpt2=0;//pour Ap
		for (std::map< pair<int,int>, double>::iterator it=M.begin(); it!=M.end(); ++it){
			if(it->first.second==j){
				pair<int,int> key;
				key.first=it->first.first;
				key.second=j;
				AI[cpt]=key.first;
				Ax[cpt]=M[it->first];
				cpt++;
				cpt2++;
			}
		}	
		Ap[j+1]=Ap[j]+cpt2;
	}
	int taille = 2*n+Th.nv;
/*	for(int i=0;i<taille;i++){
		cout<<"b["<<i<<"]="<<b[i]<<endl;
	}*/
  double *null = ( double * ) NULL;
  void *Numeric;
  int status;
  void *Symbolic;
  double x[taille];

  timestamp ( );

  status = umfpack_di_symbolic ( taille, taille, Ap, AI, Ax, &Symbolic, null, null );
  status = umfpack_di_numeric (Ap, AI, Ax, Symbolic, &Numeric, null, null );
  umfpack_di_free_symbolic ( &Symbolic );
//  Solve the linear system.
  status = umfpack_di_solve ( UMFPACK_A, Ap, AI, Ax, x, b, Numeric, null, null );
  umfpack_di_free_numeric ( &Numeric );
  cout << "\n";
  cout << "  Computed solution Stokes\n";
  cout << "\n";
	
	//ofstream file("solutionS.txt");
  for ( int i = 0; i < taille; i++ ) 
  {
  //  file << "  x[" << i << "] = " << x[i] << "\n";
		solution.push_back(x[i]);
  }

//  Terminate.
//
  cout << "\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );
	return solution;
}



//VNN = val non nuls de M; AI et Ax (de taille VNN), n degre de liberte de P2
vector<double> resolution_Navier_Stokes(Mesh2d Th,double alpha,double nu, MatMap & M,int n,vector<double> xprec){
	int nt=Th.nbt;
	double b[2*n+Th.nv]; //2nd membre
	for(int i=0;i<2*n+Th.nv;i++){
		b[i]=0; //init second membre
	}
	for(int k=0;k<nt;k++){
		double A[15][15];
		BuildMatNS(Th, alpha,nu, A,k,1);
		int i,j;
		for(int il=0;il<15;il++){//mat glob assemblage: Th(k,il) renvoie le num glob d'un sommet si il<3, et sinon d'un milieu si il<6
			for(int jl=0;jl<15;jl++){
				if(il<6){
					i=Th(k,il);
				}
				else if(il<12){
					i=Th(k,il-6)+n;
				}
				else{
					i=Th(k,il-12)+2*n;
				}
				if(jl<6){
					j=Th(k,jl);
				}
				else if(jl<12){
					j=Th(k,jl-6)+n;
				}
				else{
					j=Th(k,jl-12)+2*n;
				}
				if(fabs(A[il][jl])>1e-15){
					M[make_pair(i,j)]+=A[il][jl];
				}
			}
		}
	}
	for(int i=2*n;i<2*n+Th.nv;i++){
		M[make_pair(i,i)]=-(10e-8);
	}
	//cout<< "taille de la map "<<M.size()<<endl;
	double u1pks[3]={0};
	double u2pks[3]={0};
	double xp_as[3]={0};
	double yp_as[3]={0};
	double phi[3]={0};
	vector<R2> p;
	double u1pNvks[3]={0}; double u2pNvks[3]={0};
	int vois;
	int i;
	double c;
	vector<double> u1pk,u2pk;
	vector<double> u1pNvk; vector<double> u2pNvk;
	for(int k=0; k<nt;k++){
		u1pk.clear();u2pk.clear();u1pNvk.clear();u2pNvk.clear();p.clear();
		double areak=Th.t[k].area;
		c=0;
		calculPtInterp(Th.t[k],p);
		for(int il=0;il<15;il++){
			if(il<6){
				i=Th(k,il);
				phi[0]=Phi(il,p[0]);
				phi[1]=Phi(il,p[1]);
				phi[2]=Phi(il,p[2]);
			}
			else if(il<12){
				i=Th(k,il-6)+n;
				phi[0]=Phi(il-6,p[0]);
				phi[1]=Phi(il-6,p[1]);
				phi[2]=Phi(il-6,p[2]);
			}
			else{
				i=Th(k,il-12)+2*n;
			}
			recup(Th.t[k],xprec,u1pk,u2pk,n);//calcul de u1pk et u2pk
			
			for(int j=0;j<3;j++){

				u1pks[j]=vitesseInterpolee(u1pk,p[j]);
				u2pks[j]=vitesseInterpolee(u2pk,p[j]);
				//cout<<"vitesse interpolee"<<u1pks[j]<<" "<<u2pks[j]<<endl;
				xp_as[j]= p[j].x-(1./alpha)*u1pks[j]; //Position du point d'interpolation au pas précédent
				yp_as[j]=p[j].y-(1./alpha)*u2pks[j];
				R2 nvPt(xp_as[j],yp_as[j]);
				//cout<<"k"<<k<<" "<<p[j].x<<p[j].y<<endl;
				//cout<<xp_as[j]<<" "<<yp_as[j]<<endl;
				vois=RecupVoisins(Th,k,nvPt);
				//cout<<"voisin "<<vois<<endl;
				recup(Th.t[vois],xprec,u1pNvk,u2pNvk,n);//calcul de u1pk et u2pk
				u1pNvks[j]=vitesseInterpolee(u1pNvk,nvPt);
				u2pNvks[j]=vitesseInterpolee(u2pNvk,nvPt);
			}
			if(i<n){
				for(int j=0;j<3;j++){
					c+=phi[j]*u1pNvks[j];
				}
				b[i]+=alpha*areak*(1./3)*c;
			}
			else if(i<2*n){
				for(int j=0;j<3;j++){
					c+=phi[j]*u2pNvks[j];
				}
				b[i]+=alpha*areak*(1./3)*c;
			}
			else
				b[i]=0; //rien sur la pression
		}
	}
	/*for (std::map< pair<int,int>,double>::iterator it=M.begin(); it!=M.end(); ++it)
  	 std::cout << it->first.first<<" "<< it->first.second<<" "<<it->second << endl;*/
	int lab[6];//Condition aux limites
	for(int k=0;k<nt;k++){
		for(int il=0;il<6;il++){
			if(il<3)
				lab[il]=Th.t[k].v[il].getLab().OnGamma();
			else
				lab[il]=Th.t[k].mil[il-3].getLab().OnGamma();
			if(lab[il]==10||lab[il]==20||lab[il]==40){//BORDS (30 = sortie) 
				int i1,i2;
				i1=Th(k,il);
				i2=Th(k,il)+n;
				pair <int,int> key1=make_pair(i1,i1);
				pair <int,int> key2=make_pair(i2,i2);
				//cout<<"i"<<i<<endl;
				if(M.find(key1)!= M.end())//si on trouve dans la map
					M[key1]=tgv;
				if(M.find(key2)!= M.end())//si on trouve dans la map
					M[key2]=tgv;
				if(il<3){
					b[i1]=g(Th.t[k].v[il],lab[il])*tgv;  ///Conditions aux bords
					b[i2]=0;
				}
				else{
					b[i1]=g(Th.t[k].mil[il-3],lab[il])*tgv;
					b[i2]=0;
				}	
			}
		}
	}
	/*ofstream file3("mat2.txt");
  for ( int i = 0; i < 2*n+Th.nv; i++ ) {
 		for ( int j = 0; j < 2*n+Th.nv; j++ ){
			pair<int, int> key=make_pair(i,j);
			if(M.find(key)!= M.end()){//si on trouve dans la map
    		file3 << i<<" "<<j <<" "<< M[key]<<endl; 
  		}
		}
	}*/

	//SparseMatrix
	int VNN=M.size();
	int AI[VNN]; //indices des lignes
	double Ax[VNN]; //valeurs
	int Ap[2*n+Th.nv+1]; //nombre d'elmts non nuls par col
	Ap[0]=0;
	int cpt=0;
	for(int j=0;j<2*n+Th.nv;j++){
		int cpt2=0;//pour Ap
		for (std::map< pair<int,int>, double>::iterator it=M.begin(); it!=M.end(); ++it){
			if(it->first.second==j){
				pair<int,int> key;
				key.first=it->first.first;
				key.second=j;
				AI[cpt]=key.first;
				Ax[cpt]=M[it->first];
				cpt++;
				cpt2++;
			}
		}	
		Ap[j+1]=Ap[j]+cpt2;
	}
  int taille = 2*n+Th.nv;
	//SparseMatrix
	//cout<< "i j" << taille << " Nombre de valeurs non nulles " << VNN <<endl; 
  double *null = ( double * ) NULL;
  void *Numeric;
  int status;
  void *Symbolic;
  double x[taille];

  timestamp ( );


  status = umfpack_di_symbolic ( taille, taille, Ap, AI, Ax, &Symbolic, null, null );
  status = umfpack_di_numeric (Ap, AI, Ax, Symbolic, &Numeric, null, null );
  umfpack_di_free_symbolic ( &Symbolic );
//  Solve the linear system.
  status = umfpack_di_solve ( UMFPACK_A, Ap, AI, Ax, x, b, Numeric, null, null );
  umfpack_di_free_numeric ( &Numeric );
  cout << "\n";
  cout << "  Computed solution NS:\n";
  cout << "\n";
	//ofstream file("solutionNS.txt");
	vector<double> solution;
  for ( int i = 0; i < taille; i++ ) 
  {
  //  file << "  x[" << i << "] = " << x[i] << "\n";
		solution.push_back(x[i]);
  }
	//file << "\n";
//  Terminate.
//
  cout << "\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );
	return solution;
}

