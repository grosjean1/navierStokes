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
vector<double> resolution_Stokes(Mesh2d & Th,double alpha,double nu, MatMap & M,int n,vector<double> xprec, int NS,bool MapExiste){
	int nt=Th.nbt;
	//ofstream StokesMatElement("MaMat.txt");
	vector<double> solution;
	double b[2*n+Th.nv]; //2nd membre
	for(int i=0;i<2*n+Th.nv;i++){
		b[i]=0; //second membre
	}
	if(MapExiste==0){
		for(int k=0;k<nt;k++){
			double A[15][15];
			BuildMatNS(Th, alpha,nu, A,k);
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
	}
	if(NS==1){
		//cout<<"calcul caract"<<endl;
		CalculCaracteristique(Th,alpha,xprec,n,b);
	}	
	//cout<<"fin carac "<<endl;
	
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
				if(MapExiste==0){
					pair <int,int> key1=make_pair(i1,i1);
					pair <int,int> key2=make_pair(i2,i2);
					if(M.find(key1)!= M.end())//si on trouve dans la map
						M[key1]=tgv;
					if(M.find(key2)!= M.end())//si on trouve dans la map
						M[key2]=tgv;
				}
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

//SparseMatrix
	//cout << " build sparse mat " << endl;
	int VNN=M.size();
	int * AI=new int[VNN]; //indices des lignes
	double * Ax=new double[VNN]; //valeurs
	int * Ap=new int[2*n+Th.nv+1]; //nombre d'elmts non nuls par col
	Ap[0]=0;
	int cpt=0;
	for (std::map< pair<int,int>, double>::iterator it=M.begin(); it!=M.end(); ++it)
  {
    int j = it->first.second;
		int i= it->first.first;
		AI[cpt]=j;
		Ax[cpt]=it->second;
    Ap[i+1] = ++ cpt;
	}
	//cout << " FAC  sparse mat " << endl;
	int taille = 2*n+Th.nv;

  double *null = ( double * ) NULL;
  void *Numeric;
  int status;
  void *Symbolic;
  double x[taille];

  timestamp ( );
	//besoin seulement de solve si on passe en copie Ap AI et Ax
  status = umfpack_di_symbolic ( taille, taille, Ap, AI, Ax, &Symbolic, null, null );
  status = umfpack_di_numeric (Ap, AI, Ax, Symbolic, &Numeric, null, null );
	//cout << " SOLV  sparse mat " << endl;
  umfpack_di_free_symbolic ( &Symbolic );
	//  Solve the linear system.
  status = umfpack_di_solve ( UMFPACK_At, Ap, AI, Ax, x, b, Numeric, null, null );
  umfpack_di_free_numeric ( &Numeric );
  cout << "\n";
	if(NS==0){
 	 cout << "  Computed solution Stokes\n";
	}
	else{
		cout << "  Computed solution Navier - Stokes\n";
	}
  cout << "\n";
	
  for ( int i = 0; i < taille; i++ ) 
  {
		solution.push_back(x[i]);
  }

  cout << "\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );
	return solution;
}


