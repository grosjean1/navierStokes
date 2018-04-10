#include <cassert>
#include "mesh.hpp"
#include <fstream>
#include <iostream>
#include <math.h> 
#include <algorithm>
double tgv=10e30;
//////////////////////////////////////////////////////        FONCTIONS DE BASE      /////////////////////////////////
double lambda(int i,R2 P){
	double result;
	if(i==0)
			result=1-P.x-P.y;
	else if(i==1)
			result=P.x;
	else
			result=P.y;
	return result;
}

int PartialLamb(int i,int ind){ //ind=0 pour x; 1 pour y
	if(i==0)
		return (-1);
	else if((i==1 && ind==0)||(i==2 && ind==1))
		return 1;
	else 
		return 0;
}

double Phi(int i, R2 P){
	if(i<3)
		return (lambda(i,P)*(2*lambda(i,P)-1));
	else
		return (4*lambda((i-1)%3,P)*lambda((i-2)%3,P));
}

double PartialPhi(int i,int ind,R2 P){
	if(i<3)
		return ((4*lambda(i,P)-1)*PartialLamb(i,ind));
	else
		return (4*(lambda((i+1)%3,P)*PartialLamb((i+2)%3,ind)+lambda((i+2)%3,P)*PartialLamb((i+1)%3,ind)));
}


//////////////////////////////////////////////////////        MATRICE ELEMENTAIRE   NAVIER-STOKES  (NS=0 pour pb stationnaire) /////////////////////////////////
void BuildMatNS(const Mesh2d &Th, double alpha, double nu, double A[15][15],int cpt,int NS) //A mat de taille [15][15]
{
	//pts et poids d'integration
	double ptint1=(6-sqrt(15))/21;
	double ptint2=(9-sqrt(15)*2)/21;
	double ptint3=(6+sqrt(15))/21;
	double ptint4=(9+sqrt(15)*2)/21;
	double poids1=(155-sqrt(15))/1200;
	double poids2=(155+sqrt(15))/1200;
	R2 p1(1./3,1./3),p2(ptint1,ptint1),p3(ptint1,ptint4),p4(ptint4,ptint1),p5(ptint3,ptint3),p6(ptint3,ptint2),p7(ptint2,ptint3);
	R2 PtsRef[7]={p1,p2,p3,p4,p5,p6,p7};//Points d'intégration

	double Poids[7]={0.225,poids1,poids1,poids1,poids2,poids2,poids2};
//	R2 p1(0.5,0.5),p2(0,0.5),p3(0.5,0);
//	R2 PtsRef[3]={p1,p2,p3};
//	double Poids[3]={1./3,1./3,1./3};
  Triangle K = Th[cpt]; // triangle cpt du maillage
  // calcul de la matrice Ak
  double areak = K.area;
	double coeff=nu/(4*areak);
	double coeff1=alpha*areak;
  double J[2][2];//transformation affine
	J[0][0]=K.v[2].getY()-K.v[0].getY();
	J[0][1]=K.v[0].getY()-K.v[1].getY();
	J[1][0]=K.v[0].getX()-K.v[2].getX();
	J[1][1]=K.v[1].getX()-K.v[0].getX();

	for(int i=0;i<15;i++){//init de la mat élémentaire (C  0  B1 
		for(int j=0;j<15;j++){											//    0  C  B2
			A[i][j]=0;                               //     B1 B2 -Ieps)
		}
	}
	double acoef=pow(J[0][0],2)+pow(J[1][0],2);   //Bk'*Bk
	double bcoef=J[0][0]*J[0][1]+J[1][0]*J[1][1];
	double ccoef=pow(J[0][1],2)+pow(J[1][1],2);
	if(NS==1){ //NS
		for(int i=0;i<6;i++){ //C
			for(int j=0;j<6;j++){
				for(int k=0;k<7;k++){//les pts d'integ
					A[i][j]+=coeff1*Poids[k]*Phi(j,PtsRef[k])*Phi(i,PtsRef[k])+coeff*Poids[k]*(acoef*PartialPhi(i,0,PtsRef[k])*PartialPhi(j,0,PtsRef[k])+bcoef*(PartialPhi(i,1,PtsRef[k])*PartialPhi(j,0,PtsRef[k])+PartialPhi(i,0,PtsRef[k])*PartialPhi(j,1,PtsRef[k]))+ccoef*PartialPhi(i,1,PtsRef[k])*PartialPhi(j,1,PtsRef[k]));
				}
				A[i+6][j+6]=A[i][j];
			}
		}
	}
	else{ //Cas stationnaire
		for(int i=0;i<6;i++){ //C
			for(int j=0;j<6;j++){
				for(int k=0;k<7;k++){//les pts d'integ
					A[i][j]+=coeff*Poids[k]*(acoef*PartialPhi(i,0,PtsRef[k])*PartialPhi(j,0,PtsRef[k])+bcoef*(PartialPhi(i,1,PtsRef[k])*PartialPhi(j,0,PtsRef[k])+PartialPhi(i,0,PtsRef[k])*PartialPhi(j,1,PtsRef[k]))+ccoef*PartialPhi(i,1,PtsRef[k])*PartialPhi(j,1,PtsRef[k]));
								
				}
				A[i+6][j+6]=A[i][j];
			}
		}
	}
	for(int i=0;i<6;i++){
		for(int j=12;j<15;j++){ //B1
			for(int k=0;k<7;k++){
				A[i][j]+=Poids[k]*((J[0][0]*PartialPhi(i,0,PtsRef[k])+J[0][1]*PartialPhi(i,1,PtsRef[k]))*lambda(j-12,PtsRef[k]));
			}
			A[i][j]=-(1./2)*A[i][j];
			A[j][i]=A[i][j]; //matrice symetrique
		}
	}
	for(int i=6;i<12;i++){
		for(int j=12;j<15;j++){ //B2
			for(int k=0;k<7;k++){
				A[i][j]+=Poids[k]*(J[1][0]*PartialPhi(i-6,0,PtsRef[k])+J[1][1]*PartialPhi(i-6,1,PtsRef[k]))*lambda(j-12,PtsRef[k]);
			}
			A[i][j]=-(1./2)*A[i][j];
			A[j][i]=A[i][j];////matrice symetrique
		}
	}
	for(int i=12;i<15;i++){//-eps
		A[i][i]=-(10e-8);
	}
}

/////////////////////////////////////////////  FONCTIONS UTILES   ///////////////////////////////////////////////
//Fonction qui permet de récupérer les valeurs de u^n aux points du triangle t
void recup(Triangle t, vector<double> un, vector<double>& u1n, vector<double>& u2n, int n){ //n = degre de liberte P2
	int tab[6]={0};
	u1n.clear();
	u2n.clear();

	for(int il=0;il<6;il++){
			if(il<3)
				tab[il]=t.v[il].getNum();
			else
				tab[il]=t.mil[il-3].getNum();
			//cout<<"num glob"<<tab[il]<<endl;
	}
	for(int i=0;i<6;i++){
		u1n.push_back(un[tab[i]]);
		u2n.push_back(un[tab[i]+n]);
	}
}

double vitesseInterpolee(vector<double> un, R2 PtInterp){
	double uninterp=0;
	for(int i=0;i<6;i++){
		uninterp=uninterp+Phi(i, PtInterp)*un[i];
	}
	return uninterp;
}


// fonction CL
double g(Vertex P, int label)
{
	double x=P.getX();
	double y=P.getY();
  if(label == 10)
      return (1-y)*(y-0.5)*16;
  else 
		return 0;
}

//Fonction qui renvoie le triangle auquel le point PtInterp appartient  
int RecupVoisins(Mesh2d Th, int triangle, R2 PtInterp){
	//cout<<"MIN"<<min(1,2)<<endl;
	double result;double area0;double area2,area3,area4;
	for(int i=0;i<(int)Th.voisins[triangle].size();i++){
		int j=Th.voisins[triangle][i];
		Vertex v0=Th.t[j].v[0];
		Vertex v1=Th.t[j].v[1];
		Vertex v2=Th.t[j].v[2];
		Vertex P;
		P.setX(PtInterp.x);P.setY(PtInterp.y);
		double a=v0.getm(v1);
		double b=v1.getm(v2);
		double c=v2.getm(v0);
		double d=v0.getm(P);
		double e=v1.getm(P);
		double f=v2.getm(P);
		double p0=(a+b+c)/2; //semi-perimetre
		double p1=(a+d+e)/2;
		double p2=(b+e+f)/2; 
		double p3=(c+d+f)/2; 
		area0=sqrt(p0*(p0-a)*(p0-b)*(p0-c));//aire du triangle
		area2=sqrt(p1*(p1-a)*(p1-d)*(p1-e));
		area3=sqrt(p2*(p2-b)*(p2-e)*(p2-f));
		area4=sqrt(p3*(p3-c)*(p3-d)*(p3-f));
		result=(area2+area3+area4);
		//cout<<"aire"<<j<<" "<<"resultat"<<result<<" "<<area0<<endl;
		if(fabs(result-area0)<10e-16){
			//cout<<j<<endl;
			return j;
		}
	}
	return triangle;
}

//Fonction qui retourne les points de quadrature dans le triangle t
void calculPtInterp(Triangle t,vector<R2> &points){
	double v0x=t.v[0].getX();double v0y=t.v[0].getY();
	double v1x=t.v[1].getX();double v1y=t.v[1].getY();
	double v2x=t.v[2].getX();double v2y=t.v[2].getY();

	R2 p1(0.5*v1x+0.5*v2x,0.5*v1y+0.5*v2y);
	R2 p2(0.5*v0x+0.5*v1x,0.5*v0y+0.5*v1y);
	R2 p3(0.5*v0x+0.5*v2x,0.5*v0y+0.5*v2y);
	points.push_back(p1);points.push_back(p2);points.push_back(p3);
}

