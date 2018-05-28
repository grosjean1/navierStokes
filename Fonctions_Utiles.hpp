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
void BuildMatNS(const Mesh2d &Th, double alpha, double nu, double A[15][15],int cpt) //A mat de taille [15][15]
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
  Triangle K = Th[cpt]; // triangle cpt du maillage
  // calcul de la matrice Ak
  double areak = K.area;
	double coeff=nu/(4*areak);
	double coeff1=alpha*areak;
	assert(coeff1>=0);
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

	for(int i=0;i<6;i++){ //C
		for(int j=0;j<6;j++){
			for(int k=0;k<7;k++){//les pts d'integ
				A[i][j]+=coeff1*Poids[k]*Phi(j,PtsRef[k])*Phi(i,PtsRef[k])+coeff*Poids[k]*(acoef*PartialPhi(i,0,PtsRef[k])*PartialPhi(j,0,PtsRef[k])+bcoef*(PartialPhi(i,1,PtsRef[k])*PartialPhi(j,0,PtsRef[k])+PartialPhi(i,0,PtsRef[k])*PartialPhi(j,1,PtsRef[k]))+ccoef*PartialPhi(i,1,PtsRef[k])*PartialPhi(j,1,PtsRef[k]));
							
			}
			A[i+6][j+6]=A[i][j];
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
//Fonction qui permet de récupérer les valeurs de u1^n et u2^n aux sommets et milieu du triangle t
void recup(Triangle t, vector<double> un, double * u1n, double * u2n, int n){ //n = degre de liberte P2
	int * tab=new int[6];
	for(int il=0;il<6;il++){
			if(il<3)
				tab[il]=t.v[il].getNum();
			else
				tab[il]=t.mil[il-3].getNum();
	}

	for(int i=0;i<6;i++){
		u1n[i]=un[tab[i]];
		u2n[i]=un[tab[i]+n];
	}
	delete[] tab;
}


//calcule sum (u_i phi_i)
double vitesseInterpolee(double * un, R2 PtInterp){
	double uninterp=0;
	for(int i=0;i<6;i++){
		uninterp=uninterp+Phi(i,PtInterp)*un[i];
	}
	return uninterp;
}


// fonction CL
double g(Vertex P, int label)
{
	//double x=P.getX();
	double y=P.getY();
  if(label == 10)
      return (1-y)*(y-0.5)*16; 
	else 
		return 0; 
	//if(label == 10 ||label==20 ||label==40)
		//return (4*y*(1-y));
}


double min(double x,double y,double z){
	return min(min(x,y),z);
}

//Recupere le triangle voisin auquel appartient le point PtInterp et le projette dans le triangle de ref (PtNv)
int RecupVoisins(Mesh2d & Th, int triangle, R2 PtInterp,R2 & PtNv){
	double area0,area1,area2;
	Vertex v0t=Th.t[triangle].v[0];
	Vertex v1t=Th.t[triangle].v[1];
	Vertex v2t=Th.t[triangle].v[2];

	double v0x=v0t.getX();
	double v0y=v0t.getY();
	double v1x=v1t.getX();
	double v1y=v1t.getY();
	double v2x=v2t.getX();
	double v2y=v2t.getY();

	area0=det(PtInterp,v1t,v2t)*0.5;
	area1=det(v0t,PtInterp,v2t)*0.5;
	area2=det(v0t,v1t,PtInterp)*0.5;

	//double aire=area0+area1+area2;
	//cout<<"a0 "<<area0<<" a1 "<<area1<< " a2 "<<area2<< " aire "<<aire<<endl;
	double rm=min(area0,area1,area2);
	double d=(v1x-v0x)*(v2y-v0y)-(v1y-v0y)*(v2x-v0x);
	double x=PtInterp.x; double y=PtInterp.y;
	PtNv.x=(1/d)*((v2y-v0y)*(x-v0x)+(v0x-v2x)*(y-v0y));
	PtNv.y=(1/d)*((v0y-v1y)*(x-v0x)+(v1x-v0x)*(y-v0y));
	if(rm>=0){
		return triangle;
	}

	for(int i=0;i<(int)Th.voisins[triangle].size();i++){
		int j=Th.voisins[triangle][i];
		Vertex v0=Th.t[j].v[0];
		Vertex v1=Th.t[j].v[1];
		Vertex v2=Th.t[j].v[2];
		double v0x=v0.getX();
		double v0y=v0.getY();
		double v1x=v1.getX();
		double v1y=v1.getY();
		double v2x=v2.getX();
		double v2y=v2.getY();
		area0=det(PtInterp,v1,v2)*0.5;
		area1=det(v0,PtInterp,v2)*0.5;
		area2=det(v0,v1,PtInterp)*0.5;

		double rm=min(area0,area1,area2);
		double d=(v1x-v0x)*(v2y-v0y)-(v1y-v0y)*(v2x-v0x);
		double x=PtInterp.x; double y=PtInterp.y;
		PtNv.x=(1/d)*((v2y-v0y)*(x-v0x)+(v0x-v2x)*(y-v0y));
		PtNv.y=(1/d)*((v0y-v1y)*(x-v0x)+(v1x-v0x)*(y-v0y));
		if(rm>=0){
			return j;
		}
	}
	//cout<<"erreur"<<endl;
	PtNv.x=PtNv.y=sqrt(-1.);
	return (-1);
}

int find_triangle(R2 nvPt, Mesh2d & Th){
	for(unsigned int i=0;i<Th.triangleSortie.size();i++){
		int j=Th.triangleSortie[i];
		Vertex v0=Th.t[j].v[0];
		Vertex v1=Th.t[j].v[1];
		Vertex v2=Th.t[j].v[2];
		Vertex p;p.setX(nvPt.x);p.setY(nvPt.y);
		double area0=det(p,v1,v2)*0.5;
		double area1=det(v0,p,v2)*0.5;
		double area2=det(v0,v1,p)*0.5;
		double rm=min(area0,area1,area2);
		if(rm>=0){
			return j;
		}
	}
	cout<<"erreur: pas de triangle trouvé pour le point "<< nvPt<<"!"<<endl;
	return (-1);
}

//Fonction qui retourne les points de quadrature dans le triangle t
void PointK(Triangle t,R2 * PtsRef, R2 * points){
	double v0x=t.v[0].getX();double v0y=t.v[0].getY();
	double v1x=t.v[1].getX();double v1y=t.v[1].getY();
	double v2x=t.v[2].getX();double v2y=t.v[2].getY();
	
	for(int ps=0;ps<7;ps++){
		double px=lambda(0,PtsRef[ps])*v0x+lambda(1,PtsRef[ps])*v1x+lambda(2,PtsRef[ps])*v2x;
		double py=lambda(0,PtsRef[ps])*v0y+lambda(1,PtsRef[ps])*v1y+lambda(2,PtsRef[ps])*v2y;
		R2 p(px,py);
		points[ps]=p;
	}

}

//Fct qui calcule les caractéristiques dans le second membre
void CalculCaracteristique(Mesh2d & Th,double alpha,vector<double> xprec,int n,double * b){
	int nt=Th.nbt;
	assert(xprec.size()>0);
	//pts de quadrature
	double ptint1=(6-sqrt(15))/21;
	double ptint2=(9-sqrt(15)*2)/21;
	double ptint3=(6+sqrt(15))/21;
	double ptint4=(9+sqrt(15)*2)/21;
	R2 p1(1./3,1./3),p2(ptint1,ptint1),p3(ptint1,ptint4),p4(ptint4,ptint1),p5(ptint3,ptint3),p6(ptint3,ptint2),p7(ptint2,ptint3);
	R2 PtsRef[7]={p1,p2,p3,p4,p5,p6,p7};//Points de quadrature

	double poids1=(155-sqrt(15))/1200;//Poids des points de quadrature
	double poids2=(155+sqrt(15))/1200;
	double Poids[7]={0.225,poids1,poids1,poids1,poids2,poids2,poids2};
	double * u1pk=new double[6];double * u2pk=new double[6];
	double * PointCaractX=new double[7];double * PointCaractY=new double[7];
	R2 * Point=new R2[7];
	double * u1pInterp=new double[6]; double * u2pInterp=new double[6];
	double * u1pInterp2=new double[7]; double * u2pInterp2=new double[7];

	int i;
	int vois;
	double c=0;

	double * phi=new double[7]; //Fonction test
	bool boolRecup=1; //1 = si on est toujours dans le domaine en remontant les caracteristiques

	for(int k=0; k<nt;k++){ //boucle sur les triangles
	
		double areak=Th.t[k].area;
		PointK(Th.t[k],PtsRef, Point); //transforme les points PtsRef en points dans le triangle k
		recup(Th.t[k],xprec,u1pk,u2pk,n);//on recupere dans le triangle k les vitesses u1pk et u2pk
		for(int ps=0;ps<7;ps++){ //boucle sur les points de quadratures
			boolRecup=1; 
			u1pInterp[ps]=vitesseInterpolee(u1pk,PtsRef[ps]);
			u2pInterp[ps]=vitesseInterpolee(u2pk,PtsRef[ps]);
			assert(alpha>0);assert(Th.voisins[k].size()>0);
			assert(u1pInterp[ps]<3 && u2pInterp[ps]<3);
			PointCaractX[ps]= Point[ps].x-(1./alpha)*u1pInterp[ps]; //Position du point de quadrature au pas précédent
			PointCaractY[ps]=Point[ps].y-(1./alpha)*u2pInterp[ps];
			R2 PointCaract(PointCaractX[ps],PointCaractY[ps]);
			R2 PointCaractRef;//Point caracteristique dans le triangle de reference
			vois=RecupVoisins(Th,k,PointCaract,PointCaractRef); //triangle auquel appartient PointCaract et PointCaractRef = dans triangle ref
			//assert(vois>=0);
			//cout<<"Ref "<<PtsRef[ps]<<" P "<< Point[ps]<<" PtcaracRef "<<PointCaractRef<<" Ptcarac "<<PointCaract<<endl;
			if(vois<0){ //si on sort du domaine...
				//cout<<PointCaract<<endl;
				if(PointCaract.x<0){
					PointCaract.x=0;
					if(PointCaract.y<0.5)
						PointCaract.y=0.5;
					else if(PointCaract.y<=1)
						PointCaract.y=PointCaract.y;
					else
						PointCaract.y=1;
					Vertex PointCaractBord;PointCaractBord.setX(PointCaract.x);PointCaractBord.setY(PointCaract.y);
					u1pInterp2[ps]=g(PointCaractBord,10);
					u2pInterp2[ps]=0;
					boolRecup=0;
				}
				else if(PointCaract.x<=10 ){
					u1pInterp2[ps]=0;
					u2pInterp2[ps]=0;
					boolRecup=0;
				}
				else{
					PointCaract.x=10; //bord droit du domaine
					if(PointCaract.y<=0 || PointCaract.y>=1){
						u1pInterp2[ps]=0;
						u2pInterp2[ps]=0;
						boolRecup=0;
					}
					else{
					vois=find_triangle(PointCaract,Th);//on trouve le triangle auquel appartient PointCaract
					assert(vois>=0);
					}
				}
			}
			//assert(boolRecup==1);
			if(boolRecup==1){
				recup(Th.t[vois],xprec,u1pInterp,u2pInterp,n);//calcul de u1pk et u2pk
				u1pInterp2[ps]=vitesseInterpolee(u1pInterp,PointCaractRef);
				u2pInterp2[ps]=vitesseInterpolee(u2pInterp,PointCaractRef);
			}
		}

		for(int il=0;il<15;il++){
			if(il<6){
				i=Th(k,il);
				for(int ps=0;ps<7;ps++){	
					phi[ps]=	Phi(il, PtsRef[ps]);
				}
			}
			else if(il<12){
				i=Th(k,il-6)+n;
				for(int ps=0;ps<7;ps++){		
					phi[ps]=	Phi(il-6,PtsRef[ps]);
				}
			}
			else{
				i=Th(k,il-12)+2*n;
			}
		c=0;
		if(i<n){
			for(int ps=0;ps<7;ps++){
				c+=Poids[ps]*phi[ps]*u1pInterp2[ps];
			}
			b[i]+=alpha*areak*c;
		}
		else if(i<2*n){
			for(int ps=0;ps<7;ps++){
				c+=Poids[ps]*phi[ps]*u2pInterp2[ps];
			}
			b[i]+=alpha*areak*c;
		}
		else
			b[i]=0; //rien sur la pression
		}
	}
	delete[] Point; delete[] u1pk; delete[] u2pk; delete[] u1pInterp; delete[] u2pInterp; delete[] u1pInterp2;delete[] u2pInterp2;delete[] phi; delete[] PointCaractX; delete[] PointCaractY; 
}

