#include "R2.hpp"
#include <cassert>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <iostream>
#include <vector>

using namespace std;

class Label {
public:
  int lab; 
  int  OnGamma() const { return lab;}
  Label(int l=0) : lab(l) {}
};

class Vertex :public R2,public  Label
{
  typedef R2 Rd;
 	public:
  Vertex() : R2(),Label(){};
	void build(double * coord, int ind,Label label=0){
		x=coord[0];
		y=coord[1];
		NumGlobal_=ind;
		lab=label.lab;
	}
	R getm(Vertex & Q){
		return (sqrt((x-Q.getX())*(x-Q.getX())+(y-Q.getY())*(y-Q.getY())));
	}
	void setX(double xx){x=xx;}
	void setY(double yy){y=yy;}
	double getX(){return x;}
	double getY(){return y;}
	void setNum(int ng){NumGlobal_=ng;}
	void setLab(Label l){lab=l.lab;}
	Label getLab(){return lab;}
	int getNum(){return NumGlobal_;}
	void ajoutTri(int k){
		tri.push_back(k);
	}
	vector<int> getTri(){return tri;}
private:
	int NumGlobal_;
	vector<int> tri;
};


static int nt=0;
static int ne=0;

class Triangle {
	public:
		int numTri;
		Vertex v[3];
		Vertex mil[3];
		double area;
		Triangle(){v[0]=(v[1]=(v[2]));}
		double build(vector<Vertex> & vtot,int * I,int offset){
			numTri=nt++;
			v[0]=vtot[I[0]+offset];
			v[1]=vtot[I[1]+offset];
			v[2]=vtot[I[2]+offset];
			double a=v[0].getm(v[1]);
			double b=v[1].getm(v[2]);
			double c=v[2].getm(v[0]);
			double p=(a+b+c)/2; //semi-perimetre
			area=sqrt(p*(p-a)*(p-b)*(p-c));
			return area;
		}
};


class Edge {//bord du domaine

	public:
		int numEdg;
		Label lab;
		Vertex v[2];
		Edge(){v[0]=(v[1]);}

		void build(vector<Vertex> vtot, int * I,int offset){
			numEdg=ne++;
			v[0]=vtot[I[0]+offset];
			v[1]=vtot[I[1]+offset];
			lab=I[2];
		}
	private:
};


class Mesh2d 
{
public:
	vector<vector <int>> voisins;
  int nv,nbt,nbe;
	vector<Vertex> v;
  vector<Triangle> t;
	vector<Edge> e;
  Mesh2d(const char *  filename);
  ~Mesh2d() {};
	int PointsMil();
	int operator()(int k, int i); // num global du sommet/milieu i du triangle k
	Triangle operator[](int k)const;
	vector<int> triangleSortie;
private:
  Mesh2d(const Mesh2d &);
  void operator=(const Mesh2d&);
};
