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
		x_=coord[0];
		y_=coord[1];
		NumGlobal_=ind;//+offset;
		lab_=label;
	}
	R getm(Vertex & Q){
		return (sqrt((x_-Q.getX())*(x_-Q.getX())+(y_-Q.getY())*(y_-Q.getY())));
	}
	void setX(double x){x_=x;}
	void setY(double y){y_=y;}
	double getX(){return x_;}
	double getY(){return y_;}
	void setNum(int ng){NumGlobal_=ng;}
	void setLab(Label l){lab_=l;}
	Label getLab(){return lab_;}
	int getNum(){return NumGlobal_;}
	Vertex operator=(Vertex v2){NumGlobal_=v2.getNum();x_=v2.getX();y_=v2.getY();lab_=v2.getLab();return *this;}
	Vertex operator=(R2 & v2){x_=v2.x;y_=v2.y;lab_=0;return *this;}
 
private:
	int NumGlobal_;
	R x_;
	R y_;
	Label lab_;
};

static int nt=0;
static int ne=0;

class Triangle {
	public:
		int numTri;
		Vertex v[3];
		Vertex mil[3];
		//int MilNum[3];
		double area;
		Triangle(){v[0]=(v[1]=(v[2]));}
		double build(vector<Vertex> vtot,int * I,int offset){
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
	//	Vertex operator[](int i){return v[i];}
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
			//cout<<"lab"<<lab<<endl;
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
	//int ReturnEdge(Vertex v0,Vertex v1);
	int PointsMil();
	int operator()(int k, int i); // num global du sommet/milieu i du triangle k
	Triangle operator[](int k)const;
private:
 
};
