#include "mesh_gmsh.hpp"
#include <cassert>
#include <map>

using namespace std;

//int static num=0;

Mesh2d::Mesh2d(const char * filename)
{
  std::ifstream  f(filename); 
  assert( f); 
  double coord[3] ; 
	int ind;
	string ligne;
	while(getline(f,ligne)&& ligne!="$Nodes"){};//passer à la ligne suivante
	f>>nv;

	for(int i=0;i<nv;++i)
  { 
		f>>ind;
		for(int l=0;l<3;l++){
		f>>coord[l];
		}
		Vertex Vt;
    Vt.build(coord,ind,-1); 
    assert( f.good());
		v.push_back(Vt);
  }
	
	f>>ligne;
	f>>ligne;
	int NbEl;
	f>>NbEl;
	vector<double> area;
	int I[8];
	for(int i=0;i<NbEl;i++){
		for(int k=0;k<7;k++){
			f>>I[k];
		}
		if (I[1]==2){
			f>>I[7];
			Triangle Trg;
			double areak=Trg.build(v,I,-1);
			area.push_back(areak);
			t.push_back(Trg);
		}
		else{
			Edge Edg;
			Edg.build(v,I,-1);
			e.push_back(Edg);
		}
	}
	nbt=t.size();
	this->voisins.resize(nbt);
	nbe=e.size();
	double a=0;
	for(int i=0;i<(int)area.size();i++){
		a=a+area[i];
	}
	cout<<"nbv: "<<nv<<" nbt: "<<nbt<<" ne: "<<nbe<<" aire: "<<a<<endl;
}


int Mesh2d::operator()(int k, int i){
	if(i<3)
		return (this->t[k].v[i].getNum());
	else
		return (this->t[k].mil[i-3].getNum());
}// num global du sommet/milieu i du triangle k

Triangle Mesh2d::operator[](int k)const{
	return (this->t[k]);
}

int Mesh2d::PointsMil(){//num des pts milieux

	map< pair<int,int>,int> ME; //s1,s2,label (edge)

	for(int k=0;k<nbe;k++){//aretes du bord
			int s1 = e[k].v[0].getNum();
			int s2 = e[k].v[1].getNum();
			if( s2>s1) 
				swap(s1,s2);
			pair<int,int> key(s1,s2);
			int val=e[k].lab.lab;
			ME[key]=val;
			v[s1].setLab(val);
			v[s2].setLab(val);
	}

	int n=nv;
	bool cree=false;
 	map< pair<int,int>, pair<int,int>> M;
	for(int k=0;k<nbt;k++){
		for(int a=0; a< 3; ++a){// Arete oppose au sommet
			cree=false;
			int s1 = t[k].v[(a+1)%3].getNum();
			int s2 = t[k].v[(a+2)%3].getNum();
			if( s2>s1) 
				std::swap(s1,s2);
			double x0=v[s1].getX();
			double x1=v[s2].getX();
			double y0=v[s1].getY();
			double y1=v[s2].getY();
			pair<int,int> key(s1,s2);
			Vertex milieu;
			milieu.setX((x0+x1)/2);
			milieu.setY((y0+y1)/2);
			if(M.find(key)== M.end()){//si on ne trouve pas dans la map on ajoute le milieu
				M[key].first=n++;
				M[key].second=k;
				cree=true;
			}
			milieu.setNum(M[key].first);
			t[k].mil[a]=milieu;
			if(ME.find(key)!=ME.end()){//si on trouve la clé correspond au bord du domaine on ajoute le label
				if(ME[key]==30){//bord sortie pour l'interpolation dans la methode des caracteristiques
					triangleSortie.push_back(k);
				}
				t[k].v[(a+1)%3].setLab(ME[key]);
				t[k].v[(a+2)%3].setLab(ME[key]);
				t[k].mil[a].setLab(ME[key]);
				milieu.setLab(ME[key]);
			}
			if(cree==true)//Si le point milieu vient d être créer, on l'ajoute au vecteur des points
				this->v.push_back(milieu);
		}
		for(int a=0;a<3;a++){
			int num=t[k].v[a].getNum();
			vector<int> tri=v[num].getTri();
			for(unsigned int i=0;i<tri.size();i++){
				if(tri[i]!=k)
					voisins[k].push_back(tri[i]);
			}
		}
	}
	std::vector<int>::iterator it1;
	for(int k=0;k<nbt;k++){
		sort(voisins[k].begin(),voisins[k].end());
		it1 = std::unique (voisins[k].begin(), voisins[k].end()); 
  	voisins[k].resize( std::distance(voisins[k].begin(),it1) );
	}
	/*for(int k=0;k<nbt;k++){
		cout<<"triangle  "<<k<<": ";
		for(int i=0;i<voisins[k].size();i++){
			cout<<voisins[k][i]<<" ";
		}
		cout<<endl;
	}*/
	//cout<<"nb de triangles au bord"<<triangleSortie.size()<<endl;
	return n;
}

