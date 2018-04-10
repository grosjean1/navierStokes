#include "mesh.hpp"
#include <cassert>
#include <map>
#include <algorithm> 

using namespace std;

//int static num=0;

Mesh2d::Mesh2d(const char * filename)
{
  std::ifstream  f(filename); 
  assert( f); 
  double coord[3] ; 
	int inu;
	int ind;
	string ligne;
	f>>nv;f>>nbt;f>>nbe;
	for(int i=0;i<nv;++i)
  { 
		ind=i;
		for(int l=0;l<2;l++){
		f>>coord[l];
		}
		f>>inu;
		Vertex Vt;
    Vt.build(coord,ind); 
    assert( f.good());
		v.push_back(Vt);
  }
	vector<double> area;
	int I1[4];int I2[3];
	for(int i=0;i<nbt;i++){
		for(int k=0;k<3;k++){
			f>>I1[k];
		}
		f>>inu;
		Triangle Trg;
		double areak=Trg.build(v,I1,-1);
		area.push_back(areak);
		t.push_back(Trg);
	}
	for(int i=0;i<nbe;i++){
		for(int k=0;k<3;k++){
			f>>I2[k];
		}
		Edge Edg;
		Edg.build(v,I2,-1);
		e.push_back(Edg);
	}
	this->voisins.resize(nbt);
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
			//cout<<"s1 "<<key.first<<" s2 "<<key.second<<" val "<<val<<endl;
	}

	//cout<<"edges"<<nbe<<endl;
	int n=nv;
	bool cree=false;
 	map< pair<int,int>, pair<int,int>> M;
	bool VoisinExiste=false;
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
				for(map<pair<int,int>,pair<int,int>>::iterator it=M.begin();it!=M.end();++it){			
					if(it->first.first==s1||it->first.first==s2||it->first.second==s1||it->first.second==s2){
					//	cout<<"k "<<k<<" s1 "<<s1<<" s2 "<<s2<<" "<<it->first.first<<" "<<it->first.second<<" "<<it->second.second<<endl;
						if(it->second.second!=k){
							this->voisins[k].push_back(it->second.second);
							this->voisins[it->second.second].push_back(k);
						}
					}
				}
				std::vector<int>::iterator it1;
				for(int k=0;k<nbt;k++){
					sort(voisins[k].begin(),voisins[k].end());
				  it1 = std::unique (voisins[k].begin(), voisins[k].end()); 
  				voisins[k].resize( std::distance(voisins[k].begin(),it1) );
				}
			/*else{
				//cout<<m<<" "<<k<<" "<<M[key]<<" "<<s1<<" "<<s2<<endl;
				this->voisins[k].push_back(M[key].second);
				this->voisins[M[key].second].push_back(k);
			}*/
			milieu.setNum(M[key].first);
			t[k].mil[a]=milieu;
			if(ME.find(key)!=ME.end()){//si on trouve la clé correspond au bord du domaine on ajoute le label
				t[k].v[(a+1)%3].setLab(ME[key]);
				t[k].v[(a+2)%3].setLab(ME[key]);
				t[k].mil[a].setLab(ME[key]);
				milieu.setLab(ME[key]);
			}
			if(cree==true)//Si le point milieu vient d être créer, on l'ajoute au vecteur des points
				this->v.push_back(milieu);
		}
	}

	//cout<<"taille des sommets!" << v.size()<<endl;
	return n;
}

