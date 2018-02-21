#include <vector>
#include <list>
#include <map>
#include <set>
#include <queue>
#include <deque>
#include <stack>
#include <bitset>
#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <string.h>
#include <stdio.h>
#define all(v) (v).begin(),(v).end()
#define rall(v) (v).rbegin(),(v).rend()
using namespace std;  // H A M L E T


string tos(int t){stringstream st;st<<t;return st.str();}
int toi(string s){istringstream is(s);int t;is>>t;return t;}
vector<string >Bcambio;
vector<vector<string> >moves;
string turno;
int dnx[8]={-1,-1,1,1,-2,-2,2,2};
int dny[8]={2,-2,2,-2,1,-1,1,-1};	

string namechange(vector<string>c,int x,int y,int X,int Y,bool comido){
	string dev="";
	if(c[x][y]!='P' && c[x][y]!='p')dev+=c[x][y];
	dev+=('a'+y);
	dev+=('0'+(8-x));
	if(c[x][y]!='P' && c[x][y]!='p')dev+=c[x][y];
	dev+=('a'+Y);
	dev+=('0'+(8-X));
	return dev;
}

bool orden(string a,string b){
	int r1=a[1]-'0';
	int p1=a[0]-'0';
	
	int r2=b[1]-'0';
	int p2=b[0]-'0';
	
	
	if(r1!=r2)return r1>r2;
	if(p1!=p2)return p1<p2;
	
	r1=a[3]-'0';
	p1=a[2]-'0';
	
	r2=b[3]-'0';
	p2=b[2]-'0';

	if(r1!=r2)return r1>r2;
	return p1<p2;	
}


bool llegaNegro(int x,int y,int dx,int dy,int X,int Y,vector<string>c){
	if(x<0 || y<0 || x>=8 || y>=8)return 0;
	if(x==X && y==Y)return 1;
	if(c[x][y]>='A' && c[x][y]<='Z')return 0;
	if(c[x][y]>='a' && c[x][y]<='z')return 0;
	return llegaNegro(x+dx,y+dy,dx,dy,X,Y,c);
}	

bool negrohackeado(vector<string>c){
	int x,y;
	
	for(int i=0;i<8;i++)
		for(int j=0;j<8;j++)
			if(c[i][j]=='k'){
				x=i;y=j;
			}
	
	for(int i=0;i<8;i++)
		for(int j=0;j<8;j++){
			if(c[i][j]>='A' && c[i][j]<='Z'){
				char ch=c[i][j];
				if(ch=='R'){
					if( llegaNegro(i-1,j,-1,0,x,y,c) || llegaNegro(i+1,j,1,0,x,y,c) || 
						llegaNegro(i,j-1,0,-1,x,y,c) || llegaNegro(i,j+1,0,1,x,y,c) ) return 1;
				}
				
				if(ch=='B'){
					if( llegaNegro(i+1,j+1,1,1,x,y,c) || llegaNegro(i+1,j-1,1,-1,x,y,c) ||
						llegaNegro(i-1,j+1,-1,1,x,y,c) || llegaNegro(i-1,j-1,-1,-1,x,y,c) )return 1;
				}
	
				if(ch=='Q'){
					if( llegaNegro(i-1,j,-1,0,x,y,c) || llegaNegro(i+1,j,1,0,x,y,c) || 
						llegaNegro(i,j-1,0,-1,x,y,c) || llegaNegro(i,j+1,0,1,x,y,c) || 
						llegaNegro(i+1,j+1,1,1,x,y,c) || llegaNegro(i+1,j-1,1,-1,x,y,c) ||
						llegaNegro(i-1,j+1,-1,1,x,y,c) ||  llegaNegro(i-1,j-1,-1,-1,x,y,c) ) return 1;
				}
	
				if(ch=='N'){
					for(int k=0;k<8;k++){
						int X=i+dnx[k];
						int Y=j+dny[k];
						if(X==x && Y==y)return 1;
					}
				}
				
				if(ch=='K'){
					if(abs(i-x)<=1 && abs(j-y)<=1)return 1;
				}
				
				if(ch=='P'){
					if(i-1==x && abs(j-y)==1 )return 1;
				}	
			}
		}
	return 0;
}

bool llegaBlanco(int x,int y,int dx,int dy,int X,int Y,vector<string>c){
	if(x<0 || y<0 || x>=8 || y>=8)return 0;
	if(x==X && y==Y)return 1;
	if(c[x][y]>='A' && c[x][y]<='Z')return 0;
	if(c[x][y]>='a' && c[x][y]<='z')return 0;
	return llegaBlanco(x+dx,y+dy,dx,dy,X,Y,c);
}	


bool blancohackeado(vector<string>c){
	int x,y;
	
	for(int i=0;i<8;i++)
		for(int j=0;j<8;j++)
			if(c[i][j]=='K'){
				x=i;y=j;
			}
	
	for(int i=0;i<8;i++)
		for(int j=0;j<8;j++){
			if(c[i][j]>='a' && c[i][j]<='z'){
				char ch=c[i][j];
				if(ch=='r'){
					if( llegaBlanco(i-1,j,-1,0,x,y,c) || llegaBlanco(i+1,j,1,0,x,y,c) || 
						llegaBlanco(i,j-1,0,-1,x,y,c) || llegaBlanco(i,j+1,0,1,x,y,c) ) return 1;
				}
				
				if(ch=='b'){
					if( llegaBlanco(i+1,j+1,1,1,x,y,c) || llegaBlanco(i+1,j-1,1,-1,x,y,c) ||
						llegaBlanco(i-1,j+1,-1,1,x,y,c) ||  llegaBlanco(i-1,j-1,-1,-1,x,y,c) )return 1;
				}
	
				if(ch=='q'){
					if( llegaBlanco(i-1,j,-1,0,x,y,c) || llegaBlanco(i+1,j,1,0,x,y,c) || 
						llegaBlanco(i,j-1,0,-1,x,y,c) || llegaBlanco(i,j+1,0,1,x,y,c) || 
						llegaBlanco(i+1,j+1,1,1,x,y,c) || llegaBlanco(i+1,j-1,1,-1,x,y,c) ||
						llegaBlanco(i-1,j+1,-1,1,x,y,c) ||  llegaBlanco(i-1,j-1,-1,-1,x,y,c) ) return 1;
				}
	
				if(ch=='n'){
					for(int k=0;k<8;k++){
						int X=i+dnx[k];
						int Y=j+dny[k];
						if(X==x && Y==y)return 1;
					}
				}
				
				if(ch=='k'){
					if(abs(i-x)<=1 && abs(j-y)<=1)return 1;
				}
				
				if(ch=='p'){
					if(i+1==x && abs(j-y)==1 )return 1;
				}	
			}
		}
	return 0;
}

void Bposible(int x,int y,int dx,int dy,vector<string>c,int X,int Y){
	if(x<0 || y<0 || x>=8 || y>=8)return;
	if(c[x][y]>='A' && c[x][y]<='Z')return;
	vector<string>aux=c;
	aux[x][y]=c[X][Y];
	aux[X][Y]=' ';
		
	if(c[x][y]>='a' && c[x][y]<='z'){
		if(!blancohackeado(aux)){
			Bcambio.push_back( namechange(c,X,Y,x,y,1) );
			moves.push_back(aux);
		}
		return ;		
	}
	
	if(!blancohackeado(aux)){
		Bcambio.push_back( namechange(c,X,Y,x,y,0)  );
		moves.push_back(aux);
	}
	Bposible(x+dx,y+dy,dx,dy,c,X,Y);
}

void Bsolouno(int x,int y,vector<string>c,int X,int Y){
	if(x<0 || y<0 || x>=8 || y>=8)return;
	if(c[x][y]>='A' && c[x][y]<='Z')return;
	vector<string>aux=c;
	aux[x][y]=c[X][Y];
	aux[X][Y]=' ';
	
	
	if(c[x][y]>='a' && c[x][y]<='z'){
		if(!blancohackeado(aux)){
			Bcambio.push_back( namechange(c,X,Y,x,y,1)  );
			moves.push_back(aux);
		}
		return;
	}
	
	if(!blancohackeado(aux)){		
		Bcambio.push_back( namechange(c,X,Y,x,y,0)  );
		moves.push_back(aux);
	}
}

vector<vector<string> > mover(vector<string> c,int x,int y){
	string move="";
	vector<vector<string> > ans;
	char ch=c[x][y];
	moves.clear();
	
	if(ch=='P'){
		//if(x==7 && c[5][y]==' ')
		//	v.push_back(make_pair(5,y));
		if(x==6 && c[x-2][y]==' ' && c[x-1][y]==' '){
			vector<string>aux=c;
			aux[x][y]=' ';
			aux[x-2][y]='P';
			if(!blancohackeado(aux)){
				Bcambio.push_back(namechange(c,x,y,x-2,y,0) );
				moves.push_back(aux);
			}
		}
		
		
		if(x>0 && c[x-1][y]==' '){
			vector<string>aux=c;
			aux[x][y]=' ';
			aux[x-1][y]='P';
			if(!blancohackeado(aux)){
				Bcambio.push_back(namechange(c,x,y,x-1,y,0) );
				moves.push_back(aux);
			}
		}
			
		if(x>0 && y>0 && c[x-1][y-1]>='a' && c[x-1][y-1]<='z'){
			vector<string>aux=c;
			aux[x][y]=' ';
			aux[x-1][y-1]='P';
			if(!blancohackeado(aux)){
				Bcambio.push_back(namechange(c,x,y,x-1,y-1,1) );
				moves.push_back(aux);
			}
		}
			
		if(x>0 && y<7 && c[x-1][y+1]>='a' && c[x-1][y+1]<='z'){
			vector<string>aux=c;
			aux[x][y]=' ';
			aux[x-1][y+1]='P';
			if(!blancohackeado(aux)){
				Bcambio.push_back(namechange(c,x,y,x-1,y+1,1) );
				moves.push_back(aux);
			}
		}
	}
	
	return moves;
}




vector<vector<string> > blancas(vector<string>c){
	vector<vector<string> >dev;
	
	for(int i=0;i<8;i++){
		for(int j=0;j<8;j++){
			if(c[i][j]>='A' && c[i][j]<='Z'){
				vector<vector<string> >x  =mover(c,i,j);
				for(int ii=0;ii<(int)x.size();ii++)dev.push_back(x[ii]);			
			}
		}
	}
	
	return dev;
}



vector<vector<string> > movernegras(vector<string> c,int x,int y){
	string move="";
	vector<vector<string> > ans;
	char ch=c[x][y];
	moves.clear();
	
	if(ch=='p'){
		//if(x==7 && c[5][y]==' ')
		//	v.push_back(make_pair(5,y));
		if(x==1 && c[x+2][y]==' ' && c[x+1][y]==' '){
			vector<string>aux=c;
			aux[x][y]=' ';
			aux[x+2][y]='p';
			if(!negrohackeado(aux)){
				Bcambio.push_back(namechange(c,x,y,x+2,y,0) );
				moves.push_back(aux);
			}
		}
		
		
		if(x<7 && c[x+1][y]==' '){
			vector<string>aux=c;
			aux[x][y]=' ';
			aux[x+1][y]='p';
			if(!negrohackeado(aux)){
				Bcambio.push_back(namechange(c,x,y,x+1,y,0) );
				moves.push_back(aux);
			}
		}
			
		if(x<7 && y>0 && c[x+1][y-1]>='A' && c[x+1][y-1]<='Z'){
			vector<string>aux=c;
			aux[x][y]=' ';
			aux[x+1][y-1]='p';
			if(!negrohackeado(aux)){
				Bcambio.push_back(namechange(c,x,y,x+1,y-1,1) );
				moves.push_back(aux);
			}
		}
			
		if(x<7 && y<7 && c[x+1][y+1]>='A' && c[x+1][y+1]<='Z'){
			vector<string>aux=c;
			aux[x][y]=' ';
			aux[x+1][y+1]='p';
			if(!negrohackeado(aux)){
				Bcambio.push_back(namechange(c,x,y,x+1,y+1,1) );
				moves.push_back(aux);
			}
		}
	}
	
	return moves;
}




vector<vector<string> > negras(vector<string>c){
	vector<vector<string> >dev;
	
	for(int i=0;i<8;i++){
		for(int j=0;j<8;j++){
			if(c[i][j]>='a' && c[i][j]<='z'){
				vector<vector<string> >x  =movernegras(c,i,j);
				for(int ii=0;ii<(int)x.size();ii++)dev.push_back(x[ii]);			
			}
		}
	}
	
	return dev;
}


vector<string> check(vector<string>c){
	
	//mover negras
	vector<vector<string> >chum;
	if(turno=="b"){
		chum=negras(c);
	}else{
		chum=blancas(c);
	}
	
	
	return Bcambio;
}

int main(){
    //freopen("in.txt","r",stdin);
    //freopen("out.txt","w",stdout);
	
	string s;
	
	while(cin>>s){
		cin>>turno;
		if(s=="-1")break;
		Bcambio.clear();
		moves.clear();
		for(int i=0;i<(int)s.size();i++)if(s[i]=='/')s[i]=' ';
		
		istringstream is(s);
		string aux;
		int cont=0;
		vector<string>c(8,"        ");
		
		while(is>>aux){
			int col=0;
			for(int i=0;i<(int)aux.size();i++){
				if(aux[i]>='0' && aux[i]<='9'){
					col+=(aux[i]-'0');
				}else{
					c[cont][col++]=aux[i];
				}
			}
			cont++;
		}
		
		/*
		for(int i=0;i<8;i++){
			for(int j=0;j<8;j++)
				cout<<c[i][j];
			cout<<endl;
		}*/
		
		vector<string>chum=check(c);
		cout<<"[";
		sort(chum.begin(),chum.end(),orden);
		
		for(int i=0;i<(int)chum.size();i++){
			cout<<chum[i];
			if(i+1!=(int)chum.size())cout<<", ";
		}
		
		cout<<"]"<<endl;
	}
	
	
	return 0;
}

