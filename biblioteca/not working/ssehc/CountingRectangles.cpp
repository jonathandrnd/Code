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
using namespace std;

int f(int n){
	
	set<vector<pair<int,int> > >S;
	
	for(int i=0;i<=n;i++)
		for(int j=0;j<=n;j++)
			for(int k=0;k<=n;k++)
				for(int m=0;m<=n;m++){
					int dx=k-i;
					int dy=m-j;
					if(dx==0 && dy==0)continue;
					int auxx=dx;
					int auxy=dy;
					dx=-auxy;
					dy= auxx;
					if(dx<0){
						dx=-dx;
						dy=-dy;
					}
					
					int ham=__gcd(abs(dx),abs(dy));
					dx/=ham;
					dy/=ham;
					
					
					int wi=i;
					int wj=j;
					int wk=k;
					int wm=m;
					
					while(true){
						wi+=dx;
						wj+=dy;
						
						wk+=dx;
						wm+=dy;
						
						if(wi>=0 && wi<=n && wj>=0 &&wj<=n &&
						   wk>=0 && wk<=n && wm>=0 &&wm<=n ){
								vector<pair<int,int> >v;
								v.push_back(make_pair(i,j));
								v.push_back(make_pair(k,m));
								
								v.push_back(make_pair(wi,wj));
								v.push_back(make_pair(wk,wm));
								sort(v.begin(),v.end());
								
								S.insert(v);
						   }else break;
					}
				}
	
	return S.size();
}

int main(){
	//freopen("in.txt","r",stdin);
	//freopen("out.txt","w",stdout);

int c[56];	
c[0]=0;
c[1]=1;
c[2]=10;
c[3]=44;
c[4]=130;
c[5]=313;
c[6]=640;
c[7]=1192;
c[8]=2044;
c[9]=3305;
c[10]=5078;
c[11]=7524;
c[12]=10750;
c[13]=14993;
c[14]=20388;
c[15]=27128;
c[16]=35448;
c[17]=45665;
c[18]=57922;
c[19]=72636;
c[20]=89970;
c[21]=110297;
c[22]=133976;
c[23]=161440;
c[24]=192860;
c[25]=228857;
c[26]=269758;
c[27]=316012;
c[28]=367974;
c[29]=426417;
c[30]=491468;
c[31]=564120;
c[32]=644640;
c[33]=733633;
c[34]=831674;
c[35]=939292;
c[36]=1056962;
c[37]=1186057;
c[38]=1326880;
c[39]=1479968;
c[40]=1645820;
c[41]=1825945;
c[42]=2020350;
c[43]=2230660;
c[44]=2457046;
c[45]=2700305;
c[46]=2961708;
c[47]=3242368;
c[48]=3542232;
c[49]=3863217;
c[50]=4205666;

	int n;
	while(cin>>n){
		cout<<c[n]<<endl;
	}

/*
	
	for(int tam=1;tam<=50;tam++){
		long long val=f(tam);
		int i=tam;
		cout<<"c["<<tam<<"]="<<val<<";"<<endl;
		//cout<<tam<<" "<<val<<" "<<val-( (i*(i+1))/2 )*( (i*(i+1))/2 )<<endl;
		
	}
*/
	
	return 0;
}
