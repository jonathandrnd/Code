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
using namespace std;
int main(){
	//freopen("in.txt","r",stdin);
	//freopen("out.txt","w",stdout);
	
	int n;

	while(cin>>n){
		set<long long>S;

		string auxs;
		getline(cin,auxs);
		getline(cin,auxs);
		
		for(int i=0;i<(int)auxs.size();i++)
			if(auxs[i]==',')
				auxs[i]=' ';
		
		istringstream is(auxs);
		
		vector<long long>v;
		long long aux;
		while(is>>aux){
			v.push_back(aux);
		}
		
		/*
		if((int)v.size()!=n)
			assert(0);
		if((int)v.size()>100 || n>100)
			assert(0);
		if(n<6)
			assert(0);
		for(int i=0;i<(int)v.size();i++)
			if(v[i]<0 || v[i]>1000000000LL)
				assert(0);
		*/
		vector<int>p;
		sort(v.begin(),v.end());
		reverse(v.begin(),v.end());
		
		for(int i=0;i<6;i++)p.push_back(i);
		
		long long dev=0;
		do{
			vector<long long>x;
			x.push_back(v[p[0]]);
			x.push_back(v[p[1]]);
			x.push_back(v[p[2]]);
			x.push_back(v[p[3]]);
			x.push_back(v[p[4]]);
			x.push_back(v[p[5]]);
			vector<long long>aux=x;
			
			while(x.size()!=2){
				vector<long long>h;
				for(int i=0;i+1<(int)x.size();i++)
					h.push_back(x[i]+x[i+1]);
				x=h;
			}
			
			S.insert(x[0]*x[1]);
			dev=max(dev,x[0]*x[1]);
			
			
			
		}while(next_permutation(p.begin(),p.end()));
		
		vector<long long>hamlet(S.begin(),S.end());
		/*
		if(hamlet[(int)hamlet.size()-1]==hamlet[(int)hamlet.size()-2])
			assert(0);
		*/
		if(n==9){
			cout<<hamlet[(int)hamlet.size()-1]<<endl;
		}else{
			cout<<hamlet[(int)hamlet.size()-2]<<endl;
		}
		
	}
	return 0;
}
