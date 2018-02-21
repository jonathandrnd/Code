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
vector<pair<int,int> >v;

int toi(string s){
	int dev=0;
	for(int i=0;i<(int)s.size();i++)
		dev=dev*10+(s[i]-'0');
	return dev;
}

int main(){
	freopen("in.txt","r",stdin);
	freopen("out.txt","w",stdout);
	
	int n;
	while(cin>>n){
		v=vector<pair<int,int> >(n,make_pair(0,0));
		string s;
		getline(cin,s);
		for(int i=0;i<n;i++){
			getline(cin,s);
			for(int j=0;j<(int)s.size();j++)if(s[j]==',' || s[j]==';')s[j]=' ';
			istringstream is(s);
			string a,b;
			is>>a>>b;
			
			char ch1=a[(int)a.size()-1];
			char ch2=b[(int)b.size()-1];
			int val1=toi(a.substr(0,(int)a.size()-1));
			int val2=toi(b.substr(0,(int)b.size()-1));
			if(ch1=='N'){
				v[i].first=val1;
			}else{
				v[i].first=360-val1;
			}
			
			if(ch2=='E'){
				v[i].second=val2;
			}else{
				v[i].second=360-val2;
			}	
		}
		
		double sum=0;
		double pi=acos(-1)/180.0;
		
		for(int i=0;i+1<n;i++){
			double lat1=v[i].first*pi;
			double lat2=v[i+1].first*pi;
			double lon1=v[i].second*pi;
			double lon2=v[i+1].second*pi;
			double x=2*asin(sqrt((  sin((lat1-lat2)/2.0))* sin((lat1-lat2)/2.0)) +cos(lat1)*cos(lat2)*(  sin((lon1-lon2)/2.0))*sin((lon1-lon2)/2.0))   ;
			sum+=2*asin(sqrt((  sin((lat1-lat2)/2.0))* sin((lat1-lat2)/2.0)) +cos(lat1)*cos(lat2)*(  sin((lon1-lon2)/2.0))*sin((lon1-lon2)/2.0))   ;
		}	
		
		int h=(int)(sum*6400+0.5);
		cout<<h<<endl;
	}
	
	
    return 0;
}
