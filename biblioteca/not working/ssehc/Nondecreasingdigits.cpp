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
	string s;
	
	while(cin>>s){
		string maxi=string((int)s.size(),'0');
		int last=0;
		string x="";
		bool ok=1;
		
		for(int i=0;i<(int)s.size();i++){
			int val=s[i]-'0';
			if(val>last){
				string w=x;
				w+='0'+(val-1);
				w+=string((int)s.size()-(int)w.size(),'9');
				maxi=max(maxi,w);
				x+='0'+(val);
				last=val;
			}else 
			if(val==last){
				x+='0'+val;
			}else if(val<last){
				ok=0;
				break;
			}
			
		}
		
		if(ok)maxi=max(maxi,x);
		
		
		long long ans=0;
		for(int i=0;i<(int)maxi.size();i++)
			ans=ans*10+(maxi[i]-'0');
		
		cout<<ans<<endl;
	}	
	
	// 448
	
	return 0;
}
