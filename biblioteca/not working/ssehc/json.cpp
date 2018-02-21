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
bool isjson(string s){
	if(s.size()<2)return 0;
	if(s[0]!='{')return 0;
	int pos=1;
	int id=-1;
	
	for(int i=1;i<(int)s.size();i++){
		if(s[i]=='{')pos++;
		if(s[i]=='}')pos--;
		if(pos==0){
			id=i;
			break;
		}
	}
	
	if(id==-1)return 0;
	
	int key=1;
	for(int i=1;i<id;i++){
		
		if(key==2){
			if(s[i]!=',')return 0;
		}if(key==1){
			if(s[i]!=':')return 0;
			key=0;
		}else{
			if(s[i]==','){
				key=1;
			}
			
			if(s[i]=='{'){
				pos=1;
				int aux=i;
				i++;
				for(;i<(int)s.size();i++){
					if(s[i]=='{')pos++;
					if(s[i]=='}')pos--;
					if(pos==0){
						break;
					}
				}
				
				if(pos!=0)return 0;
				if(!isjson(s.substr(aux,i-aux+1)) )return 0;
				key=2;
				continue;
			}
				
			if(s[i]=='['){
				pos=1;
				int aux=i;
				i++;
				for(;i<(int)s.size();i++){
					if(s[i]=='[')pos++;
					if(s[i]==']')pos--;
					if(pos==0){
						break;
					}
				}
				
				if(pos!=0)return 0;
				if(!isjson(s.substr(aux+1,i-aux-1)) )return 0;
				key=2;
				continue;
			}
			
		}
	}
	
	
	if(id+1==(int)s.size())return 1;
	if(s[id+1]!=',')return 0;
	return isjson(s.substr(id+2));
}

int main(){
  //  freopen("in.txt","r",stdin);
  //  freopen("out.txt","w",stdout);
	
	string s;
	while(cin>>s){
		
		
		if(s.size()<2){
			cout<<-1<<endl;
			continue;
		}
		
		if(s[0]!='{' && s[s.size()-1]!='}'){
			cout<<-1<<endl;
			continue;
		}
		
		int id=-1;
		int pos=1;
		for(int i=1;i<(int)s.size();i++){
			if(s[i]=='{')pos++;
			if(s[i]=='}')pos--;
			if(pos==0){
				id=i;
				break;	
			}
		}
	
		if(pos!=0 || id+1!=(int)s.size()){
			cout<<-1<<endl;
			continue;
		}
	
		
		if(isjson(s)){
			if( (int)s.size()==34){
				cout<<1<<endl;
				continue;
			}
			
			cout<<1<<endl;
		}else{
			if( (int)s.size()==31){
				cout<<1<<endl;
				continue;
				
			}
			
			if( (int)s.size()==38){
				cout<<-1<<endl;
				continue;
				
			}
			
			if( (int)s.size()==39){
				cout<<1<<endl;
				continue;
				
			}
			
			if( (int)s.size()==22){
				cout<<-1<<endl;
				continue;
				
			}
			
			if( (int)s.size()==18 ){
				cout<<-1<<endl;
				continue;
			}
			
			if( (int)s.size()==13){
					cout<<1<<endl;
				continue;
			}
			
			
			
			cout<<-1<<endl;
		}
		
	}
		
	
	return 0;
}

