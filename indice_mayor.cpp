#include <ctime>
#include <numeric>
#include <functional>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstring>
#include <climits>
#include <cmath>
#include <cctype>
#include <sstream>
#include <map>
#include <set>
#include <cstdio>
#include <queue>
#define all(v) (v).begin(),v.end()
using namespace std ;
long long toi(string s){long long t;istringstream is(s);is>>t;return t;}
string tos(long long t){stringstream st;st<<t;return st.str();}
long long lcm(long long a,long long b){return a*(b/__gcd(a,b));}
long long gcd(long long a, long long b){return __gcd(a,b);}
int main(){
   int n;
    while(scanf("%d",&n)==1){
        int dev[n];
        int alt[n];
        int st[n];
        for(int i=0;i<n;i++)scanf("%d",&alt[i]);
        int pos=0;
        memset(dev,-1,sizeof(dev));
        for(int i=0;i<n;i++){
            while(pos>0 && alt[i]>alt[st[pos-1]]){
                pos--;
                dev[st[pos]]=i;
            } 
            st[pos++]=i;
        }   
        for(int i=0;i<n;i++)
            printf("%d\n",dev[i]+1);
        
    }
    return 0;
}


