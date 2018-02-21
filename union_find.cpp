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
using namespace std ;
int parent[100001],rank[10001];
void create(int x){
    parent[x] = x; rank[x] = 0;
}
int find(int x){
    if(parent[x]!=x) parent[x]=find(parent[x]);
    return parent[x];
}
void merge(int x,int y){
    int PX = find(x),PY =find(y);
    if(rank[PX]>rank[PY]) parent[PY] = PX;
    else{
        parent[PX] = PY;
        if(rank[PX]==rank[PY]) ++rank[PY];
    }
}
int main(){
    int n;
    while(scanf("%d",&n)){
        if(n==0)break;
        for(int i=0;i<n;i++)create(i);
        for(int i=0;i<n-1;i++){
            int a,b;
            scanf("%d %d",&a,&b);    
            merge(a,b);
        }   
        for(int i=0;i<n;i++)cout<<i<<" "<<find(i)<<" "<<rank[i]<<endl;
    }
    system("pause");
    return 0;
}


