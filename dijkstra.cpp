#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
using namespace std;
struct edge{
    int x,y;
  edge(int X,int Y){
    x=X;
    y=Y;      
  }  
};
struct nodo{
    int x,y;
  nodo(int X,int Y){
    x=X;
    y=Y;      
  }  
};
bool operator<(nodo p,nodo q){
    return p.y>q.y;    
}
int d[5][5];
int n=5;int oo=1<<30;
vector<edge> v[5];
void dijkstra(int x){
    for(int i=0;i<n;++i)d[x][i]=oo;
    d[x][x]=0;
    priority_queue<nodo> Q;
    Q.push(nodo(x,0));
    while(!Q.empty()){
        nodo q=Q.top();
        Q.pop();
        for(int i=0;i<v[q.x].size();++i){
            edge aux=v[q.x][i];
            int temp=q.y+aux.y;
            if(temp<d[x][aux.x]){
                d[x][aux.x]=temp;
                Q.push(nodo(aux.x,temp));
            }
        }
    }
}
int main()
{
 //   freopen("in.txt", "r", stdin);
   // freopen("out.txt", "w", stdout);
    v[0].push_back(edge(1,1));
    v[0].push_back(edge(2,3));
    v[0].push_back(edge(3,7));
    v[1].push_back(edge(0,1));
    v[1].push_back(edge(4,10));
    v[2].push_back(edge(0,3));
    v[2].push_back(edge(4,5));
    v[3].push_back(edge(0,7));
    v[3].push_back(edge(4,8));
    v[4].push_back(edge(1,10));
    v[4].push_back(edge(2,5));
    v[4].push_back(edge(3,8));
    for(int i=0;i<n;i++)dijkstra(i);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
        cout<<d[i][j]<<" ";cout<<endl;}
    system("pause");
    return 0;
}
