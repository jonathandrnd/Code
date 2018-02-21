#include <stack>
#include <queue>
#include <vector>
#include <cstring>
#include <iostream>
#include <map>
#include <cstdio>

using namespace std;

#define MAX_V 100001

int n, num_scc, scc[MAX_V];
vector< vector<int> > G;
vector< vector<int> > GT;
bool visited[MAX_V];
stack<int> S;
queue<int> Q;

void dfs(int v){
    visited[v] = true;
    
    for(int i=G[v].size()-1;i>=0;--i)
        if(!visited[G[v][i]])
            dfs(G[v][i]);
    
    S.push(v);
}

void bfs(int v){
    Q.push(v);
    visited[v] = true;
    
    int aux;
    
    while(!Q.empty()){
        aux = Q.front();
        scc[aux] = num_scc;
        Q.pop();
        
        for(int i=GT[aux].size()-1;i>=0;i--){
            if(!visited[GT[aux][i]]){
                Q.push(GT[aux][i]);
                visited[GT[aux][i]] = true;
            }
        }
    }
}
void SCC(){
    memset(visited,false,sizeof(visited));
    
    for(int i=0;i<n;++i) if(!visited[i]) dfs(i);
    
    num_scc = 0;
    int aux;
    
    memset(visited,false,sizeof(visited));
    
    while(!S.empty()){
        aux = S.top();
        S.pop();
        
        if(!visited[aux]){
            bfs(aux);
            ++num_scc;
        }
    }
}

int main()
{
    //freopen("in.txt","r",stdin);
    //freopen("out.txt","w",stdout);
    
    int t,x;
    cin>>t;
    for(int i=0;i<t;i++){
        cin>>n>>x;
        G.clear();GT.clear();
        G.resize(n);
        GT.resize(n);
        for(int j=0;j<x;j++){
            int a,b;
            cin>>a>>b;
            G[a-1].push_back(b-1);
            GT[b-1].push_back(a-1);  
        }
        SCC();
        vector<bool>xuxa(n,true);
        for(int j=0;j<G.size();j++){
            for(int k=0;k<G[j].size();k++){
                if(scc[G[j][k]]!=scc[j]){
                    xuxa[scc[G[j][k]]]=false;           
                }    
            }
        }
        cout<<count(xuxa.begin(),xuxa.begin()+num_scc,true)<<endl;
    }
    
    //system("pause");
    return 0;
}


