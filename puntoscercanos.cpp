#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <queue>
#include <map>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <numeric>
using namespace std;
typedef double TipoPunto;
#define INF 1e50
struct point{
    TipoPunto x,y;
    int id;
}

X[10020],Y[10020];

bool lr[10005];
bool sortx(const point &a,const point &b){return a.x<b.x;}
bool sorty(const point &a,const point &b){return a.y<b.y;}

TipoPunto dist(point &a,point &b){
    TipoPunto dx=a.x-b.x;
    TipoPunto dy=a.y-b.y;
    return dx*dx+dy*dy;
}

TipoPunto solve(point X[],point Y[],int t){
    if(t==2)return dist(X[0],X[1]);
    if(t==3){
        TipoPunto res=INF;
        res<?=dist(X[0],X[1]);
        res<?=dist(X[0],X[2]);
        res<?=dist(X[1],X[2]);
        return res;
    }
    
    point Xl[t/2+1],Xr[t/2+1],Yl[t/2+1],Yr[t/2+1];
    memset(lr,0,sizeof(lr));
    
    int cl=0,cr=0;
    for(int i=0;i<t;i++)
        if(i<(t+1)/2)
            Xl[cl++]=X[i];
        else {
            Xr[cr++]=X[i];
            lr[X[i].id]=true;
        }

    cl=0; cr=0;
    for(int i=0;i<t;i++)
        if(!lr[Y[i].id])
            Yl[cl++]=Y[i];
        else
            Yr[cr++]=Y[i];
    
    TipoPunto dl=solve(Xl,Yl,cl);
    TipoPunto dr=solve(Xr,Yr,cr);
    TipoPunto d=min(dl,dr);
    
    point YP[t];
    double med=(X[t/2-1].x+X[t/2].x)/2;
    
    int ty=0;
    for(int i=0;i<t;i++){
        double v=fabs(Y[i].x-med);
        if(v*v<=d)
            YP[ty++]=Y[i];
    }
    
    for(int i=0;i<ty;i++)
        for(int j=i+1;j<ty&&j<i+7;j++)
            d<?=dist(YP[i],YP[j]);
    
    return d;
}

int main(){
    int n;
    cin>>n;

    for(int i=0;i<n;i++) {
        cin>>X[i].x>>X[i].y;
        Y[i].x=X[i].x; 
        Y[i].y=X[i].y;
        X[i].id=Y[i].id=i;
    }
    
    if(n==1){
        printf("INFINITY\n"); 
    }
    
    sort(X,X+n,sortx);
    sort(Y,Y+n,sorty);
    double dmin=sqrt(solve(X,Y,n));
    
    if(dmin<10000.0) 
        printf("%.4lf",dmin);
    else 
        printf("INFINITY");
    printf("\n");
    
    system("pause");
    return 0;
}
