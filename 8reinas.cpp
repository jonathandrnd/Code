#include<iostream>
#include<string>
#include<vector>
#include<sstream>
#include<cmath>
#include<algorithm>
#include<numeric>
#include<map>
#include<set>
#include<queue>
using namespace std ;
int c[15][15];
int columna[15];
int diagonal[15][15];
int cont=0;
int n=8;
void imprimir()
{
    cont++;
    cout<<"escenario "<<cont<<endl;
    for(int i=0;i<n;i++)
    {   for(int j=0;j<n;j++)
            cout<<c[i][j];    
        cout<<endl;
    }
    cout<<endl;
}
void reina(int x)
{
    if(x==n)imprimir();
    for(int i=0;i<n;i++)
    {
        if(diagonal[x][i]==0 && columna[i]==0)
        {
            columna[i]=1;
            for(int j=0;j<min(n-x,n-i);j++)diagonal[x+j][i+j]++;
            for(int j=0;j<min(n-x,i+1);j++)diagonal[x+j][i-j]++;
            c[x][i]=1;
            reina(x+1);
            columna[i]=0;
            for(int j=0;j<min(n-x,n-i);j++)diagonal[x+j][i+j]--;
            for(int j=0;j<min(n-x,i+1);j++)diagonal[x+j][i-j]--;
            c[x][i]=0;
        }
    }
        
}
int main()
{
    memset(c,0,sizeof(c));
    memset(columna,0,sizeof(columna));
    memset(diagonal,0,sizeof(diagonal));
    reina(0);
    
   system("pause");
    return 0;
}


