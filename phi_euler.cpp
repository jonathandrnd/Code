#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
using namespace std;
int memo[3000001];
int sum[3000001]={0};

long long phi(long long x){
	long long dev=1;
	for(long long i=2;i*i<=x;i++){
		if(x%i==0){
			long long izq=1;
			while(x%i==0)
				izq*=i,x/=i;
			dev*=izq/i*(i-1);
		}
	}
	if(x>1)dev*=x-1;
	return dev;
}
int main(){
    memset(memo,-1,sizeof(memo));
    int a,b;
    sum[0]=0;
    sum[1]=0;
    sum[2]=0;
    for(int i=3;i<=3000000;i++)sum[i]=phi(i)+sum[i-1];
    while(scanf("%d %d",&a,&b)==2){
        printf("%d\n",sum[a]-sum[b-1]);   
    }
    return 0;
}
