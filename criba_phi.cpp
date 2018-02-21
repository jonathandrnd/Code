#include <stdio.h>
#include <algorithm>
#include <vector>
using namespace std;
#define N 3000000
long long cnt[N] = { 0 };
bool mark[N] = { 0 };
vector<int> prime;
int main( ){
	int i, j;
	for ( i = 2; i < N; i++ ){
		if ( !mark[i] ){
			cnt[i] = i-1;
			prime.push_back( i );
		}
		for ( j=0;j<prime.size() && prime[j]*i < N; j++ ){
			mark[prime[j]*i]=1;
			if(i%prime[j]==0){
				cnt[i*prime[j]] = cnt[i]*prime[j];
				break;
			}
			else 
            cnt[i*prime[j]] = cnt[i]*(prime[j]-1 );
		}
	}
	for(i=3;i<N;i++)cnt[i]+=cnt[i-1];
	while(scanf("%d%d",&i,&j)!=EOF){
		printf("%lld\n",cnt[j]-cnt[i-1] );
	}
	return 0;
}
