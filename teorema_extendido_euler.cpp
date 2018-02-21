#include <iostream>
#include <stdio.h>
#include <algorithm>
using namespace std;
long long n, m, p;
bool flag=0;
long long phi(long long N){
    long long dev=1;
	for(long long i=2;i*i<=N;i++){
		int a=0;
		while(N%i==0){
			N/=i;
			a++;
		}
		if(a!=0){
            for(int j=0;j<a-1;j++)
                dev*=i;
            dev*=i-1;
        }
    }
	if(N>1)dev*=N-1;	
	return dev;
}
long long power( long long x, long long y, long long z ){
	if ( x == 0 && y == 0 ) return 1;
	else if ( x == 0 ) return 0;
	else if ( x == 1 ) return 1;
	long long ans = 1, t = x;
	if ( x >= z ) flag = 1;
	t = x%z;
	while ( y > 0 ){
		if ( y & 1 ) ans = ans*t;
		if ( ans >= z ){
			flag = 1;
			ans %= z;
		}
		t *= t;
		if ( t >= z ){
			flag = 1;
			t %= z;
		}
		y /= 2;
	}
	return ans;
}
long long cal_f( long long a, long long b ){
	long long ans;
	if ( a == 0 ){
		flag = 0;
		return 1;
	}
	long long x = a%10, y = a/10, tmp, p;
	p=phi(b);
	tmp = cal_f(y,p);
	if (flag)tmp+=p;
	flag = 0;
	ans = power( x, tmp, b );
	//printf( "%lld\n", ans );
	return ans;
}
int main(){
    int ca;
	long long dev;
	cin>>ca;
    while(ca--){
		cin>>n>>m;
		flag=0;
        dev=cal_f(n,m);
		cout<<dev<<endl;
    }
}

