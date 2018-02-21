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
using namespace std;
int main()
{
	    int n;cin>>n;	
        vector<int>dev(10,0);
		for(int i=1;i<=n;i*=10)
		{
			int a=(n/i)/10;
			for(int j=0;j<10;j++)dev[j]+=a*i;
			dev[0]-=i;
			for(int j=0;j<(n/i)%10;j++)dev[j]+=i;
			dev[(n/i)%10]+= (n%i)+1;
		}
		for(int i=0;i<=9;i++)
		cout<<i<<" "<<dev[i]<<endl;  
		
    system("pause");
   return 0;
}
