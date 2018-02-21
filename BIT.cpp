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
int tree[200002];
int MaxVal=200002;
int read(int idx){
    if(idx==0)return 0;
	int sum = 0;
	while (idx > 0){
		sum += tree[idx];
		idx -= (idx & -idx);
	}
	return sum;
}

int get(int uno){
	int low=1;int hi=200001;
	for(int i=0;i<20;i++){
		int me=(low+hi)/2;
		if(read(me)>=uno)hi=me;
		else low=me+1;
	}
	return low;
}
void update(int idx ,int val){
	while (idx <= MaxVal){
		tree[idx] += val;
		idx += (idx & -idx);
	}
}
int main()
{
    freopen("in.txt","r",stdin);
    freopen("out.txt","w",stdout);
    int cases;
    //scanf("%d",&cases);
    cin>>cases;
    for(int casos=0;casos<cases;casos++){
        memset(tree,0,sizeof(tree));
        int n;
        cin>>n;
        //scanf("%d",&n);
        int c[n];
        for(int i=0;i<n;i++){
            int besito;cin>>besito;
            c[i]=besito;
            //scanf("%d",&c[i]);
        }
        for(int i=0;i<n;i++){
            //update(i+1,1);
        }
        int val=n;
        int dev[n];
        for(int i=n-1;i>=0;i--){
            //dev[i]=get(val-c[i]);
            val--;
        }
        for(int i=0;i<n;i++){
            if(i!=n-1)printf("%d ",dev[i]);    
            else printf("%d",dev[i]);
        }
        printf("%d\n");
    }
    //system("pause");
    return 0;
}


