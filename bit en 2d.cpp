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
int tree[1030][1030];
int max_x=20;int max_y=20;
int read(int idx,int idy){
    if(idx==0 ||idy==0)return 0;
	int sum = 0;
	int y1;
	while (idx >0){
		y1=idy;
        while(y1>0){
            sum += tree[idx][y1];
		    y1 -= (y1 & -y1);
        }
        idx -= (idx & -idx);
	}
	return sum;
}
void update(int x , int y , int val){
	int y1;
	while (x <= max_x){
		y1 = y;
		while (y1 <= max_y){
			tree[x][y1] += val;
			y1 += (y1 & -y1); 
		}
		x += (x & -x); 
	}
}
int main()
{
    int t;scanf("%d",&t);
    for(int caso=0;caso<t;caso++){
        int n;scanf("%d",&n);  
        n++;
        int h[n][n];
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)tree[i][j]=h[i][j]=0;
        max_x=n+1;max_y=n+1;
        char s[12];
        while(true){
            scanf("%s",&s);
            if(strcmp(s,"END")==0)break;
            int a,b,c,d;
            if(strcmp(s,"SET")==0){
                scanf("%d",&a);scanf("%d",&b);scanf("%d",&c);
                a++;b++;
                update(a,b,c-h[a][b]);
                h[a][b]=c;
            }
            if(strcmp(s,"SUM")==0){
                scanf("%d",&a);scanf("%d",&b);scanf("%d",&c);scanf("%d",&d);
                a++;b++;c++;d++;
                printf("%d",read(c,d)+read(a-1,b-1)-read(a-1,d)-read(c,b-1));
                printf("\n");
            }
        }
        printf("\n");
    }
    return 0;
}
