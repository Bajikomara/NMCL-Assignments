#include<iostream>
#include<algorithm>
#include<math.h>
using namespace std;
void FVM(){
	int n,i,j,q2,T;
	double lbr,rbr;
	cout<<"Enter the Number of Iterations"<<endl;
	cin>>n;
	cout<<"Enter the which part of 2nd question"<<endl;
	cin>>q2;
	cout<<"Enter the Time"<<endl;
	cin>>T;
	cout<<"Enter the Left Boundary"<<endl;
	cin>>lbr;
	cout<<"Enter the Right Boundary"<<endl;
	cin>>rbr;
	double dx=(rbr-lbr)/(n-1),CFL=0.8,dt=CFL*dx,df1_fvm,df2_fvm,sigma1,sigma2,a_fvm,b_fvm,ac1,ac2,a1,a2,lb,rb;
	int nt=int(T/dt)+1;	
	double rho_fvm[n][nt],x[n];
	FILE * f;
	x[0]=lbr;
	for(i=1;i<=n;i++){
		x[i]=x[0]+i*dx;
	}
	switch (q2){
		case 1:{
			
			for(i=0;i<n;i++){
				rho_fvm[i][0]=0.25+(0.75*exp(-0.25*x[i]*x[i]));
			}
				lb=0.25;
				rb=0.25;
		
			break;
		}
		case 2 :{
			for(i=0;i<n;i++){
				if(x[i]<0){
					rho_fvm[i][0]=0.25;	
					
				}if(x[i]>=0){
					rho_fvm[i][0]=1.0;
				}
			}
			lb=0.25;
			rb=1.0;
			break;
		}
		case 3:{
			for(i=0;i<n;i++){
				if(x[i]<0){
					rho_fvm[i][0]=1;
				}if(x[i]>=0){
					rho_fvm[i][0]=0;
				}
			}
			break;
		}
	}
	/*for(i=0;i<n;i++){
		if(x[i]<0){
			rho_fvm[i][0]=0.25;	
			
		}if(x[i]>=0){
			rho_fvm[i][0]=1.0;
		}
	}*/
	
		
	for(j=1;j<nt;j++){
		for(i=1;i<n;i++){
			rho_fvm[0][j-1]=lb;
			rho_fvm[n-1][j-1]=rb;	
			df1_fvm=(rho_fvm[i][j-1]*(1-rho_fvm[i][j-1]))-(rho_fvm[i-1][j-1]*(1-rho_fvm[i-1][j-1]))	;
			df2_fvm=(rho_fvm[i+1][j-1]*(1-rho_fvm[i+1][j-1]))-(rho_fvm[i][j-1]*(1-rho_fvm[i][j-1]))	;
			if((rho_fvm[i][j-1]-rho_fvm[i-1][j-1])==0){
				a1=0;
			}else {
				a1=df1_fvm/(rho_fvm[i][j-1]-rho_fvm[i-1][j-1]);
			}
			if((rho_fvm[i+1][j-1]-rho_fvm[i][j-1])==0){
				a2=0;
			}else{
				a2=df2_fvm/(rho_fvm[i+1][j-1]-rho_fvm[i][j-1]);
			}
			
			ac1=std::min(a1, 0.0);
			ac2=std::min(a2, 0.0);
			a_fvm=(rho_fvm[i][j-1]*(1-rho_fvm[i][j-1]))+(ac2*(rho_fvm[i+1][j-1]-rho_fvm[i][j-1]));
			b_fvm=(rho_fvm[i-1][j-1]*(1-rho_fvm[i-1][j-1]))+(ac1*(rho_fvm[i][j-1]-rho_fvm[i-1][j-1]));
			rho_fvm[i][j]=rho_fvm[i][j-1]-(dt/dx)*(a_fvm-b_fvm);
		}
	}
	for(i=1;i<(n-1);i++){
		cout<<x[i]<<","<<rho_fvm[i][0]<<","<<rho_fvm[i][nt-1]<<endl;
	}
}
int main(){
	FVM();
	return 0;
}
	
