#include<iostream>
#include<math.h>
using namespace std;
int sign(double val){
	if(val >0){
		return 1;
	}if (val<0){
		return-1;
	}else{
		return 0;
	} 
	return 0;
}
void FDS(){
	int n,i,j,q2,T,nt;
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
	double dx=(rbr-lbr)/(n-1),CFL=0.8,dt=CFL*dx,df1,df2,sigma1,sigma2,a,b,lb,rb;
	nt= (T/dt)+1;
	double rho[n][nt],x[n];
	FILE * f;
	x[0]=lbr;
	for(i=1;i<=n;i++){
		x[i]=x[0]+i*dx;
	}
	switch (q2){
		case 1 :{
			for(i=0;i<n;i++){
				rho[i][0]=0.25+(0.75*exp(-0.25*x[i]*x[i]));
			}
			lb=0.25;
			rb=0.25;
			break;
		}
		case 2 :{
			for(i=0;i<n;i++){
				if(x[i]<0){
					rho[i][0]=0.25;
				}else if(x[i]>=0){
					rho[i][0]=1;
				}
			}
			lb=0.25;
			rb=1.0;
			break;
		}
		case 3 :{
			for(i=0;i<n;i++){
				if(x[i]<0){
					rho[i][0]=1.0;
				}if(x[i]>=0){
					rho[i][0]=0.0;
				}
			}
			lb=1.0;
			rb=0.0;
			break;
		}
	}
	for(j=1;j<nt;j++){
		rho[0][j-1] =lb; // BC
		rho[n-1][j-1] =rb; //BC
		for(i=0;i<n;i++){
			df1=((rho[i][j-1]*(1-rho[i][j-1]))-(rho[i-1][j-1]*(1-rho[i-1][j-1])));
			df2=((rho[i+1][j-1]*(1-rho[i+1][j-1]))-(rho[i][j-1]*(1-rho[i][j-1])));
			if((rho[i][j-1]-rho[i-1][j-1])==0){
				sigma1=0;
			}else{
				sigma1=sign(df1/(rho[i][j-1]-rho[i-1][j-1]));
			}
			if((rho[i+1][j-1]-rho[i][j-1])==0){
				sigma2=0;
			}else{
				sigma2=sign(df2/(rho[i+1][j-1]-rho[i][j-1]));
			}
			a=((1+sigma1)/2)*df1;
			b=((1-sigma2)/2)*df2;
			rho[i][j]=rho[i][j-1]-(dt/dx)*(a+b);	
		}
	}
	f=fopen("Question 2 FDs .csv","w");
	fprintf(f,"Location,Initial,FDS\n");
	for(int i=1;i<n;i++){
		cout<<x[i]<<","<<rho[i][0]<<","<<rho[i][nt-1]<<endl;
		fprintf(f," %lf,%lf,%lf\n",x[i],rho[i][0],rho[i][nt-1]);
	}
}

int main(){
	FDS();
	return 0;
}
