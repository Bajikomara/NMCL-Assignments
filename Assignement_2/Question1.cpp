#include<iostream>
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
	int n=100,i,j,ul=0,ur=0,um=1;
	double dx=4.0/n,CFL=0.8,dt=CFL*dx,T=1.33,df1,df2,sigma1,sigma2,a,b,df1_fvm,df2_fvm,a_fvm,b_fvm;
	double a1,a2,ac1,ac2,s=0.5;
	int nt=int(T/dt)+1;	
	double u[n][nt],x[n],u_fvm[n][nt],u_exact[n];
	FILE * f;
	x[0]=-2;
	for(i=1;i<=n;i++){
		x[i]=x[0]+i*dx;
		//cout<<x[i]<<endl;
	}
	for(i=1;i<n;i++){
		if(x[i]>(-0.33) && x[i]<(0.33)){
			u[i][0]=1;
		}if (x[i]<(-0.33) && x[i]>(0.33)){
			u[i][0]=0;
		}
		//cout<<x[i]<<","<<u[i][0]<<endl;
	}
	for(j=1;j<nt;j++){
		for(i=1;i<n;i++){
			u[1][j-1] =0;
			u[n-1][j-1] =0;
			df1=((u[i][j-1]*u[i][j-1])-(u[i-1][j-1]*u[i-1][j-1]))*0.5;
			df2=((u[i+1][j-1]*u[i+1][j-1])-(u[i][j-1]*u[i][j-1]))*0.5;
			if((u[i][j-1]-u[i-1][j-1])==0){
				sigma1=0;
			}else{
				sigma1=sign(df1/(u[i][j-1]-u[i-1][j-1]));
			}
			if((u[i+1][j-1]-u[i][j-1])==0){
				sigma2=0;
			}else{
				sigma2=sign(df2/(u[i+1][j-1]-u[i][j-1]));
			}
			a=((1+sigma1)/2)*df1;
			b=((1-sigma2)/2)*df2;
			u[i][j]=u[i][j-1]-(dt/dx)*(a+b);	
		}
		
	}
	
	//----------------FVM---------------------------------------------
	for(i=1;i<n;i++){
		if(x[i]>(-0.33) && x[i]<(0.33)){
			u_fvm[i][0]=1;
		}if (x[i]<(-0.33) && x[i]>(0.33)){
			u_fvm[i][0]=0;
		}
		//cout<<x[i]<<","<<u[i][0]<<endl;
	}
	for(j=1;j<nt;j++){
		for(i=1;i<n;i++){
			u_fvm[1][j-1] =0;
			u_fvm[n-1][j-1] =0;
			df1_fvm=((u_fvm[i][j-1]*u_fvm[i][j-1])-(u_fvm[i-1][j-1]*u_fvm[i-1][j-1]))*0.5;
			df2_fvm=((u_fvm[i+1][j-1]*u_fvm[i+1][j-1])-(u_fvm[i][j-1]*u_fvm[i][j-1]))*0.5;
			if((u_fvm[i][j-1]-u_fvm[i-1][j-1])==0){
				a1=0;
			}else{
				a1=(df1_fvm/(u_fvm[i][j-1]-u_fvm[i-1][j-1]));
			}
			if((u_fvm[i+1][j-1]-u_fvm[i][j-1])==0){
				a2=0;
			}else{
				a2=(df2_fvm/(u_fvm[i+1][j-1]-u_fvm[i][j-1]));
			}
			
			ac1=std::min(a1, 0.0);
			ac2=std::min(a2, 0.0);
			//cout<<a1<<","<<ac1<<endl;
			a_fvm=((u_fvm[i][j-1]*u_fvm[i][j-1])*0.5)+(ac1*(u_fvm[i+1][j-1]-u_fvm[i][j-1]));
			b_fvm=((u_fvm[i-1][j-1]*u_fvm[i-1][j-1])*0.5)+(ac2*(u_fvm[i][j-1]-u_fvm[i-1][j-1]));
			u_fvm[i][j]=u_fvm[i][j-1]-(dt/dx)*(a_fvm-b_fvm);	
		}
	}
	for(i=1;i<n;i++){
		if(x[i] < ((ul*T)-(0.33))){
			//u[i]=x[i]+(0.33/(um*t));
			u_exact[i]=0;
		}if(x[i] >= ((ul*T)-(0.33)) && x[i] <= ((um*T)-(0.33))){
			u_exact[i]=(x[i]+(0.33))/(um*T);
		}if(x[i] > ((um*T)-(0.33)) && x[i] <= ((s*T)+(0.33))){
			u_exact[i]=1;
		}else if((x[i] > ((s*T)+(0.33)))){
			u_exact[i]=0;
		}
			
	}
	f=fopen("Question 1 .csv","w");
	fprintf(f,"Location,Initial,FDS,FVM,Exact\n");
	for(int i=1;i<n;i++){
		cout<<x[i]<<","<<u[i][0]<<","<<u[i][nt-1]<<","<<u_fvm[i][nt-1]<<","<<u_exact[i]<<endl;
		fprintf(f," %lf,%lf,%lf,%lf,%lf\n",x[i],u[i][0],u[i][nt-1],u_fvm[i][nt-1],u_exact[i]);
	}
}

int main(){
	FDS();
	return 0;
}
