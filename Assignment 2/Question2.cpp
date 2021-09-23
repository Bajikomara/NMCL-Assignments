#include<bits/stdc++.h>
using namespace std;
int main(){
    FILE * f;
    int R_l=1,L_l=0;
    
    double T_stp;
    cout<<"Enter the Time Step:";
    cin>>T_stp;
    cout<<endl;
    double dx=0.005,dt=0.004,x_o=0.4,x_bar=0.075,p_bar=0.2;
    int J= int((T_stp/dt)+1);
	// cout<<J<<endl;
    double L_L,R_L;
	L_L=L_l-(0.5*dx*J),R_L=R_l+(0.5*dx*J);
    
    int I=int ((-L_L+R_L)/dx);
    // cout<<L_L<<" "<<R_L<<" "<<I<<endl;
    
    double X[I];
    double U[I][J],P[I][J],alpha1,alpha2,K_x[I],rho_x[I],c_o[I],lam_p_p,lam_p_n,Z_o[I],df1,df2,dnm1,dnm2;
    double U2[I][J],P2[I][J],nu[I],alpha1p,alpha1n,alpha2n,alpha2p,a1,a2,a3,a4,K_x2[I],rho_x2[I],c_o2[I],lam_p2,Z_o2[I],dnm1_2,dnm2_2;
    X[0]=L_L;
    for(int i=1;i<I;i++){
        X[i]=X[i-1]+dx;
    }
    //First Order======================================================================================
    //Initial Conditions 
    for(int i=0;i<I;i++){
        U[i][0]=0.0;
        
        if (abs(X[i]-x_o)<x_bar){
            P[i][0]=p_bar*(sqrt(1-(pow(((X[i]-x_o)/x_bar),2))));
        }else{
            P[i][0]=0.0;
        }
    }
    // Constants
    for (int j=0;j<J;j++){
        for(int i=0;i<I;i++ ){
            K_x[i]=1; 
            if(X[i]<0.6)
            {
                rho_x[i]=1.0;
            }else if(X[i]>0.6)
            {
                rho_x[i]=4.0;
            }
            c_o[i]=sqrt(K_x[i]/rho_x[i]);
            
        }
    }
   
	for (int j=1;j<J;j++) // time loop
    {   
        for(int i=1;i<(I-1);i++)  // grid loop
		{
            Z_o[i]=rho_x[i]*c_o[i];
            dnm1=Z_o[i]+Z_o[i+1];// denominator for alpha1
            dnm2=Z_o[i]+Z_o[i-1];// denominator for alpha2
        	alpha1=( -( P[i+1][j-1] - P[i][j-1] ) + Z_o[i+1]*( U[i+1][j-1] - U[i][j-1] ) )/dnm1;
            alpha2=( ( P[i][j-1] - P[i-1][j-1] ) + Z_o[i-1]*( U[i][j-1] - U[i-1][j-1] ) )/dnm2;
            
            //Eigen Values============================
            lam_p_p=sqrt(K_x[i]/rho_x[i]); 
            lam_p_n=-1*lam_p_p;
            //========================================
            df1=lam_p_p * alpha1 * Z_o[i];//df(i-1/2)
            df2=lam_p_n * alpha2 * (-1*Z_o[i]);//df(i+1/2)
            //Pressure and velocity Calculations at each grid point and time step
            P[i][j]=P[i][j-1]-(dt/dx)*(df1+df2);
            U[i][j]=U[i][j-1]-(dt/dx)*( (lam_p_p * alpha2) +  (lam_p_n * alpha1)   );
        }        
    }
    
    // Second Order=====================================================================================
    // Initial COnditions
    for(int i=0;i<I;i++){
        U2[i][0]=0.0;
        
        if (abs(X[i]-x_o)<x_bar){
            P2[i][0]=p_bar*(sqrt(1-(pow(((X[i]-x_o)/x_bar),2))));
        }else{
            P2[i][0]=0.0;
        }
    }
    
   
    
    // Constants
    for(int i=0;i<I;i++ ){
        K_x2[i]=1;
        if(X[i]<0.6)
        {
            rho_x2[i]=1.0;
        }else if(X[i]>0.6)
        {
            rho_x2[i]=4.0;
        }
        c_o2[i]=sqrt(K_x2[i]/rho_x2[i]);
        nu[i]=c_o2[i]*(dt/dx);
        
    }
   
    
      	
    //============================================
	for (int j=1;j<J;j++) // time loop
    {   
        for(int i=1;i<(I-1);i++)  // grid loop
		{
            lam_p2=c_o2[i];
            Z_o2[i]=rho_x2[i]*c_o2[i];
            dnm1_2=Z_o2[i]+Z_o2[i+1];
            dnm2_2=Z_o2[i]+Z_o2[i-1];
            
            
        	alpha1n=( -( P2[i][j-1] - P2[i-1][j-1] ) + Z_o2[i]*( U2[i][j-1] - U2[i-1][j-1] ) )/dnm2_2;
            alpha1p=( -( P2[i+1][j-1] - P2[i][j-1] ) + Z_o2[i+1]*( U2[i+1][j-1] - U2[i][j-1] ) )/dnm1_2;

            alpha2n=( ( P2[i][j-1] - P2[i-1][j-1] ) + Z_o2[i-1]*( U2[i][j-1] - U2[i-1][j-1] ) )/dnm2_2;
            alpha2p=( ( P2[i+1][j-1] - P2[i][j-1] ) + Z_o2[i]*( U2[i+1][j-1] - U2[i][j-1] ) )/dnm1_2;

            a1=(0.5*(1-nu[i])*(-lam_p2)*alpha1n);
            a2=(0.5*(1+nu[i])*(-lam_p2)*alpha1p);
            a3=(0.5*(1+nu[i])*(lam_p2)*alpha2n);
            a4=(0.5*(1-nu[i])*(lam_p2)*alpha2p);

            P2[i][j]=P2[i][j-1]-(dt/dx)*(a1*(-Z_o2[i])+a2*(-Z_o2[i])+a3*(Z_o2[i])+a4*(Z_o2[i]));
            U2[i][j]=U2[i][j-1]-(dt/dx)*(a1+a2+a3+a4); 
        }        
    }
	
	
	//Output File
    // Writing Output To file 
    
    f=fopen("Question_1.txt","w");
    // fprintf(f,"Time Step =");
    fprintf(f,"%lf\n",T_stp);
    f=fopen("Question_1.csv","w");
	fprintf(f,"Location,init_P,P_1st,P_2nd,U_1st,U_2nd\n");
    // fprintf(f,"Location,init_P,U,P\n");
    for(int i=0;i<I;i++){
        fprintf(f," %lf,%lf,%lf,%lf,%lf,%lf\n",X[i],P[i][0],P[i][J-1],P2[i][J-1],U[i][J-1],U2[i][J-1]);
        // fprintf(f," %lf,%lf,%lf,%lf\n",X[i],P[i][0],U[i][J-1],P[i][J-1]);
        // fprintf(f," %lf,%lf,%lf,%lf\n",X[i],P2[i][0],U2[i][J-1],P2[i][J-1]);
        //cout<<i<<setw(15)<<X[i]<<setw(15)<<P[i][0]<<setw(15)<<U[i][J-1]<<setw(15)<<P[i][J-1]<<endl;
	}

    
    return 0;
}
