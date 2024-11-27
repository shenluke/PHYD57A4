#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define N 1000
#define PI 3.14159265389793
#define rho 1.225
#define Uinf 
#define c1
#define d1 
#define c2 
#define d2 
#define c3 
#define d3
#define c4
void cross(double a[3], double b[3],double output[3]){
    output[0]=a[1]*b[2]-a[2]*b[1];
    output[1]=a[2]*b[0]-a[0]*b[2];
    output[2]=a[0]*b[1]-a[1]*b[0];
}
void dot(double a[3], double b[3], double output){
    output=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

void LU_solve(double A[N][N], double B[N], double X[N]){
    double L[N][N], U[N][N];
    for(j=0; j<n; j++)
    {
        for(i=0; i<n; i++)
        {
            if(i<=j)
            {
                U[i][j]=A[i][j];
                for(k=0; k<i-1; k++){
                    U[i][j]-=L[i][k]*U[k][j];}
                if(i==j){
                    L[i][j]=1;}
                else{
                    L[i][j]=0;}
            }
            else
            {
                L[i][j]=A[i][j];
                for(k=0; k<=j-1; k++){
                    L[i][j]-=L[i][k]*U[k][j];}
                L[i][j]/=U[j][j];
                U[i][j]=0;}
        }
    }
    for(i=0; i<n; i++)
    {
        Y[i]=B[i];
        for(j=0; j<i; j++)
        {
            Y[i]-=L[i][j]*Y[j];
        }
    }
    for(i=n-1; i>=0; i--)
    {
        X[i]= Y[i];
        for(j=i+1; j<n; j++)
        {
            X[i]-=U[i][j]*X[j];
        }
        X[i]/=U[i][i];
    }
}

void Lift_Force_Distribution(double Alpha[N], double Gamma[N], double chord[N], double dy_l, double dy_r, double L[N]){
//Use dF=CL(alpha(y))*dA=CL(alpha(y))*chord(y)*dy
    
}
void Tau(double L[N], double CM[3], double LL[N][3],double Tau){
    double Tau=0;
    double r[3];
    double tau[3];
    for (int i=0;i<N;i++){
        r[0]=LL[i][0]-CM[0];
        r[1]=LL[i][1]-CM[1];
        r[2]=LL[i][2]-CM[2];
        cross(r,L,tau);
        Tau+=tau[0];//x component
    }
}
//w_i induced velocity
void step(double Gamma, double dGamma[N], double Torque, double omega_x,double theta, double vz, double wi){
    //Calculate using 4th order symplectic method/Leapfrog
    //Calculate lift force
    //Calculate torque
    //Update velocity and angular velocity
    //Update position and angle theta
    
}
int main(){
    double CM[3];
    CM[0]=;
    CM[1]=0;
    CM[2]=;
    double BC[N][3];
    double LL[N][3];
    double chord[N],dGamma[N];
    double tip_frac=0.02;
    double cosb1,cosb2,cosbm,dp1,dp2,dpm;//m refer ti the middle point
    double vz=0;
    double az,Tau;
    double LL_slope=tan(40*PI/180)+0.75*(tan(63*PI/180)-tan(40*PI/180));
    double BC_slope=tan(40*PI/180)+0.25*(tan(63*PI/180)-tan(40*PI/180));//Slope refer to dx/dy
    double theta=0;
    double omega_x=0;
    double dx_l,dy_l,dz_l,dx_r,dy_r,dz_r,dz_t;
    double dy_l=2*(13.05)/N;
    double dy_r=2*(18.75)/N;
    double dz_l=-dy_l*3*PI/180;
    //Small Angle approximation
    double dz_r=dy_r*3*PI/180;
    double dz_t=-dy_r*1*PI/180;
    double dt=0.01;

    //Initialize LL BC for the undamaged wing
    for (int i=0;i<tip_frac*N;i++){
        BC[N-i-1][0]=(i+0.5)*dx_r;
        BC[N-i-1][1]=wing_span/2-(i+0.50)*dy_r;
        BC[N-i-1][2]=i*dz_t+0.5*dz_t;
        LL[N-i-1][0]=(i+1)*dx_r;
        LL[N-i-1][1]=wing_span/2-(i+1)*dy_r;
        LL[N-i-1][2]=i*dz_t;
    }
    for (int i=tip_frac*N;i<N/2;i++){
        BC[N-i-1][0]=(i+0.5)*dx_r;
        BC[N-i-1][1]=wing_span/2-(i+0.50)*dy_r;
        BC[N-i-1][2]=BC[N-tip_frac*N-2][2]+(i-tip_frac*N)*dz_r;
        LL[N-i-1][0]=(i+1)*dx_r;
        LL[N-i-1][1]=wing_span/2-(i+1)*dy_r;
        LL[N-i-1][2]=LL[N-tip_frac*N-2][2]+(i-tip_frac*N)*dz_r;
    }
    //Points for damaged wing
    for (int i=N/2;i<N;i++){
        BC[N-i-1][0]=BC[N/2][0]+(i-N/2+1)*dx_l;
        BC[N-i-1][1]=BC[N/2][1]+(i-N/2+1)*dy_l;
        BC[N-i-1][2]=BC[N/2][2]+(i-N/2+1)*dz_l;
        LL[N-i-1][0]=BC[N/2][0]+(i-N/2+1)*dx_l;
        LL[N-i-1][1]=BC[N/2][1]+(i-N/2+1)*dy_l;
        LL[N-i-1][2]=BC[N/2][2]+(i-N/2+1)*dz_l;
    }
    // Compute matrix A, relative to the plane this should not change
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){

        }
    }
    //loop through steps 
    for (int iter=0;iter<end;iter++){
        step();//4th order symplectic
    }
}
