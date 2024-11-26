#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define N 1000
#define PI 3.14159265389793
#define rho 1.225



void cross(double a[3], double b[3],double output[3]){
    output[0]=a[1]*b[2]-a[2]*b[1];
    output[1]=a[2]*b[0]-a[0]*b[2];
    output[2]=a[0]*b[1]-a[1]*b[0];
}
double find_max_abs(int n,double v[n]){
    double max=abs(v[0]);
    for (int i=0;i<n-1;i++){
        if (abs(v[i+1])>max){
            max=abs(v[i+1]);
        }
    }
    return max;
}
void LU_solve(double A[N][N], double b[N], double x[N]){
    double L[N][N], U[N][N];
    double interchanger;
}

void Lift_Force_Distribution(double Alpha[N], double Gamma[N], double chord[N], double dy_l, double dy_r, double L[N]){

}
void Tau(){

}
//w_i induced velocity
void step(double Gamma, double dGamma[N], double Torque, double omega_x,double theta, double vz, double Uinf, double wi){

}
int main(){
    double BC[N][3];
    double LL[N][3];
    double LL_slope=tan(40*PI/180)+0.75*(tan(63*PI/180)-tan(40*PI/180));
    double BC_slope=tan(40*PI/180)+0.25*(tan(63*PI/180)-tan(40*PI/180));//Slope refer to dx/dy
    double dx_l,dy_l,dz_l,dx_r,dy_r,dz_r,dz_t;
    double dy_l=2*(13.05-3.35/2)/N;
    double dy_r=2*(18.75-3.35/2)/N;
    double dz_l=-dy_l*3*PI/180;
    double dz_r=dy_r*3*PI/180;
    double dz_t=-dy_r*1*PI/180;
    //Initialize LL BC
    for (int i=0;i<N/2;i++){
        BC[i][0]=;
        BC[i][1]=;
        BC[i][2]=;
        LL[i][0]=;
        LL[i][1]=;
        LL[i][2]=;
    }
    for (int i=N/2;i<N;i++){
    }
    // Compute matrix A
    //loop through steps 
    
}