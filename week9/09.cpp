#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "json.hpp"

using namespace std;
using json = nlohmann::json;

//эйлер
void euler(double phi[],double w[],double epsilon[],int N,double dt,double phi0, double omega, double gamma,double A,double W){
    phi[0]=phi0;
    w[0]=0;
    epsilon[0]=-omega*phi0-2*gamma*0;
    for(int i=0;i<N;i++){
        epsilon[i+1]=-phi[i]*omega-2*gamma*w[i]+A*sin(W*i*dt);
        w[i+1]=w[i]+epsilon[i]*dt;
        phi[i+1]=phi[i]+w[i]*dt;
    }
}
//хойн
void heun(double phi[],double w[],double epsilon[],int N,double dt,double phi0, double omega, double gamma,double A,double W){
    phi[0]=phi0;
    w[0]=0;
    epsilon[0]=-omega*phi0-2*gamma*0;
    for(int i=0;i<N;i++){
        epsilon[i+1]=-phi[i]*omega-2*gamma*w[i]+A*sin(W*i*dt);
        w[i+1]=w[i]+(epsilon[i]+epsilon[i+1])*dt/2;
        phi[i+1]=phi[i]+(w[i]+w[i+1])*dt/2;
    }
}
//обратный эйлер
void euler1(double phi[],double w[],double epsilon[],int N,double dt,double phi0, double omega, double gamma,double A,double W){
    phi[0]=phi0;
    w[0]=0;
    epsilon[0]=-omega*phi0-2*gamma*0;
    for(int i=0;i<N;i++){
        epsilon[i+1]=-phi[i]*omega-2*gamma*w[i]+A*sin(W*i*dt);
        w[i+1]=w[i]+epsilon[i+1]*dt;
        phi[i+1]=phi[i]+w[i+1]*dt;
    }
}
//рунге-кутте
void rk(double phi[],double w[],double epsilon[],int N,double dt,double phi0, double omega, double gamma,double A,double W){
    phi[0]=phi0;
    w[0]=0;
    for(int i=0;i<N;i++){
        double e1=-phi[i]*omega-2*gamma*w[i]+A*sin(W*i*dt);
        double wk1=w[i];
        double e2=-(phi[i]+wk1*dt/2)*omega-2*gamma*(w[i]+e1*dt/2)+A*sin(W*i*dt);
        double wk2=w[i]+e1*dt/2;
        double e3=-(phi[i]+wk2*dt/2)*omega-2*gamma*(w[i]+e2*dt/2)+A*sin(W*i*dt);
        double wk3=w[i]+e2*dt/2;
        double e4=-(phi[i]+wk3*dt)*omega-2*gamma*(w[i]+e3*dt)+A*sin(W*i*dt);
        double wk4=w[i]+e3*dt;
        w[i+1]=w[i]+(e1+2*e2+2*e3+e4)*dt/6;
        phi[i+1]=phi[i]+(wk1+2*wk2+2*wk3+wk4)*dt/6;
    }
}
int main(){
    json params;
    cin>>params;
    int N=params.value("N",10000);
    double dt=params.value("dt",0.01);
    double phi0=params.value("phi0",1.0);
    double omega=params.value("omega",1.0);
    double gamma=params.value("gamma",0.1);
    double A=params.value("A",1.0);
    double W=params.value("W",1);

    double tau[N+1];
    double E1[N+1];
    double E2[N+1];

    double phi1[N+1];  double w1[N+1];  double epsilon1[N+1];
    double phi2[N+1];  double w2[N+1];  double epsilon2[N+1];
    double phi3[N+1];  double w3[N+1];  double epsilon3[N+1];
    double phi4[N+1];  double w4[N+1];  double epsilon4[N+1];

    euler(phi1,w1,epsilon1,N,dt,phi0,omega,gamma,A,W);
    heun (phi2,w2,epsilon2,N,dt,phi0,omega,gamma,A,W);
    rk(phi3,w3,epsilon3,N,dt,phi0,omega,gamma,A,W);
    euler1(phi4,w4,epsilon4,N,dt,phi0,omega,gamma,A,W);

    tau[0]=0;
    E1[0]=phi0*phi0;
    E2[0]=phi0*phi0;
    for(int i=0;i<N;i++){
        tau[i+1]=tau[i]+dt;
        E1[i+1]=phi1[i]*phi1[i]+w1[i]*w1[i];
        E2[i+1]=phi2[i]*phi2[i]+w2[i]*w2[i];
    }

    ofstream f("09.dat");
    for(int i=0;i<=N;i++)
        f<<tau[i]<<" "
         <<phi1[i]<<" "<<phi2[i]<<" "<<phi3[i]<<" "
         <<w1[i]<<" "<<w2[i]<<" "<<w3[i]<<" "
         <<E1[i]<<" "<<E2[i]<<" "
         <<phi4[i]<<"\n";
    f.close();

    return 0;
}
