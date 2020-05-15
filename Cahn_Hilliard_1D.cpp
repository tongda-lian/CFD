//
//  main.cpp
//  1D_Cahn_Hilliard_Nov13
//
//  Created by LIAN TONGDA on 2019/11/13.
//  Copyright © 2019 LIAN TONGDA. All rights reserved.
//

/* 1D Cahn-Hilliard Equation Numerical Simulation
 LIAN Tongda @ Tokyo Institute of Technology
 */

/* periodic condition   i = 2 to i = N-3
 c[0]=c[N-4];  c[1]=c[N-3];c[N-2]=c[1];c[N-1]=c[1]
 */

//Ver0.1 Nov20
//128.cpp
//1D_Cahn-Hilliard_Equation
//The relationship between Stability and Grid size(dx) or Time Step(dt)
//Assuming the characteristic length epsilon(e) is constant

//Ver0.2 Nov21
//RK.cpp
//Using function phi to compute spacial difference to make the form of code clear.
//Update RK3 and RK4 method

#include<iostream>
#include<cmath>
# include<cstdlib>
#include<cstring>
#include <string>
#include <stdlib.h>
#include<stdio.h>
#include<fstream>
# include<sstream>
#include<istream>
#include<iomanip>

using namespace std;


const double pi = 3.141592653589793626;
// N is number of points/mesh
const int N = 1600;
const double dx = (2*pi/N);

// e = 0.25*dx, 0.5*dx, dx, 2.0*dx, 4.0*dx etc. ;
// e is physically constant , to restrain the variable to single.
const double e = (2*pi/100);

// dt = 1.0e-5, 1.0e-4, 1.0e-3 etc.;
const double dt = 1.0e-7;

double c[N + 4],c1[N+4],c2[N+4],c3[N+4],cn[N+4];
double dc[N+4];
double error;
int i;
int n,m;

void init();
void Eu1();
void Error();
void output(int m);
void Errorout(int a,double b);
void RK2();
void RK3();
void RK4();
double phi(double dm2,double dm1, double d, double d1, double d2);


int main()
{
    using namespace std;
    init();
    for (n = 0;;n++)
    {
        
// Keep one scheme and annotate all others.
//(Euler Forward and Runge-Kutta 2nd, 3rd, 4th Order )
        
        Eu1();
//        RK2();
//        RK3();
//        RK4();
        
        
        if (n % 10000== 0)
        {
            //Error();
            printf("The%dth computation result\n", n);
            //printf("the max relative error of uv is:%f\n", error);
            cout << "the max relative error of uv is:";
            cout << setprecision(3) << error << endl;
            
            
            if (n >= 10000)
            {
                if (n % 100000 == 0)
                {
                    output(n);
                }
                if (error<1.0e-12)
                {
                    break;
                }
                //Check if the calculation is divergent
                // If error is nan, error == error becomes false
                if(!(error == error))
                {
                    cout<<"Divergent!!!"<<endl;
                break;
                }

            }
         
        }
        if (n % 10== 0)
        {
            Error();
            Errorout(n,error);
        }
    }
    return 0;
}

void init()
{
    //initial condition c(x) = (1+sin(x))/2
    for(i=2;i<N+2;i++)
    {
        c[i]=cos((i-2+0.5)*dx);
    }
    //periodic boundary condition
    c[0]=c[N];
    c[1]=c[N+1];
    c[N+2]=c[2];
    c[N+3]=c[3];
}

//Time Integration Euler First Order
void Eu1()
{
    for(i=2;i<N+2;i++)
    {
//        dc[i] = ((c[i - 1]* c[i - 1]* c[i - 1] - 2* c[i]* c[i]* c[i] + c[i + 1]* c[i + 1]* c[i + 1])/(dx*dx)
//            -(c[i-1] - 2*c[i] + c[i+1])/(dx*dx)
//            - e*e*(c[i-2] - 4*c[i-1] + 6*c[i] - 4*c[i+1] + c[i+2])/(dx*dx*dx*dx))*dt;

        dc[i] = phi(c[i-2],c[i-1],c[i],c[i+1],c[i+2])*dt;
//        cout << phi(c[i-2],c[i-1],c[i],c[i+1],c[i+2]) << endl;
    }
    for(i = 2; i < N+2;i++)
    {
        c[i] += dc[i];
    }
    //periodic boundary condition
    c[0]=c[N];
    c[1]=c[N+1];
    c[N+2]=c[2];
    c[N+3]=c[3];
}

//spacial discretization,  dm2 = c[i-2], dm1 = c[i-1], d = c[i], d1 = c[i+1], d2 = c[i+2]
double phi(double dm2,double dm1, double d, double d1, double d2)
{
    return (dm1*dm1*dm1 - 2*d*d*d + d1*d1*d1) / (dx*dx)
            -(dm1 - 2*d + d1) / (dx*dx)
            -e*e*(dm2 - 4*dm1 + 6*d - 4*d1 + d2) / (dx*dx*dx*dx);
}


//Runge-Kutta Second Order
void RK2()
{
    // (dc/dt)*dt|t=n
    for(i=2;i<N+2;i++)
    {
//        dc[i] = ((c[i - 1]* c[i - 1]* c[i - 1] - 2* c[i]* c[i]* c[i] + c[i + 1]* c[i + 1]* c[i + 1])/(dx*dx)
//            -(c[i-1] - 2*c[i] + c[i+1])/(dx*dx)
//            - e*e*(c[i-2] - 4*c[i-1] + 6*c[i] - 4*c[i+1] + c[i+2])/(dx*dx*dx*dx))*dt;
         dc[i] = phi(c[i-2],c[i-1],c[i],c[i+1],c[i+2])*dt;
    }
    //c(n+1)
    for(i = 2; i < N+2;i++)
    {
        c1[i] = c[i]+dc[i];
    }
    //t = n+1, periodic boundary condition
    c1[0]=c1[N];
    c1[1]=c1[N+1];
    c1[N+2]=c1[2];
    c1[N+3]=c1[3];
    
    // (dc/dt)*dt|t=n+1
    for(i=2;i<N+2;i++)
    {
//        dc[i] = ((c1[i - 1]* c1[i - 1]* c1[i - 1] - 2* c1[i]* c1[i]* c1[i] + c1[i + 1]* c1[i + 1]* c1[i + 1])/(dx*dx)
//            -(c1[i-1] - 2*c1[i] + c1[i+1])/(dx*dx)
//            - e*e*(c1[i-2] - 4*c1[i-1] + 6*c1[i] - 4*c1[i+1] + c1[i+2])/(dx*dx*dx*dx))*dt;
           dc[i] = phi(c1[i-2],c1[i-1],c1[i],c1[i+1],c1[i+2])*dt;
    }
    //c(n+2)
    for(i = 2; i < N+2;i++)
    {
        c1[i] += dc[i];
    }
    //t = n+2, periodic boundary condition
    c1[0]=c1[N];
    c1[1]=c1[N+1];
    c1[N+2]=c1[2];
    c1[N+3]=c1[3];
    
    //evolution
    for(i = 2; i < N+2;i++)
    {
        //dc needs to be modified
        dc[i] = (c1[i] - c[i])/2;
        c[i] = (c[i] + c1[i]) / 2 ;
    }
    //periodic boundary condition
    c[0]=c[N];
    c[1]=c[N+1];
    c[N+2]=c[2];
    c[N+3]=c[3];
}

//Runge-Kutta Third Order
void RK3()
{
    for(i = 2; i< N+2;i++)
    {
        cn[i] = phi(c[i-2],c[i-1],c[i],c[i+1],c[i+2]);
    }
    //periodic boundary condition
    cn[0]=cn[N];
    cn[1]=cn[N+1];
    cn[N+2]=cn[2];
    cn[N+3]=cn[3];
    
    for(i = 2; i< N+2;i++)
    {
        c1[i] = phi(c[i-2] + (2 * dt * cn[i-2] / 3),c[i-1]+ (2 * dt * cn[i-1] / 3),
                    c[i]+ (2 * dt * cn[i] / 3),c[i+1]+ (2 * dt * cn[i+1] / 3),
                    c[i+2]+ (2 * dt * cn[i+2] / 3));
    }
    c1[0]=c1[N];
    c1[1]=c1[N+1];
    c1[N+2]=c1[2];
    c1[N+3]=c1[3];
    
    for(i = 2; i< N+2;i++)
       {
           c2[i] = phi(c[i-2] + (2 * dt * c1[i-2] / 3),c[i-1]+ (2 * dt * c1[i-1] / 3),
                       c[i]+ (2 * dt * c1[i] / 3),c[i+1]+ (2 * dt * c1[i+1] / 3),
                       c[i+2]+ (2 * dt * c1[i+2] / 3));
       }
       c2[0]=c2[N];
       c2[1]=c2[N+1];
       c2[N+2]=c2[2];
       c2[N+3]=c2[3];
    
    for( i = 2 ; i < N + 2 ; i++)
    {
        dc[i] = dt * (2 * cn[i] + 3 * c1[i] + 3 * c2[i])/8;
    }
    
    for(i = 2; i < N+2;i++)
    {
        c[i] += dc[i];
    }
    //periodic boundary condition
    c[0]=c[N];
    c[1]=c[N+1];
    c[N+2]=c[2];
    c[N+3]=c[3];
}


//Runge-Kutta Fourth Order
void RK4()
{
    for(i = 2; i< N+2;i++)
    {
        cn[i] = phi(c[i-2],c[i-1],c[i],c[i+1],c[i+2]);
    }
    //periodic boundary condition
    cn[0]=cn[N];
    cn[1]=cn[N+1];
    cn[N+2]=cn[2];
    cn[N+3]=cn[3];
    
    for(i = 2; i< N+2;i++)
    {
        c1[i] = phi(c[i-2] + (dt * cn[i-2] / 2),c[i-1]+ (dt * cn[i-1] / 2),
                    c[i]+ (dt * cn[i] / 2),c[i+1]+ (dt * cn[i+1] / 2),
                    c[i+2]+ (dt * cn[i+2] / 2));
    }
    c1[0]=c1[N];
    c1[1]=c1[N+1];
    c1[N+2]=c1[2];
    c1[N+3]=c1[3];
    
    for(i = 2; i< N+2;i++)
       {
           c2[i] = phi(c[i-2] + (dt * c1[i-2] / 2),c[i-1]+ ( dt * c1[i-1] / 2),
                       c[i]+ (dt * c1[i] / 2),c[i+1]+ (dt * c1[i+1] / 2),
                       c[i+2]+ (dt * c1[i+2] / 2));
       }
       c2[0]=c2[N];
       c2[1]=c2[N+1];
       c2[N+2]=c2[2];
       c2[N+3]=c2[3];
    
    for(i = 2; i< N+2;i++)
    {
        c3[i] = phi(c[i-2] + (dt * c2[i-2]),c[i-1]+ ( dt * c2[i-1]),
                    c[i]+ (dt * c2[i]),c[i+1]+ (dt * c2[i+1]),
                    c[i+2]+ (dt * c2[i+2]));
    }
    c3[0]=c3[N];
    c3[1]=c3[N+1];
    c3[N+2]=c3[2];
    c3[N+3]=c3[3];
    
    
    for( i = 2 ; i < N + 2 ; i++)
    {
        dc[i] = dt * (cn[i] + 2 * c1[i] + 2 * c2[i] + c3[i])/6;
    }
    
    for(i = 2; i < N+2;i++)
    {
        c[i] += dc[i];
    }
    //periodic boundary condition
    c[0]=c[N];
    c[1]=c[N+1];
    c[N+2]=c[2];
    c[N+3]=c[3];
}



//Relative error
void Error()
{
    double tempt1, tempt2;
    tempt1 = 0;
    tempt2 = 0;
    for (i=2;i<N+2;i++)
    {
        tempt1 = tempt1 + dc[i]*dc[i];
        tempt2 = tempt2 + c[i]*c[i];
    }
    tempt1 = sqrt(tempt1);
    tempt2 = sqrt(tempt2);
    error = tempt1 / (tempt2 + 1e-30);
}

// Output for Tecplot360 in *.plt format
/*void output(int m)
{
    ostringstream name;
    name << "Cahn-Hilliard_" << m << ".plt";
    ofstream out(name.str().c_str());
    out<<"Title=\"Cahn-Hilliard\"\n"<<"VARIABLES=\"X\",\"C\"\n"<<"ZONE T=\"BOX\",I="<<N<<",F=POINT"<<endl;
    out<<"0 "<<(c[1]+c[2])/2<<endl;
    for(i=2;i<N+2;i++)
    {
        out << (i - 2+0.5)*dx << " " << c[i] << endl;
    }
    out<<2*pi<<" "<<(c[N+1]+c[N+2])/2<<endl;
    out.close();
}*/

//output for *.txt
void output(int m)
{
    ostringstream name;
    name << "CH" << m << ".txt";
    ofstream out(name.str().c_str());
    
  /*  out<<"Title=\"Cahn-Hilliard\"\n"<<"VARIABLES=\"X\",\"C\"\n"<<"ZONE T=\"BOX\",I="<<N<<",F=POINT"<<endl;
   */
    
    out<<"0 "<<(c[1]+c[2])/2<<endl;
    for(i=2;i<N+2;i++)
    {
        out << (i - 2+0.5)*dx << " " << c[i] << endl;
    }
    out<<2*pi<<" "<<(c[N+1]+c[N+2])/2<<endl;
    out.close();
}

//Output of Relative error
void Errorout(int a,double b)
{
    //ofstream eout("error.txt");
    ofstream eout;
    eout.open("error.txt",ios_base::app);
    eout<<a<<" "<<b<<endl;
    eout.close();
}
