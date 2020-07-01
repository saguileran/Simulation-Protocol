<!-- "## Simulation Protocol ## -->

Lattice Boltzmann Method (LBM) is a Computational Fluid Dynamis (CFD) originated from Automata Celular Method (LGA) and developed in 1990's by Hardy–Pomeau–de Pazzis and Frisch–Hasslacher–Pomeau.
Although this method is not recent, it has become more popular in the last decade between scientists and engineering due to its versatility to parallelize and model several phenomena. Now exits several LBM software (public and private) that can be used to make CFD, some of them are: Palabos, ProLB, pybm, OpenLB, etc.


# Description

This protocol describe how to implement and how works a LBM in two and three dimensions with absorbing boundaries and simple geometries. The idea is explore how LBM simulate acoustical waves in rectangular geometries and analyse the data generetad, using Fourier transform. This is just a component of the project, to visit the entire project check our main page on GitHub, [Acoustical Instruments.](https://github.com/saguileran/Acoustics-Instruments)

# Requirements

This code is made in c++ and python language, to run it you have to install g++ compiler; on ubuntu just execute # sudo apt install g++ # in a terminal, on windows you need install VisualC++. Also the OpenMP library is needed but in the majority of cases it come with the compiler.

To visualaize data three options are evaluted: GNUPLOT, ParaView and OpenCL. The first two just need its corresponding software, while OpenCL need FreeGlut library.
We recommend use linux system to execute the code or at lest a virutal machine.
#

# Lattice Boltzmann Method

The simplest code consits in three parts:

* Macroscopic quantities:



* Colide

* Stream

![](https://ars.els-cdn.com/content/image/1-s2.0-S0965997816301855-gr3.jpg)


Calling necessary libraries

```c++
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "omp.h"
```

```c++
const int proportion = 1;
const int Lx = 501*proportion, Ly = 50*proportion;
const int LFx = 330*(proportion), LFy = 14*(proportion);

const double ke = 0, kF = 1; //k_e =  k enviorment = 0 (absortion wall)
//const double Aperture_x = 2*proportion;
//const double Hole_pos = LFx/3;

const int Q = 5;
const double W0 = 1.0 / 3;

const double C = 0.5; // C<0.707 celdas/click
const double TresC2 = 3 * C * C;
const double AUX0 = 1 - TresC2 * (1 - W0);

const double tau = 0.5;
const double Utau = 1.0 / tau;
const double UmUtau = 1 - Utau;

```


Class
```c++
class LatticeBoltzmann
{
private:
  double w[Q];
  int V[2][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; // f[ix][iy][i]
public:
  LatticeBoltzmann(void);
  double rho(int ix, int iy, bool UseNew);
  double  Jx(int ix, int iy, bool UseNew);
  double  Jy(int ix, int iy, bool UseNew);
  double  Jz(int ix, int iy, bool UseNew);
  double feq(double rho0, double Jx0, double Jy0, int i);
  void Colide(void);
  void Stream(void);
  void Initialize(double rho0, double Jx0, double Jy0);
  void ImposeField(int t);
  void PrintGrid(const char * NombreArchivo, int t);
  void Print(int t, int ix, int iy, const char * NombreArchivo);
  void Microphone(int t, int ix, int iy, const char * NombreArchivo);
};

```


Velocity and weights vectors
```c++
LatticeBoltzmann::LatticeBoltzmann(void){
  //Cargar los pesos
  w[0] = W0;   w[1] = w[2] = w[3] = w[4] = (1 - W0) / 4.0;
  //Cargar los vectores
  V[0][0] = 0; V[0][1] = 1;  V[0][2] = 0;  V[0][3] = -1;  V[0][4]=  0;
  V[1][0] = 0; V[1][1] = 0;  V[1][2] = 1;  V[1][3] =  0;  V[1][4]= -1;
}

```

Macroscopic quantities

```c++
double LatticeBoltzmann::rho(int ix, int iy, bool UseNew){
  int i; double suma = 0;
    for(i=0;i<Q;i++){if(UseNew) suma += fnew[ix][iy][i]; else suma += f[ix][iy][i]; }
  return suma;
}

double LatticeBoltzmann::Jx(int ix, int iy, bool UseNew){
  int i; double suma;
    for(suma=0,i=0;i<Q;i++){if(UseNew) suma += fnew[ix][iy][i]*V[0][i]; else suma += f[ix][iy][i]*V[0][i]; }
    return suma;
}

double LatticeBoltzmann::Jy(int ix, int iy, bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++){if(UseNew) suma += fnew[ix][iy][i]*V[1][i]; else suma += f[ix][iy][i]*V[1][i];}
  return suma;
}


double LatticeBoltzmann::feq(double rho0, double Jx0, double Jy0, int i){
  if(i==0){ return rho0 * AUX0;}
  else{     return w[i] * (TresC2 * rho0 + 3* (V[0][i] *Jx0 + V[1][i] * Jy0));}
}

```


Colide operator
```c++
void LatticeBoltzmann::Colide(void){
  int ix, iy, iz, i; double rho0, Jx0, Jy0;  //for all cell

#pragma omp paralel for
  {
  for(iy=0;iy<Ly;iy++){
    for(ix=0;ix<Lx;ix++){
        //Calcular las cantidades macroscópicas
	rho0 = rho(ix, iy, false);  Jx0 = Jx(ix, iy, false);  Jy0 = Jy(ix, iy, false);
	
	fnew[ix][iy][0] = UmUtau*f[ix][iy][0] + Utau*feq(rho0, Jx0, Jy0, 0);
	//for(i=0; i<Q; i++){ fnew[ix][iy][i] = UmUtau*f[ix][iy][i] + Utau*feq(rho0, Jx0, Jy0, i);}

	//Left flute wall
	if(ix == 20 &&  iy >= Ly/2 - LFy/2 && iy <= Ly/2 + LFy/2 ) {fnew[ix][iy][1] = kF * f[ix][iy][2]; fnew[ix][iy][2] = kF * f[ix][iy][1];}
	else if(ix == Lx - 1 || ix == 1){ fnew[ix][iy][1] = ke * fnew[ix][iy][3]; fnew[ix][iy][3] = ke * fnew[ix][iy][1]; } 
	else{ fnew[ix][iy][1] = UmUtau*f[ix][iy][1] + Utau*feq(rho0, Jx0, Jy0, 1);
    	      fnew[ix][iy][3] = UmUtau*f[ix][iy][3] + Utau*feq(rho0, Jx0, Jy0, 3); }

	//hirizontall fulte walls, the comment lines are the holes
	if(((iy == Ly/2 - LFy/2 && ix >= 20 && ix <= 20 + LFx) || (iy == Ly/2 + LFy/2  &&  ix >= 20 && ix <= 20 + LFx)
	    //&&  not(ix >= 20 + Hole_pos + LFx/3 - Aperture_x/2 && ix <= 20 + Hole_pos + LFx/3 + Aperture_x/2 && iy == Ly/2 + LFy/2)
	    //&&  not(ix >= 20 + Hole_pos - Aperture_x/2 && ix <= 20 + Hole_pos + Aperture_x/2 && iy == Ly/2 + LFy/2)
	    )
	   )
	  {fnew[ix][iy][4] =  kF * fnew[ix][iy][2]; fnew[ix][iy][4] =  kF * fnew[ix][iy][4];}
	else if(iy == Ly - 2 || iy == 1){ fnew[ix][iy][2] = ke *  f[ix][iy][4]; fnew[ix][iy][4] = ke *  f[ix][iy][2]; } 
	else{ fnew[ix][iy][2] = UmUtau*f[ix][iy][2] + Utau*feq(rho0, Jx0, Jy0, 2);
	      fnew[ix][iy][4] = UmUtau*f[ix][iy][4] + Utau*feq(rho0, Jx0, Jy0, 4); }
	
	//std::cout << ix << " " << iy << " " << rho0 << std::endl; //microphone place
	//if(ix == Lx-1 and iy == Ly-1){std::cout << " " << std::endl;}
	
    }
  }
 }
}

```


Stream

```c++
void LatticeBoltzmann::Stream(void){
  #pragma omp paralel for
  {
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)	
	//This condition disable lattice periodicity
	if( ix + V[0][i] < Lx && ix + V[0][i] >= 0 && iy + V[1][i] < Ly && iy + V[1][i] >= 0){ 
	  f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i] = fnew[ix][iy][i];}
  }
}

```


Initialize
```c++
void LatticeBoltzmann::Initialize(double rho0,double Jx0,double Jy0){
  #pragma omp paralel for
  {
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
	f[ix][iy][i] = feq(rho0, Jx0, Jy0, i);
}
}

```

Source

```c++
void LatticeBoltzmann::ImposeField(int t){
  int i, ix, iy; double lambda, omega, rho0, Jx0, Jy0;
  
  //sin(omega * t), declare initial function variables
  lambda = 10; omega = 2 * M_PI / lambda; ix = 22 ; iy = (Ly) / 2;

  //Initialize macroscopic cuantities
  rho0 = 10 * sin(omega*t); //Source function
  Jx0 = Jx(ix, iy, false);  Jy0 = Jy(ix, iy, false);

  //make a pulse of 500 steps
  //  if(t < 2000){
    //for(iy=(Ly-LFy+2)/2;iy<(Ly+LFy)/2;iy++){
  for(i=0;i<Q;i++)
    fnew[ix][iy][i] = feq(rho0, Jx0, Jy0, i);
    //    }
    //}
}

```
Data exportaton
```c++
void LatticeBoltzmann::PrintGrid(const char * NombreArchivo, int t){
  std::ofstream MiArchivo(NombreArchivo + std::to_string(t));
  double rho0, Jx0, Jy0;
  #pragma omp paralel for
  {
  MiArchivo << "X,Y,rho\n";
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
        rho0 = rho(ix, iy, false);  //Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
	//export data to paraview visualization
	MiArchivo << ix << "," << iy << "," << rho0 << '\n';
    }
  }
  }
  MiArchivo.close();
}

void LatticeBoltzmann::Print(int t, int ix, int iy, const char * NombreArchivo){
  double rho0 = rho(ix, iy, false);
  std::ofstream ofs;
  ofs.open(NombreArchivo, std::ofstream::out | std::ofstream::app);
  ofs << t << ',' << rho0 << '\n';
  ofs.close();
}

void LatticeBoltzmann::Microphone(int t, int ix, int iy, const char * NombreArchivo){
  double suma = 0; 

  for(int i = 1; i < 4; i++){
    double rho0 = rho(ix + i, iy - 1, false);
    double rho1 = rho(ix + i, iy,     false);
    double rho2 = rho(ix + i, iy + 1, false);

    suma += rho0 + rho1 + rho2;
    }
  suma = suma / 9;
  std::ofstream ofs;
  ofs.open(NombreArchivo, std::ofstream::out | std::ofstream::app);
  ofs << t << '\t' << suma / 9 << '\n';
  ofs.close();
}
```

Main Function
```c++
int main(void){
  
  LatticeBoltzmann Ondas;  
  int t,tmax=1000;
  
  //GNUPLOT LINES
  
  Ondas.Initialize(0,0,0);
  for(t=0;t<tmax;t++){

    //print actual step
    std::cout << t << std::endl;

    //Making lattice steps
    Ondas.Colide();
    Ondas.ImposeField(t);
    Ondas.Stream();

    //Export microphoes data, time vs pressure
    Ondas.Print(t, 20+LFx,    Ly/2, "LongPulse10k-0mm.dat");

    //Commands to make data animation
    if(t%5 == 0){Ondas.PrintGrid("LongPulse10k.csv.", t);}

    //Uncomment to use GNUPLOT animation
    //std::cout << "splot 'Ondas.dat' using 1:2:3  with points palette pointsize 3 pointtype 7 " << std::endl;
  }
  return 0;
}

```

GNUPlot commands
```c++
  std::cout << "set terminal gif animate" << std::endl;
  std::cout << "set output 'pelicula0.gif'" << std::endl;
 
  std::cout << "set pm3d map; set palette color positive" << std::endl;
  std::cout << "set palette defined (-1 \"red\", 0 \"white\", 1 \"blue\")" << std::endl;
  std::cout << "set cbrange[-1:1]" << std::endl;
  std::cout << "set xrange[-1:501]; set yrange[-1:51]; set zrange[-1:1]" << std::endl;
  std::cout << "set view map scale 1 " << std::endl;

  std::cout << "set view map;  set size ratio .9 " << std::endl;
  std::cout << "set object 1 rect from graph 0, graph 0 to graph 1, graph 1 back " << std::endl;
  std::cout << "set object 1 rect fc rgb 'black' fillstyle solid 1.0 " << std::endl;
  ```
