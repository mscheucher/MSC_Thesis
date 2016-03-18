//
//  KHI.cpp
//
//  Program: solving set of MHD equations in 2dim. using Total Variation Diminishing - Lax Friedrich scheme
//
//  Created by Markus Scheucher on 03/12/13.
//
//
//v1.0: 03.12.2013 - starting with 1dimension continuum equation: Stepfunction solved
//v1.1: 04.12.2013 - External Datainput - file: "input.txt";variable filenames
//v1.2: 04.12.2013 - Gauss function solver; end of 1dimensional
//v2.0: 05.12.2013 - first working 2D implementation
//v2.1: 06.12.2013 - computing time reduced by for-loop changes; input file structured - now: "input_2D"
//v2.2: 18.12.2013 - first running TVD Lax friedrich scheme implementation
//v3.0: 19.12.2013 - start with MHD equations
//v3.1: 21.01.2014 - first running code with implemented MHD equations
//v3.2: 07.02.2014 - divergence cleaning with field CD method implemented
//v3.3: 27.03.2014 - first test setup: 2 tanh boundaries implemented and input file expanded for parameter studies
//v3.4: 07.04.2014 - random disturbances, limiter selection and final setup for EGU - parameter study implemented

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <new>
#include <cstdlib>

using namespace std;

#include "input.h"

#define PI 3.14159265

//global variables ##################################################################

double ***u;				//conserved Variables U
double ***w;				//primitive Variables W
double ***du;				//limiter delta U quer
double ***up;				//U^(n+1/2)
double ***cmax;				//Cmax in X[0] and Y[1]
double ***fLn;				//flux Left at t=n for U^(n+1/2)
double ***fRn;				//flux Right at t=n for U^(n+1/2)
double ***fLp;				//flux Left at t=n+1/2 for Lx/Ly - Operator
double ***fRp;				//flux Right at t=n+1/2 for Lx/Ly - Operator
double ***f;				//f^(LR) n+1/2
double ***old;				//old B-Field for divergence cleaning

double dt;				//timer - to be calculated
double dx;				//delta in x
double dy;				//delta in y

double **randnum;			//Random number between 0. and 1.
char buffer[1024];          //buffer for filenames
char dir[1024];             //buffer for data directory

/*
int nx;					//# of gridpoints in x
int ny;					//# of gridpoints in y
double boxx;				//Length of Box
double boxy;				//Length of Box
int na;					//# of x-boundary layers in y direction
int nout;				//# of output files per parameter - output every tmax/nout
double tmax;				//Endtime of simulation
int x_boundaries,y_boundaries;		//boundary conditions: 0=transmissive, 1=periodic
int cleaning;				//divergence cleaning method: 0=no_cleaning, 1=fieldCD, 2=fluxCT, 3=fluxCD
int cmethod;				//method of calculating cmax: 0=after Barmin et.al, 1=max cmax_i over grid
int limiter;				//cell-jump limiter: 0=woodward, 1=minmod
double cfl_input;			//courant number after first time steps
double kappa;				//kappa
double mu;				//mu_0
double thick;				//boundary Layer thickness

double iR0;				//input parameter: rho - lowspeed
double iR1;				//input parameter: rho - highspeed
double iVx0;				//input parameter: vx - lowspeed
double iVx1;				//input parameter: vx - highspeed
double iVy0;				//input parameter: vy - lowspeed
double iVy1;				//input parameter: vy - highspeed
double iVz0;				//input parameter: vz - lowspeed
double iVz1;				//input parameter: vz - highspeed
double iP0;				//input parameter: p - lowspeed
double iP1;				//input parameter: p - highspeed
double iBx0;				//input parameter: bx - lowspeed
double iBx1;				//input parameter: bx - highspeed
double iBy0;				//input parameter: by - lowspeed
double iBy1;				//input parameter: by - highspeed
double iBz0;				//input parameter: bz - lowspeed
double iBz1;				//input parameter: bz - highspeed

double kx;					//Wavenumber of Boundary Layer disturbance
double dist;				//Amplitude of Boundary Layer disturbance
double plim;				//lower limit for thermal pressure
double pset;				//value for pressure to be set at when lower plim
*/

//Prototypes ########################################################################

void init_rho1 (double *[]);	void init_rho2 (double *[]);		//initialize starting variables rho, p, vx for 1 or 2 boundary layers in Parameter Study
void init_p1 (double *[]);	void init_p2 (double *[]);
void init_vx1 (double *[]);	void init_vx2 (double *[]);
void init_Bx1 (double *[]);	void init_Bx2 (double *[]);
void init_By1 (double *[]);	void init_By2 (double *[]);
void init_Bz1 (double *[]);	void init_Bz2 (double *[]);

void init_rho (double *[]); 		//initialize starting variables U,W
void init_rvx (double *[]);
void init_rvy (double *[]);
void init_rvz (double *[]); 
void init_e (double *[]);
void init_Bx (double *[]);
void init_By (double *[]);
void init_Bz (double *[]);
void init_p (double *[]);
void init_vx (double *[]);
void init_vy (double *[]);
void init_vz (double *[]);		

void cmaxX (double *[], double **[]);
void cmaxY (double *[], double **[]);


void Ls_rho (double *[]);		//calculate source terms
void Ls_rvx (double *[]);
void Ls_rvy (double *[]);
void Ls_rvz (double *[]); 
void Ls_e (double *[]);
void Ls_Bx (double *[]);
void Ls_By (double *[]);
void Ls_Bz (double *[]);

void Lx (double **[],double **[],double **[],double **[],double **[],double **[],double **[],double **[],double **[]);
void calcX_du_ww (double **[]);
void calcX_du_mm (double **[]);
void calcX_fL (double **[], double **[]);
void calcX_fR (double **[], double **[]);
void calcX_up (double **[]);
void calcX_cmax (double *[], double **[]);

void Ly (double **[],double **[],double **[],double **[],double **[],double **[],double **[],double **[],double **[]);
void calcY_du_ww (double **[]);
void calcY_du_mm (double **[]);
void calcY_fL (double **[], double **[]);
void calcY_fR (double **[], double **[]);
void calcY_up (double **[]);
void calcY_cmax (double *[], double **[]);

void clean_fieldCT (double *[], double *[]);
void clean_fluxCT (double *[], double *[]);
void clean_fluxCD (double *[], double *[]);
void clean_fieldCD (double *[], double *[]);

void upd_W (double **[],double);
void upd_old (double **[]);
void updBoundX_transmissive (double **[]);
void updBoundY_transmissive (double **[]);
void updBoundX_periodic (double **[]);
void updBoundY_periodic (double **[]);
void divergenceB (double *[]);

string intToString(int t);		//convert time to string for use in output files
string doubleToString(double t);		//convert time to string for use in output files




//main program ############################################################################


int main(void){

  cout<<endl<<"Solver for set of 2D MHD equations using TVD Lax-Friedrich scheme\n\n";
  
//---------- variables init ------
double vx;				//x Velocity of advection
double vy;				//y Velocity of advection
double x_0;				//for wave function: x -peak
double y_0;				//for wave function: y -peak
double cfl;				//courant number
double cfl_start;			//courant number for first time steps
double alphax;				//special factor for algorithm - to be calculated
double alphax2;				//special factor for algorithm - to be calculated
double alphay;				//special factor for algorithm - to be calculated
double alphay2;				//special factor for algorithm - to be calculated
double **divB;				//control for divB=0
    
//---------- file input ----------  

/*
  ifstream input ("input.txt");
  vector<double> variables;
  string var, is, remark;
  double value;
  while (input >> var >> is >> value >> remark) {
  cout << var << " = " << value << " " << remark << endl;
  variables.push_back(value);
  }
*/

//---------- write variables ----------  

/*
  tmax=variables[0];
  nx=variables[1];
  ny=variables[2];
  boxx=variables[3];
  boxy=variables[4];
  cfl_input=variables[5];
  x_boundaries=variables[6];
  y_boundaries=variables[7];
  cleaning=variables[8];
  cmethod=variables[9];
  kappa=variables[10];
  mu=variables[11];
  na=variables[12];
  thick=variables[13];
  iR0=variables[14];
  iR1=variables[15];
  iVx0=variables[16];
  iVx1=variables[17];
  iVy0=variables[18];
  iVy1=variables[19];
  iVz0=variables[20];
  iVz1=variables[21];
  iP0=variables[22];
  iP1=variables[23];
  iBx0=variables[24];
  iBx1=variables[25];
  iBy0=variables[26];
  iBy1=variables[27];
  iBz0=variables[28];
  iBz1=variables[29];
  kx=variables[30];
  dist=variables[31];
  nout=variables[32];
  limiter=variables[33];
// plim=variables[34];
// pset=variables[35];
*/

  dx=boxx/nx;
  dy=boxy/ny;
//  plim=0.;
//  pset=0.01;
  cfl_start=0.2;
  
  u=new double **[8];
  w=new double **[5];
  du=new double **[8];  
  fLn=new double **[8];  
  fRn=new double **[8];  
  up=new double **[8];  
  fLp=new double **[8];  
  fRp=new double **[8];  
  cmax=new double **[2];
  f=new double **[2];
  old=new double **[12];
  
  divB=new double *[nx+4];
  randnum=new double *[nx+4];
  for (int j=0;j<nx+4;++j){divB[j]= new double [ny+4]; randnum[j]= new double [ny+4];}

  for (int i=0;i<8;++i){ 
    u[i]= new double *[nx+4]; du[i]= new double *[nx+4]; fLn[i]= new double *[nx+4]; fRn[i]= new double *[nx+4]; up[i]= new double *[nx+4]; fLp[i]= new double *[nx+4]; fRp[i]= new double *[nx+4];  
  }	
  for (int i=0;i<8;++i) for (int j=0;j<nx+4;++j){ 
    u[i][j]= new double [ny+4]; du[i][j]= new double [ny+4]; fLn[i][j]= new double [ny+4]; fRn[i][j]= new double [ny+4]; up[i][j]= new double [ny+4]; fLp[i][j]= new double [ny+4]; fRp[i][j]= new double [ny+4];
  }
  for (int i=0;i<5;++i) w[i]= new double *[nx+4];
  for (int i=0;i<5;++i) for (int j=0;j<nx+4;++j) w[i][j]= new double [ny+4];
  for (int i=0;i<2;++i){ cmax[i]= new double *[nx+4]; f[i]= new double *[nx+4];}
  for (int i=0;i<2;++i) for (int j=0;j<nx+4;++j){ cmax[i][j]= new double [ny+4]; f[i][j]= new double [ny+4];}
  for (int i=0;i<12;++i) old[i]= new double *[nx+4];
  for (int i=0;i<12;++i) for (int j=0;j<nx+4;++j) old[i][j]= new double [ny+4];
  
 

//---------- initialization of starting parameters U[j][i], W[j][i] and output files ----------

  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) randnum[j][i] = (double)rand() / (double)RAND_MAX ;	//random number generator 

  (na==1)? init_rho1(u[0]):init_rho2(u[0]) ;
  (na==1)? init_p1(w[3]):init_p2(w[3]) ;
  (na==1)? init_vx1(w[0]):init_vx2(w[0]) ;
  (na==1)? init_Bx1(u[5]):init_Bx2(u[5]) ;
  (na==1)? init_By1(u[6]):init_By2(u[6]) ;
  (na==1)? init_Bz1(u[7]):init_Bz2(u[7]) ;

// for test cases:
//   init_rho(u[0]);   init_vx(w[0]);  init_p(w[3]);  init_Bx(u[5]);  init_By(u[6]);  init_Bz(u[7]); 
  init_vy(w[1]);  init_vz(w[2]);
  init_rvx(u[1]);  init_rvy(u[2]);  init_rvz(u[3]);  init_e(u[4]); 
  
  sprintf(dir,"data_rho%.2f_vx%.2f_phiB%.0f_kx%.2f_Lx%.0f_Ly%.0f", iR0, iVx0, atan (iBx0/iBz0) * 180/PI, kx, boxx, boxy);
  sprintf(buffer, "mkdir %s", dir);
  system(buffer);

  ofstream initial;

  for (int f=0; f<8; f++) {
      sprintf(buffer, "%s/initial_MHD_u%i.txt", dir, f);
      initial.open (buffer);
      for (int j=2;j<nx+2;++j){
          for (int i=2;i<ny+2;++i) initial<<u[f][j][i]<<" ";
          initial<<endl;
      }
      initial.close();
  }
  for (int f=0; f<5; f++) {
      sprintf(buffer, "%s/initial_MHD_w%i.txt", dir, f);
      initial.open (buffer);
      for (int j=2;j<nx+2;++j){
          for (int i=2;i<ny+2;++i) initial<<w[f][j][i]<<" ";
          initial<<endl;
      }
      initial.close();
  }

  

//---------- main loops: calling main functions and write files ----------

  ofstream data;
  int n=0;
  cfl=cfl_start;
  for (double t=0.;t<=tmax;t+=2.*dt){
    int tint=t;double cx=0.,cy=0.,divBmax=0.,pmin=(iP0<iP1) ? iP0:iP1,rmin=(iR0<iR1) ? iR0:iR1;
      cmaxX (cmax[0],u);
      cmaxY (cmax[1],u);
      divergenceB (divB);

  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i){ 	//calc in X
	cx=(cx<cmax[0][j][i]) ? cmax[0][j][i]:cx;
	cy=(cy<cmax[1][j][i]) ? cmax[1][j][i]:cy;
	divBmax=(divBmax<divB[j][i]) ? divB[j][i]:divBmax;
	pmin=(pmin<w[3][j][i]) ? pmin:w[3][j][i];
	rmin=(rmin<u[0][j][i]) ? rmin:u[0][j][i];
  }
    dt=cfl*min(dx/cx,dy/cy);
    cout << "n= "<< n << "	t= "<< t << "	dt= " << dt <<  "	divB= " << divBmax <<  "	p_min= " << pmin <<  "	rho_min= " << rmin << endl;
      
    if (pmin < plim){
        sprintf(buffer, "%s/pressure_problem_t%f_pmin_%f.txt", dir, t, pmin);
        data.open (buffer);
	for (int j=2;j<nx+2;++j){ 
	  for (int i=2;i<ny+2;++i) data<<w[3][j][i]<<" ";
	  data<<endl;
	}
	data.close();
    }
      
    upd_old(old);
      
// for test cases:
//     Ls_rho(u[0]);    Ls_rvx(u[1]);    Ls_rvy(u[2]);    Ls_rvz(u[3]);    Ls_e(u[4]);    Ls_Bx(u[5]);    Ls_By(u[6]);    Ls_Bz(u[7]);
    Lx (u,du,fLn,fRn,up,fLp,fRp,cmax,f); Ly (u,du,fLn,fRn,up,fLp,fRp,cmax,f); upd_W(w,t);

    if (cleaning==4) clean_fieldCT (u[5],u[6]);
    if (cleaning==2) clean_fluxCT (u[5],u[6]);
    if (cleaning==3) clean_fluxCD (u[5],u[6]);
    if (cleaning==1) clean_fieldCD (u[5],u[6]);
    if (cleaning>0){ 
      (x_boundaries==0)? updBoundX_transmissive (u):updBoundX_periodic (u) ;
      (y_boundaries==0)? updBoundY_transmissive (u):updBoundY_periodic (u) ;
      }
    upd_old(old);

// for test cases:
//     Ls_rho(u[0]);    Ls_rvx(u[1]);    Ls_rvy(u[2]);    Ls_rvz(u[3]);    Ls_e(u[4]);    Ls_Bx(u[5]);    Ls_By(u[6]);    Ls_Bz(u[7]);
    Ly (u,du,fLn,fRn,up,fLp,fRp,cmax,f); Lx (u,du,fLn,fRn,up,fLp,fRp,cmax,f); upd_W(w,t);


    if (cleaning==4) clean_fieldCT (u[5],u[6]);
    if (cleaning==2) clean_fluxCT (u[5],u[6]);
    if (cleaning==3) clean_fluxCD (u[5],u[6]);
    if (cleaning==1) clean_fieldCD (u[5],u[6]);
    if (cleaning>0){ 
      (x_boundaries==0)? updBoundX_transmissive (u):updBoundX_periodic (u) ;
      (y_boundaries==0)? updBoundY_transmissive (u):updBoundY_periodic (u) ;
      }
    Ls_rho(u[0]);    Ls_rvx(u[1]);    Ls_rvy(u[2]);    Ls_rvz(u[3]);    Ls_e(u[4]);    Ls_Bx(u[5]);    Ls_By(u[6]);    Ls_Bz(u[7]);

    
    n = n+2;
    cfl = (n>=101) ? cfl_input:cfl;
    
//         cout << "t= "<< t << "	dt= " << dt << " writing fluxes, etc..." << endl;    
   
    
    if ( fmod(t,tmax/nout) < 2*dt )
    {
        cout << "n= "<< n << "	t= "<< t << "	dt= " << dt << " writing file..." << endl;

        for (int f=0; f<8; f++) {
            sprintf(buffer, "%s/data_MHD_u%i_t%f_n%i.txt", dir, f, t, nx);
            data.open (buffer);
            for (int j=2;j<nx+2;++j){
                for (int i=2;i<ny+2;++i) data<<u[f][j][i]<<" ";
                data<<endl;
            }
            data.close();
        }
        for (int f=0; f<5; f++) {
            sprintf(buffer, "%s/data_MHD_w%i_t%f_n%i.txt", dir, f, t, nx);
            data.open (buffer);
            for (int j=2;j<nx+2;++j){
                for (int i=2;i<ny+2;++i) data<<w[f][j][i]<<" ";
                data<<endl;
            }
            data.close();
        }
    }

  }

  for (int f=0; f<8; f++) {
    sprintf(buffer, "%s/data_MHD_u%i_t%f_n%i.txt", dir, f, tmax, nx);
    data.open (buffer);
    for (int j=2;j<nx+2;++j){
        for (int i=2;i<ny+2;++i) data<<u[f][j][i]<<" ";
        data<<endl;
    }
    data.close();
  }
  for (int f=0; f<5; f++) {
    sprintf(buffer, "%s/data_MHD_w%i_t%f_n%i.txt", dir, f, tmax, nx);
    data.open (buffer);
    for (int j=2;j<nx+2;++j){
        for (int i=2;i<ny+2;++i) data<<w[f][j][i]<<" ";
        data<<endl;
    }
    data.close();
  }
    

	
for(int i = 0; i < 8; ++i){
  for(int j = 0; j < nx+4; ++j){ 
    delete [] u[i][j];
    delete [] du[i][j];
    delete [] fLn[i][j];
    delete [] fRn[i][j];
    delete [] up[i][j];
    delete [] fLp[i][j];
    delete [] fRp[i][j];
  }
  delete [] u[i];
  delete [] du[i];
  delete [] fLn[i];
  delete [] fRn[i];
  delete [] up[i];
  delete [] fLp[i];
  delete [] fRp[i];
} 
for(int i = 0; i < 5; ++i){
  for(int j = 0; j < nx+4; ++j) delete [] w[i][j];
  delete [] w[i];
}
for(int i = 0; i < 2; ++i){
  for(int j = 0; j < nx+4; ++j){ delete [] cmax[i][j];delete [] f[i][j];}
  delete [] cmax[i];delete [] f[i];
}
for(int i = 0; i < 12; ++i){ for(int j = 0; j < nx+4; ++j) delete [] old[i][j]; delete [] old[i];}

delete [] u;
delete [] w;
delete [] du;
delete [] fLn;
delete [] fRn;
delete [] up;
delete [] fLp;
delete [] fRp;
delete [] cmax;
delete [] f;
delete [] old;


u = NULL;
w = NULL;
du = NULL;
fLn = NULL;
fRn = NULL;
up = NULL;
fLp = NULL;
fRp = NULL;
cmax = NULL;
f = NULL;
old = NULL;

cout<<endl<<"##### MHD Simulation successful #####\n\n";
return 0;
}




//functions ######################################################################################


string intToString(int t){
std::string ch;
ostringstream out;
out << t; // Convert value into a string.
ch = out.str();
return ch;
} 

string doubleToString(double t){
std::string ch;
ostringstream out;
out << t; // Convert value into a string.
ch = out.str();
return ch;
} 


//---------- init functions: fill arrays of variables with initial conditions----------

void init_rho1 (double *u[]){

  for (int j=0;j<(nx+4);++j) for (int i=0;i<(ny+4);++i) u[j][i]=0.5*(iR0+iR1) + tanh((ny+4-2*i)/(2*thick))*0.5*(iR1-iR0) ;
  
}

void init_vx1 (double *u[]){   

  for (int j=0;j<(nx+4);++j) for (int i=0;i<(ny+4);++i) u[j][i]=(0.5*(iVx0+iVx1) + tanh((ny+4-2*i)/(2*thick))*0.5*(iVx1-iVx0));// + 2.*randnum[j][i]*dist - dist;

}

void init_p1 (double *u[]){  

  for (int j=0;j<(nx+4);++j) for (int i=0;i<(ny+4);++i) u[j][i]=0.5*(iP0+iP1) + tanh((ny+4-2*i)/(2*thick))*0.5*(iP1-iP0) ;
  
}

void init_Bx1 (double *u[]){  

  for (int j=0;j<(nx+4);++j) for (int i=0;i<(ny+4);++i) u[j][i]=0.5*(iBx0+iBx1) + tanh((ny+4-2*i)/(2*thick))*0.5*(iBx1-iBx0) ;
  
}

void init_By1 (double *u[]){  

  for (int j=0;j<(nx+4);++j) for (int i=0;i<(ny+4);++i) u[j][i]=0.5*(iBy0+iBy1) + tanh((ny+4-2*i)/(2*thick))*0.5*(iBy1-iBy0) ;
  
}

void init_Bz1 (double *u[]){  

  for (int j=0;j<(nx+4);++j) for (int i=0;i<(ny+4);++i) u[j][i]=0.5*(iBz0+iBz1) + tanh((ny+4-2*i)/(2*thick))*0.5*(iBz1-iBz0) ;
  
}


void init_rho2 (double *u[]){

  for (int j=0;j<(nx+4);++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.5*(iR0+iR1) + tanh((4*i-ny-4)/(4*thick))*0.5*(iR1-iR0) ;
  for (int j=0;j<(nx+4);++j) for (int i=(ny+4)/2;i<(ny+4);++i) u[j][i]=0.5*(iR0+iR1) + tanh((3*ny+12-4*i)/(4*thick))*0.5*(iR1-iR0) ;
  
}

void init_vx2 (double *u[]){   

  for (int j=0;j<(nx+4);++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=(0.5*(iVx0+iVx1) + tanh((4*i-ny-4)/(4*thick))*0.5*(iVx1-iVx0)) ;
  for (int j=0;j<(nx+4);++j) for (int i=(ny+4)/2;i<(ny+4);++i) u[j][i]=(0.5*(iVx0+iVx1) + tanh((3*ny+12-4*i)/(4*thick))*0.5*(iVx1-iVx0)) ;

}

void init_p2 (double *u[]){  

  for (int j=0;j<(nx+4);++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.5*(iP0+iP1) + tanh((4*i-ny-4)/(4*thick))*0.5*(iP1-iP0) ;
  for (int j=0;j<(nx+4);++j) for (int i=(ny+4)/2;i<(ny+4);++i) u[j][i]=0.5*(iP0+iP1) + tanh((3*ny+12-4*i)/(4*thick))*0.5*(iP1-iP0) ;
  
}

void init_Bx2 (double *u[]){  

  for (int j=0;j<(nx+4);++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.5*(iBx0+iBx1) + tanh((4*i-ny-4)/(4*thick))*0.5*(iBx1-iBx0) ;
  for (int j=0;j<(nx+4);++j) for (int i=(ny+4)/2;i<(ny+4);++i) u[j][i]=0.5*(iBx0+iBx1) + tanh((3*ny+12-4*i)/(4*thick))*0.5*(iBx1-iBx0) ;
  
}

void init_By2 (double *u[]){  

  for (int j=0;j<(nx+4);++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.5*(iBy0+iBy1) + tanh((4*i-ny-4)/(4*thick))*0.5*(iBy1-iBy0) ;
  for (int j=0;j<(nx+4);++j) for (int i=(ny+4)/2;i<(ny+4);++i) u[j][i]=0.5*(iBy0+iBy1) + tanh((3*ny+12-4*i)/(4*thick))*0.5*(iBy1-iBy0) ;
  
}

void init_Bz2 (double *u[]){  

  for (int j=0;j<(nx+4);++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.5*(iBz0+iBz1) + tanh((4*i-ny-4)/(4*thick))*0.5*(iBz1-iBz0) ;
  for (int j=0;j<(nx+4);++j) for (int i=(ny+4)/2;i<(ny+4);++i) u[j][i]=0.5*(iBz0+iBz1) + tanh((3*ny+12-4*i)/(4*thick))*0.5*(iBz1-iBz0) ;
  
}


void init_rho (double *u[]){

//------------------------- Tests -------------------------------------------  
    //Square function Test: Toth
//   for (int j=0;j<(nx+4);++j) for (int i=2;i<23;++i) u[j][i]=2.;
//   for (int j=0;j<(nx+4);++j) for (int i=23;i<(ny+4);++i) u[j][i]=0.5;

  //Semicircle function Test: Toth
/*  for (int j=0;j<(nx+4);++j) for (int i=0;i<5;++i) u[j][i]=1.;
  for (int j=0;j<(nx+4);++j) for (int i=5;i<35;++i) u[j][i]=1.+2.*sqrt(1.-(i-20.)*(i-20.)/(15.*15.));
  for (int j=0;j<(nx+4);++j) for (int i=35;i<(ny+4);++i) u[j][i]=1.;*/

  //Toro Test: 1
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.3*(ny+4);++i) u[j][i]=1.;
//   for (int j=0;j<(nx+4);++j) for (int i=0.3*(ny+4);i<(ny+4);++i) u[j][i]=0.125;

  //Toro Test: 2
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<(ny+4);++i) u[j][i]=1.0;

  //Toro Test: 3
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.5*(ny+4);++i) u[j][i]=1.;
//   for (int j=0;j<(nx+4);++j) for (int i=0.5*(ny+4);i<(ny+4);++i) u[j][i]=1.;

  //Toro Test: 4
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.4*(ny+4);++i) u[j][i]=5.99924;
//   for (int j=0;j<(nx+4);++j) for (int i=0.4*(ny+4);i<(ny+4);++i) u[j][i]=5.99242;

  //Toro Test: 5
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.8*(ny+4);++i) u[j][i]=1.;
//   for (int j=0;j<(nx+4);++j) for (int i=0.8*(ny+4);i<(ny+4);++i) u[j][i]=1.;


  //HD-2DTest: Liska 3
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.138;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.5323;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.5323;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=1.5;

  //HD-2DTest: Liska 4
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=1.1;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.5065;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.5065;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=1.1;

  //HD-2DTest: Liska 6
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=1.0;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=2.0;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=3.0;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=1.0;

  //HD-2DTest: Liska 12
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.8;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=1.0;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=1.0;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.5313;

  //HD-2DTest: Liska 15
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.8;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.5197;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.5313;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=1.0;

  //HD-2DTest: Liska 17
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=1.0625;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=2.0;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.5197;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=1.0;

  //HD-2DTest: Euler 4-contacts
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=1.0;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=2.0;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=3.0;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=1.0;

  //Toth MHD Shock Tube Test:
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.5*(ny+4);++i) u[j][i]=1.;
//   for (int j=0;j<(nx+4);++j) for (int i=0.5*(ny+4);i<(ny+4);++i) u[j][i]=0.125;

  //Toth MHD Shear Alf�n Waves Test:
//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=1.;

  //Toth MHD Vortex Test:
//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=25./9.;
  
//   //Ryu Jones 1D MHD Test: 1a
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.5*(ny+4);++i) u[j][i]=1.;
//   for (int j=0;j<(nx+4);++j) for (int i=0.5*(ny+4);i<(ny+4);++i) u[j][i]=1.;

  
}

void init_rvx (double *rvx[]){ 
  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) rvx[j][i]=u[0][j][i]*w[0][j][i];
}

void init_rvy (double *rvy[]){  
  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) rvy[j][i]=u[0][j][i]*w[1][j][i];
}

void init_rvz (double *rvz[]){  
  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) rvz[j][i]=u[0][j][i]*w[2][j][i];
}

void init_e (double *e[]){  
  //for given thermal pressure p:
  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) e[j][i]=w[3][j][i]/(kappa-1) + 0.5*u[0][j][i]*(w[0][j][i]*w[0][j][i]+w[1][j][i]*w[1][j][i]+w[2][j][i]*w[2][j][i]) + 0.5/mu*(u[5][j][i]*u[5][j][i]+u[6][j][i]*u[6][j][i]+u[7][j][i]*u[7][j][i]);

}

void init_Bx (double *u[]){  

  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=iBx0;

  
//------------------------- Tests -------------------------------------------  
  
  //   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=0.;

  //Toth MHD Shock Tube Test:
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.5*(ny+4);++i) u[j][i]=1.;
//   for (int j=0;j<(nx+4);++j) for (int i=0.5*(ny+4);i<(ny+4);++i) u[j][i]=-1.;

  //Toth MHD Shear Alf�n Waves Test:
//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=0.;

  //Toth MHD Vortex Test:
//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=sin(2.*(i-2)*dy);

  //Ryu Jones 1D MHD Test: 1a
//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=5./sqrt(4.*M_PI);


}

void init_By (double *u[]){  

  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=iBy0;

  
//------------------------- Tests -------------------------------------------  

//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=0.;

  //Toth MHD Shock Tube Test:
//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=0.75;

  //Toth MHD Shear Alf�n Waves Test:
//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=1.;

  //Toth MHD Vortex Test:
//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=-sin((j-2)*dx);

  //Ryu Jones 1D MHD Test: 1a
/*  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=5./sqrt(4.*M_PI);*/
  
}

void init_Bz (double *u[]){  

  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=iBz0;

  
//------------------------- Tests -------------------------------------------  

// for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=0.;
}

void init_vx (double *u[]){   

//------------------------- Tests -------------------------------------------  
  
  //   //Square function Test: Toth
//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=0.;

//HD-2DTest: Liska 3
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=1.206;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=1.206;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.;

//HD-2DTest: Liska 4
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.8939;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.8939;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.;

//HD-2DTest: Liska 6
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.5;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.5;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=-0.5;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=-0.5;

//HD-2DTest: Liska 12
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.7276;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.;

//HD-2DTest: Liska 15
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=-0.3;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=-0.3;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.4276;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=-0.3;

//HD-2DTest: Liska 17
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.2145;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=-0.3;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=-1.1259;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=-0.4;

  //HD-2DTest: Euler 4-contacts
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.5;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.5;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=-0.5;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=-0.5;

  //Toth MHD Shear Alf�n Waves Test:
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<(ny+4)/3.;++i) u[j][i]=0.;
//   for (int j=0;j<(nx+4);++j) for (int i=(ny+4)/3.;i<(ny+4)*2./3.;++i) u[j][i]=0.001;
//   for (int j=0;j<(nx+4);++j) for (int i=(ny+4)*2./3.;i<(ny+4);++i) u[j][i]=0.;

  //Toth MHD Vortex Test:
//  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=sin((i-2)*dy);
}

void init_vy (double *u[]){  

//------------------------- Random Disturbance ------------------------------  

//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=iVy0 + 2.*randnum[j][i]*dist - dist;   
  
  
//------------------------- Sinus Disturbance -------------------------------  

 for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=iVy0 + dist*sin((j-2.)/nx*(2.*M_PI))*exp(-0.25*(2.*i-ny-4.)*(2.*i-ny-4.)/thick/thick);
  
//------------------------- Tests -------------------------------------------  
  
//   //Square function Test: Toth
//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=1.;

  //Toro Test: 1
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.3*(ny+4);++i) u[j][i]=0.75;
//   for (int j=0;j<(nx+4);++j) for (int i=0.3*(ny+4);i<(ny+4);++i) u[j][i]=0.;

  //Toro Test: 2
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.5*(ny+4);++i) u[j][i]=-2.0;
//   for (int j=0;j<(nx+4);++j) for (int i=0.5*(ny+4);i<(ny+4);++i) u[j][i]=2.0;

  //Toro Test: 3
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.5*(ny+4);++i) u[j][i]=0.;
//   for (int j=0;j<(nx+4);++j) for (int i=0.5*(ny+4);i<(ny+4);++i) u[j][i]=0.;

  //Toro Test: 4
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.4*(ny+4);++i) u[j][i]=19.5975;
//   for (int j=0;j<(nx+4);++j) for (int i=0.4*(ny+4);i<(ny+4);++i) u[j][i]=-6.19633;

  //Toro Test: 5
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.8*(ny+4);++i) u[j][i]=-19.59745;
//   for (int j=0;j<(nx+4);++j) for (int i=0.8*(ny+4);i<(ny+4);++i) u[j][i]=-19.59745;

  //HD-2DTest: Liska 3
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=1.206;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=1.206;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.;

  //HD-2DTest: Liska 4
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.8939;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.8939;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.;

  //HD-2DTest: Liska 6
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=-0.75;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.75;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=-0.75;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.75;

  //HD-2DTest: Liska 12
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.0;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.7276;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.;

  //HD-2DTest: Liska 15
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.1;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=-0.6259;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.1;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.1;

  //HD-2DTest: Liska 17
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.0;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.0;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.0;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.0;

  //HD-2DTest: Euler 4-contacts
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=-0.75;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.75;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=-0.75;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.75;

  //Toth MHD Shock Tube Test:
//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=0.;

  //Toth MHD Shear Alf�n Waves Test:
//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=0.;

  //Toth MHD Vortex Test:
//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=-sin((j-2)*dx);

  //Ryu Jones 1D MHD Test: 1a
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.5*(ny+4);++i) u[j][i]=10.;
//   for (int j=0;j<(nx+4);++j) for (int i=0.5*(ny+4);i<(ny+4);++i) u[j][i]=-10.;

  
}

void init_vz (double *u[]){  

  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=iVz0;
  
//------------------------- Tests -------------------------------------------  
  
/*  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=0.;*/
}

void init_p (double *u[]){  

//------------------------- Tests -------------------------------------------  
  
  //Square function Test: Toth
//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=0.;

  //Toro Test: 1
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.3*(ny+4);++i) u[j][i]=1.;
//   for (int j=0;j<(nx+4);++j) for (int i=0.3*(ny+4);i<(ny+4);++i) u[j][i]=0.1;

  //Toro Test: 2
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<(ny+4);++i) u[j][i]=0.4;

  //Toro Test: 3
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.5*(ny+4);++i) u[j][i]=1000.;
//   for (int j=0;j<(nx+4);++j) for (int i=0.5*(ny+4);i<(ny+4);++i) u[j][i]=0.01;

  //Toro Test: 4
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.4*(ny+4);++i) u[j][i]=460.894;
//   for (int j=0;j<(nx+4);++j) for (int i=0.4*(ny+4);i<(ny+4);++i) u[j][i]=46.0950;

  //Toro Test: 5
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.8*(ny+4);++i) u[j][i]=1000.;
//   for (int j=0;j<(nx+4);++j) for (int i=0.8*(ny+4);i<(ny+4);++i) u[j][i]=0.01;

  //HD-2DTest: Liska 3
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.029;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.3;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.3;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=1.5;

  //HD-2DTest: Liska 4
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=1.1;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.35;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.35;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=1.1;

  //HD-2DTest: Liska 6
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=1.0;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=1.0;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=1.0;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=1.0;

  //HD-2DTest: Liska 12
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=1.0;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=1.0;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=1.0;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.4;

  //HD-2DTest: Liska 15
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.4;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.4;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.4;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=1.0;

  //HD-2DTest: Liska 17
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=0.4;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=1.0;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=0.4;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=1.0;

  //HD-2DTest: Euler 4-contacts
//   for (int j=0;j<(nx+4)/2;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=1.0;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=0;i<(ny+4)/2;++i) u[j][i]=1.0;
//   for (int j=0;j<(nx+4)/2;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=1.0;
//   for (int j=(nx+4)/2;j<nx+4;++j) for (int i=(ny+4)/2;i<ny+4;++i) u[j][i]=1.0;

  //Toth MHD Shock Tube Test:
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.5*(ny+4);++i) u[j][i]=1.;
//   for (int j=0;j<(nx+4);++j) for (int i=0.5*(ny+4);i<(ny+4);++i) u[j][i]=0.1;

  //Toth MHD Shear Alf�n Waves Test:
//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=10.0e-09;

  //Toth MHD Vortex Test:
//   for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) u[j][i]=1.666667;

  //Ryu Jones 1D MHD Test: 1a
//   for (int j=0;j<(nx+4);++j) for (int i=0;i<0.5*(ny+4);++i) u[j][i]=20.;
//   for (int j=0;j<(nx+4);++j) for (int i=0.5*(ny+4);i<(ny+4);++i) u[j][i]=1.;

}


//---------- update functions: calculation of boundaries and primitive variables----------

void updBoundX_transmissive (double **u[]){
//transmissive boundaries  
  for (int k=0;k<8;++k) for (int i=2;i<ny+2;++i){ u[k][0][i]=u[k][2][i]; u[k][1][i]=u[k][2][i]; u[k][nx+2][i]=u[k][nx+1][i]; u[k][nx+3][i]=u[k][nx+1][i]; }	//set boundaries X
/*  for (int k=0;k<8;++k){
  u[k][0][0]=u[k][2][2];u[k][0][1]=u[k][2][2];u[k][1][0]=u[k][2][2];u[k][1][1]=u[k][2][2];
  u[k][nx+2][1]=u[k][nx+1][2];u[k][nx+3][1]=u[k][nx+1][2];u[k][nx+2][0]=u[k][nx+1][2];u[k][nx+3][0]=u[k][nx+1][2];
  u[k][0][ny+2]=u[k][2][ny+1];u[k][1][ny+2]=u[k][2][ny+1];u[k][0][ny+3]=u[k][2][ny+1];u[k][1][ny+3]=u[k][2][ny+1];
  u[k][nx+2][ny+2]=u[k][nx+1][ny+1];u[k][nx+3][ny+2]=u[k][nx+1][ny+1];u[k][nx+2][ny+3]=u[k][nx+1][ny+1];u[k][nx+3][ny+3]=u[k][nx+1][ny+1];
	}*/
}

void updBoundY_transmissive (double **u[]){
//transmissive boundaries  
  for (int k=0;k<8;++k) for (int j=0;j<nx+4;++j){ u[k][j][0]=u[k][j][2]; u[k][j][1]=u[k][j][2]; u[k][j][ny+2]=u[k][j][ny+1]; u[k][j][ny+3]=u[k][j][ny+1]; }	//set boundaries Y
/*  for (int k=0;k<8;++k){
  u[k][0][0]=u[k][2][2];u[k][0][1]=u[k][2][2];u[k][1][0]=u[k][2][2];u[k][1][1]=u[k][2][2];
  u[k][nx+2][1]=u[k][nx+1][2];u[k][nx+3][1]=u[k][nx+1][2];u[k][nx+2][0]=u[k][nx+1][2];u[k][nx+3][0]=u[k][nx+1][2];
  u[k][0][ny+2]=u[k][2][ny+1];u[k][1][ny+2]=u[k][2][ny+1];u[k][0][ny+3]=u[k][2][ny+1];u[k][1][ny+3]=u[k][2][ny+1];
  u[k][nx+2][ny+2]=u[k][nx+1][ny+1];u[k][nx+3][ny+2]=u[k][nx+1][ny+1];u[k][nx+2][ny+3]=u[k][nx+1][ny+1];u[k][nx+3][ny+3]=u[k][nx+1][ny+1];
	}*/
}

void updBoundX_periodic (double **u[]){
//periodic boundaries
  for (int k=0;k<8;++k) for (int i=2;i<ny+2;++i){ u[k][0][i]=u[k][nx][i]; u[k][1][i]=u[k][nx+1][i]; u[k][nx+2][i]=u[k][2][i]; u[k][nx+3][i]=u[k][3][i]; }	//set boundaries X
}

void updBoundY_periodic (double **u[]){
//periodic boundaries
  for (int k=0;k<8;++k) for (int j=0;j<nx+4;++j){ u[k][j][0]=u[k][j][ny]; u[k][j][1]=u[k][j][ny+1]; u[k][j][ny+2]=u[k][j][2]; u[k][j][ny+3]=u[k][j][3]; }	//set boundaries Y  
}

void upd_W (double **w[],double t){  
// P update
  bool correcting=0;
  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) {
    w[3][j][i]= (kappa-1)* (u[4][j][i] - 0.5/u[0][j][i]*(u[3][j][i]*u[3][j][i]+u[1][j][i]*u[1][j][i]+u[2][j][i]*u[2][j][i]) - 0.5/mu*(u[5][j][i]*u[5][j][i]+u[6][j][i]*u[6][j][i]+u[7][j][i]*u[7][j][i]));

//------ NEGATIVE PRESSURE CORRECTION --------------------------------------------------   
    if (w[3][j][i]<plim){ 

      correcting=1;
      u[4][j][i]=u[4][j][i] -w[3][j][i]/(kappa-1) + pset/(kappa-1);
      w[3][j][i]= (kappa-1)* (u[4][j][i] - 0.5/u[0][j][i]*(u[3][j][i]*u[3][j][i]+u[1][j][i]*u[1][j][i]+u[2][j][i]*u[2][j][i]) - 0.5/mu*(u[5][j][i]*u[5][j][i]+u[6][j][i]*u[6][j][i]+u[7][j][i]*u[7][j][i]));
    }
  }
  if (correcting==1){
	ofstream data;
    sprintf(buffer, "%s/pressure_problem_t%f.txt", dir, t);
    data.open (buffer);
	for (int j=2;j<nx+2;++j){ 
	  for (int i=2;i<ny+2;++i) data<<w[3][j][i]<<" ";
	  data<<endl;
	}
	data.close(); 
	correcting=0;
  }



// vx update
  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) w[0][j][i]= u[1][j][i]/ u[0][j][i];
    // vy update
  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) w[1][j][i]= u[2][j][i]/ u[0][j][i];
    // vz update
  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) w[2][j][i]= u[3][j][i]/ u[0][j][i];
    // internal energy update
  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) w[4][j][i]= w[3][j][i]/((kappa-1.)*u[0][j][i]);
    
  
}

void upd_old (double **old[]){  
  for (int k=0;k<12;++k) for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i) old[k][j][i] = (k<8)? u[k][j][i]:w[k-8][j][i];  
}

//---------- Ls functions: calculation of source term updates----------

void Ls_rho (double *rho[]){
}

void Ls_rvx (double *rvx[]){
}

void Ls_rvy (double *rvy[]){
}

void Ls_rvz (double *rvz[]){
}

void Ls_e (double *e[]){
}

void Ls_Bx (double *Bx[]){
}

void Ls_By (double *By[]){
}

void Ls_Bz (double *Bz[]){
}

void Ls_vx (double *vx[]){  
}

void Ls_vy (double *vy[]){  
}

void Ls_vz (double *vz[]){  
}

void Ls_p (double *p[]){
}


//---------- Lx/Ly functions: calculation of U^(n+1) updates in X/Y----------

void Lx (double **u[],double **du[],double **fLn[],double **fRn[],double **up[],double **fLp[],double **fRp[],double **cmax[],double **f[]){
  ofstream data;
    
  //Limiter
    (limiter==0)? calcX_du_ww (du):calcX_du_mm (du);

  calcX_fL (fLn,u);

  calcX_fR (fRn,u);

  calcX_up (up);
    (x_boundaries==0)? updBoundX_transmissive (u):updBoundX_periodic (u);
    (y_boundaries==0)? updBoundY_transmissive (u):updBoundY_periodic (u);

  calcX_fL (fLp,up);

  calcX_fR (fRp,up);

  calcX_cmax (cmax[0],up);

  for (int j=1;j<nx+2;++j) for (int i=1;i<ny+2;++i) f[1][j][i] = 0.5* (fLp[6][j][i] + fRp[6][j][i+1]);
    
  for (int k=0;k<8;++k) for (int j=2;j<nx+2;++j) for (int i=2;i<ny+2;++i){
    u[k][j][i] = u[k][j][i] - 0.5*dt/dx*(fLp[k][j][i]+fRp[k][j+1][i]-fLp[k][j-1][i]-fRp[k][j][i] 
    -cmax[0][j][i]*(up[k][j+1][i]-0.5*du[k][j+1][i]-up[k][j][i]-0.5*du[k][j][i])
    +cmax[0][j-1][i]*(up[k][j][i]-0.5*du[k][j][i]-up[k][j-1][i]-0.5*du[k][j-1][i]));     
  }

  //Boundaries
    (x_boundaries==0)? updBoundX_transmissive (u):updBoundX_periodic (u);
    (y_boundaries==0)? updBoundY_transmissive (u):updBoundY_periodic (u);

}

void Ly (double **u[],double **du[],double **fLn[],double **fRn[],double **up[],double **fLp[],double **fRp[],double **cmax[],double **f[]){
  ofstream data;

  //Limiter
    (limiter==0)? calcY_du_ww (du):calcY_du_mm (du);

  calcY_fL (fLn,u);

  calcY_fR (fRn,u);

  calcY_up (up);
    (x_boundaries==0)? updBoundX_transmissive (u):updBoundX_periodic (u);
    (y_boundaries==0)? updBoundY_transmissive (u):updBoundY_periodic (u);

  calcY_fL (fLp,up);

  calcY_fR (fRp,up);

  calcY_cmax (cmax[1],up);

  for (int j=1;j<nx+2;++j) for (int i=1;i<ny+2;++i) f[0][j][i] = 0.5* (fLp[5][j][i] + fRp[5][j+1][i]);
  
  for (int k=0;k<8;++k) for (int j=2;j<nx+2;++j) for (int i=2;i<ny+2;++i){
    u[k][j][i] = u[k][j][i] - 0.5*dt/dy*(fLp[k][j][i]+fRp[k][j][i+1]-fLp[k][j][i-1]-fRp[k][j][i] 
    -cmax[1][j][i]*(up[k][j][i+1]-0.5*du[k][j][i+1]-up[k][j][i]-0.5*du[k][j][i])
    +cmax[1][j][i-1]*(up[k][j][i]-0.5*du[k][j][i]-up[k][j][i-1]-0.5*du[k][j][i-1]));     
  }

  //Boundaries
    (x_boundaries==0)? updBoundX_transmissive (u):updBoundX_periodic (u);
    (y_boundaries==0)? updBoundY_transmissive (u):updBoundY_periodic (u);

}


//---------- calc functions: calculation of necessary intermediate steps----------

void calcX_du_ww (double **du[]){
  double a;
  for (int k=0;k<8;++k) for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){ 	//calc in X
    a=(u[k][j-1][i]<=u[k][j][i]) ? 1.:-1.;
    du[k][j][i]=a*max(0.,min( abs(2.*(u[k][j][i]-u[k][j-1][i])) , min( a*2.*(u[k][j+1][i]-u[k][j][i]) , a*0.5*(u[k][j+1][i]-u[k][j-1][i]) ) ));
  }
}

void calcY_du_ww (double **du[]){   
  double a;
  for (int k=0;k<8;++k) for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){ 	//calc in X
    a=(u[k][j][i-1]<=u[k][j][i]) ? 1.:-1.;
    du[k][j][i]=a*max(0.,min( abs(2.*(u[k][j][i]-u[k][j][i-1])) , min( a*2.*(u[k][j][i+1]-u[k][j][i]) , a*0.5*(u[k][j][i+1]-u[k][j][i-1]) ) ));
  }
}

void calcX_du_mm (double **du[]){
  double a;
  for (int k=0;k<8;++k) for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){ 	//calc in X
    a=(u[k][j-1][i]<=u[k][j][i]) ? 1.:-1.;
    du[k][j][i]=a*max(0.,min( abs(u[k][j][i]-u[k][j-1][i]) , a*(u[k][j+1][i]-u[k][j][i]) ));
  }
}

void calcY_du_mm (double **du[]){   
  double a;
  for (int k=0;k<8;++k) for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){ 	//calc in X
    a=(u[k][j][i-1]<=u[k][j][i]) ? 1.:-1.;
    du[k][j][i]=a*max(0.,min( abs(u[k][j][i]-u[k][j][i-1]) , a*(u[k][j][i+1]-u[k][j][i]) ));
  }
}

void calcX_fL (double **f[],double **u[]){
// f[0] ... continuum eq. ------------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[0][j][i] = (u[1][j][i]+0.5*du[1][j][i]);
  }  
// f[1] ... momentum eq. in X --------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[1][j][i] = (u[1][j][i]+0.5*du[1][j][i])*(u[1][j][i]+0.5*du[1][j][i])/(u[0][j][i]+0.5*du[0][j][i])
   // +Pi (total pressure)
   +(kappa-1.)*((u[4][j][i]+0.5*du[4][j][i])
   -0.5*((u[1][j][i]+0.5*du[1][j][i])*(u[1][j][i]+0.5*du[1][j][i])+(u[2][j][i]+0.5*du[2][j][i])*(u[2][j][i]+0.5*du[2][j][i])+(u[3][j][i]+0.5*du[3][j][i])*(u[3][j][i]+0.5*du[3][j][i]))/(u[0][j][i]+0.5*du[0][j][i]))
   +0.5*(2.-kappa)*((u[5][j][i]+0.5*du[5][j][i])*(u[5][j][i]+0.5*du[5][j][i])+(u[6][j][i]+0.5*du[6][j][i])*(u[6][j][i]+0.5*du[6][j][i])+(u[7][j][i]+0.5*du[7][j][i])*(u[7][j][i]+0.5*du[7][j][i]))
   // end of Pi
   -(u[5][j][i]+0.5*du[5][j][i])*(u[5][j][i]+0.5*du[5][j][i]);    
  }
// f[2] ... momentum eq. in Y --------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[2][j][i] = (u[1][j][i]+0.5*du[1][j][i])*(u[2][j][i]+0.5*du[2][j][i])/(u[0][j][i]+0.5*du[0][j][i])
   - (u[5][j][i]+0.5*du[5][j][i])*(u[6][j][i]+0.5*du[6][j][i]);    
  }
// f[3] ... momentum eq. in Z --------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[3][j][i] = (u[1][j][i]+0.5*du[1][j][i])*(u[3][j][i]+0.5*du[3][j][i])/(u[0][j][i]+0.5*du[0][j][i])
   - (u[5][j][i]+0.5*du[5][j][i])*(u[7][j][i]+0.5*du[7][j][i]);    
  }
// f[4] ... energy eq. ---------------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[4][j][i] = ((u[1][j][i]+0.5*du[1][j][i])*
   // e+Pi
   ((u[4][j][i]+0.5*du[4][j][i])
   +(kappa-1.)*((u[4][j][i]+0.5*du[4][j][i])
   -0.5*((u[1][j][i]+0.5*du[1][j][i])*(u[1][j][i]+0.5*du[1][j][i])+(u[2][j][i]+0.5*du[2][j][i])*(u[2][j][i]+0.5*du[2][j][i])+(u[3][j][i]+0.5*du[3][j][i])*(u[3][j][i]+0.5*du[3][j][i]))/(u[0][j][i]+0.5*du[0][j][i]))
   +0.5*(2.-kappa)*((u[5][j][i]+0.5*du[5][j][i])*(u[5][j][i]+0.5*du[5][j][i])+(u[6][j][i]+0.5*du[6][j][i])*(u[6][j][i]+0.5*du[6][j][i])+(u[7][j][i]+0.5*du[7][j][i])*(u[7][j][i]+0.5*du[7][j][i])))
   // end of e+Pi
   -(u[5][j][i]+0.5*du[5][j][i])*((u[5][j][i]+0.5*du[5][j][i])*(u[1][j][i]+0.5*du[1][j][i])+(u[6][j][i]+0.5*du[6][j][i])*(u[2][j][i]+0.5*du[2][j][i])+(u[7][j][i]+0.5*du[7][j][i])*(u[3][j][i]+0.5*du[3][j][i])))
   /(u[0][j][i]+0.5*du[0][j][i]);    
  }
// f[5] ... B-field eq. in X ---------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[5][j][i] = 0.;    
  }
// f[6] ... B-field eq. in Y ---------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[6][j][i] = ((u[1][j][i]+0.5*du[1][j][i])*(u[6][j][i]+0.5*du[6][j][i]) 
   - (u[2][j][i]+0.5*du[2][j][i])*(u[5][j][i]+0.5*du[5][j][i]))
   /(u[0][j][i]+0.5*du[0][j][i]);    
  }
// f[7] ... B-field eq. in Z ---------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[7][j][i] = ((u[1][j][i]+0.5*du[1][j][i])*(u[7][j][i]+0.5*du[7][j][i]) 
   - (u[3][j][i]+0.5*du[3][j][i])*(u[5][j][i]+0.5*du[5][j][i]))
   /(u[0][j][i]+0.5*du[0][j][i]);     
  }

}

void calcY_fL (double **f[],double **u[]){ 
// f[0] ... continuum eq. ------------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[0][j][i] = (u[2][j][i]+0.5*du[2][j][i]);
  }  
// f[1] ... momentum eq. in X --------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[1][j][i] = (u[1][j][i]+0.5*du[1][j][i])*(u[2][j][i]+0.5*du[2][j][i])/(u[0][j][i]+0.5*du[0][j][i])
   -(u[5][j][i]+0.5*du[5][j][i])*(u[6][j][i]+0.5*du[6][j][i]);    
  }
// f[2] ... momentum eq. in Y --------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[2][j][i] = (u[2][j][i]+0.5*du[2][j][i])*(u[2][j][i]+0.5*du[2][j][i])/(u[0][j][i]+0.5*du[0][j][i])
   // +Pi (total pressure)
   +(kappa-1.)*((u[4][j][i]+0.5*du[4][j][i])
   -0.5*((u[1][j][i]+0.5*du[1][j][i])*(u[1][j][i]+0.5*du[1][j][i])+(u[2][j][i]+0.5*du[2][j][i])*(u[2][j][i]+0.5*du[2][j][i])+(u[3][j][i]+0.5*du[3][j][i])*(u[3][j][i]+0.5*du[3][j][i]))/(u[0][j][i]+0.5*du[0][j][i]))
   +0.5*(2.-kappa)*((u[5][j][i]+0.5*du[5][j][i])*(u[5][j][i]+0.5*du[5][j][i])+(u[6][j][i]+0.5*du[6][j][i])*(u[6][j][i]+0.5*du[6][j][i])+(u[7][j][i]+0.5*du[7][j][i])*(u[7][j][i]+0.5*du[7][j][i]))
   // end of Pi
   - (u[6][j][i]+0.5*du[6][j][i])*(u[6][j][i]+0.5*du[6][j][i]);    
  }
// f[3] ... momentum eq. in Z --------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[3][j][i] = (u[2][j][i]+0.5*du[2][j][i])*(u[3][j][i]+0.5*du[3][j][i])/(u[0][j][i]+0.5*du[0][j][i])
   - (u[6][j][i]+0.5*du[6][j][i])*(u[7][j][i]+0.5*du[7][j][i]);    
  }
// f[4] ... energy eq. ---------------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[4][j][i] = ((u[2][j][i]+0.5*du[2][j][i])*
   // e+Pi
   ((u[4][j][i]+0.5*du[4][j][i])
   +(kappa-1.)*((u[4][j][i]+0.5*du[4][j][i])
   -0.5*((u[1][j][i]+0.5*du[1][j][i])*(u[1][j][i]+0.5*du[1][j][i])+(u[2][j][i]+0.5*du[2][j][i])*(u[2][j][i]+0.5*du[2][j][i])+(u[3][j][i]+0.5*du[3][j][i])*(u[3][j][i]+0.5*du[3][j][i]))/(u[0][j][i]+0.5*du[0][j][i]))
   +0.5*(2.-kappa)*((u[5][j][i]+0.5*du[5][j][i])*(u[5][j][i]+0.5*du[5][j][i])+(u[6][j][i]+0.5*du[6][j][i])*(u[6][j][i]+0.5*du[6][j][i])+(u[7][j][i]+0.5*du[7][j][i])*(u[7][j][i]+0.5*du[7][j][i])))
   // end of e+Pi
   -(u[6][j][i]+0.5*du[6][j][i])*((u[5][j][i]+0.5*du[5][j][i])*(u[1][j][i]+0.5*du[1][j][i])+(u[6][j][i]+0.5*du[6][j][i])*(u[2][j][i]+0.5*du[2][j][i])+(u[7][j][i]+0.5*du[7][j][i])*(u[3][j][i]+0.5*du[3][j][i])))
   /(u[0][j][i]+0.5*du[0][j][i]);    
  }
// f[5] ... B-field eq. in X ---------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[5][j][i] = ((u[2][j][i]+0.5*du[2][j][i])*(u[5][j][i]+0.5*du[5][j][i]) 
   - (u[1][j][i]+0.5*du[1][j][i])*(u[6][j][i]+0.5*du[6][j][i]))
   /(u[0][j][i]+0.5*du[0][j][i]);    
  }
// f[6] ... B-field eq. in Y ---------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[6][j][i] = 0.;    
  }
// f[7] ... B-field eq. in Z ---------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[7][j][i] = ((u[2][j][i]+0.5*du[2][j][i])*(u[7][j][i]+0.5*du[7][j][i]) 
   - (u[3][j][i]+0.5*du[3][j][i])*(u[6][j][i]+0.5*du[6][j][i]))
   /(u[0][j][i]+0.5*du[0][j][i]);     
  }
}

void calcX_fR (double **f[],double **u[]){ 
// f[0] ... continuum eq. ------------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[0][j][i] = (u[1][j][i]-0.5*du[1][j][i]);
  }  
// f[1] ... momentum eq. in X --------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[1][j][i] = (u[1][j][i]-0.5*du[1][j][i])*(u[1][j][i]-0.5*du[1][j][i])/(u[0][j][i]-0.5*du[0][j][i])
   // +Pi (total pressure)
   +(kappa-1.)*((u[4][j][i]-0.5*du[4][j][i])
   -0.5*((u[1][j][i]-0.5*du[1][j][i])*(u[1][j][i]-0.5*du[1][j][i])+(u[2][j][i]-0.5*du[2][j][i])*(u[2][j][i]-0.5*du[2][j][i])+(u[3][j][i]-0.5*du[3][j][i])*(u[3][j][i]-0.5*du[3][j][i]))/(u[0][j][i]-0.5*du[0][j][i]))
   +0.5*(2.-kappa)*((u[5][j][i]-0.5*du[5][j][i])*(u[5][j][i]-0.5*du[5][j][i])+(u[6][j][i]-0.5*du[6][j][i])*(u[6][j][i]-0.5*du[6][j][i])+(u[7][j][i]-0.5*du[7][j][i])*(u[7][j][i]-0.5*du[7][j][i]))
   // end of Pi
   -(u[5][j][i]-0.5*du[5][j][i])*(u[5][j][i]-0.5*du[5][j][i]);    
  }
// f[2] ... momentum eq. in Y --------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[2][j][i] = (u[1][j][i]-0.5*du[1][j][i])*(u[2][j][i]-0.5*du[2][j][i])/(u[0][j][i]-0.5*du[0][j][i])
   - (u[5][j][i]-0.5*du[5][j][i])*(u[6][j][i]-0.5*du[6][j][i]);    
  }
// f[3] ... momentum eq. in Z --------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[3][j][i] = (u[1][j][i]-0.5*du[1][j][i])*(u[3][j][i]-0.5*du[3][j][i])/(u[0][j][i]-0.5*du[0][j][i])
   - (u[5][j][i]-0.5*du[5][j][i])*(u[7][j][i]-0.5*du[7][j][i]);    
  }
// f[4] ... energy eq. ---------------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[4][j][i] = ((u[1][j][i]-0.5*du[1][j][i])*
   // e+Pi
   ((u[4][j][i]-0.5*du[4][j][i])
   +(kappa-1.)*((u[4][j][i]-0.5*du[4][j][i])
   -0.5*((u[1][j][i]-0.5*du[1][j][i])*(u[1][j][i]-0.5*du[1][j][i])+(u[2][j][i]-0.5*du[2][j][i])*(u[2][j][i]-0.5*du[2][j][i])+(u[3][j][i]-0.5*du[3][j][i])*(u[3][j][i]-0.5*du[3][j][i]))/(u[0][j][i]-0.5*du[0][j][i]))
   +0.5*(2.-kappa)*((u[5][j][i]-0.5*du[5][j][i])*(u[5][j][i]-0.5*du[5][j][i])+(u[6][j][i]-0.5*du[6][j][i])*(u[6][j][i]-0.5*du[6][j][i])+(u[7][j][i]-0.5*du[7][j][i])*(u[7][j][i]-0.5*du[7][j][i])))
   // end of e+Pi
   -(u[5][j][i]-0.5*du[5][j][i])*((u[5][j][i]-0.5*du[5][j][i])*(u[1][j][i]-0.5*du[1][j][i])+(u[6][j][i]-0.5*du[6][j][i])*(u[2][j][i]-0.5*du[2][j][i])+(u[7][j][i]-0.5*du[7][j][i])*(u[3][j][i]-0.5*du[3][j][i])))
   /(u[0][j][i]-0.5*du[0][j][i]);    
  }
// f[5] ... B-field eq. in X ---------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[5][j][i] = 0.;    
  }
// f[6] ... B-field eq. in Y ---------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[6][j][i] = ((u[1][j][i]-0.5*du[1][j][i])*(u[6][j][i]-0.5*du[6][j][i]) 
   - (u[2][j][i]-0.5*du[2][j][i])*(u[5][j][i]-0.5*du[5][j][i]))
   /(u[0][j][i]-0.5*du[0][j][i]);    
  }
// f[7] ... B-field eq. in Z ---------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[7][j][i] = ((u[1][j][i]-0.5*du[1][j][i])*(u[7][j][i]-0.5*du[7][j][i]) 
   - (u[3][j][i]-0.5*du[3][j][i])*(u[5][j][i]-0.5*du[5][j][i]))
   /(u[0][j][i]-0.5*du[0][j][i]);     
  }

}

void calcY_fR (double **f[],double **u[]){ 
// f[0] ... continuum eq. ------------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[0][j][i] = (u[2][j][i]-0.5*du[2][j][i]);
  }  
// f[1] ... momentum eq. in X --------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[1][j][i] = (u[1][j][i]-0.5*du[1][j][i])*(u[2][j][i]-0.5*du[2][j][i])/(u[0][j][i]-0.5*du[0][j][i])
   -(u[5][j][i]-0.5*du[5][j][i])*(u[6][j][i]-0.5*du[6][j][i]);    
  }
// f[2] ... momentum eq. in Y --------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[2][j][i] = (u[2][j][i]-0.5*du[2][j][i])*(u[2][j][i]-0.5*du[2][j][i])/(u[0][j][i]-0.5*du[0][j][i])
   // +Pi (total pressure)
   +(kappa-1.)*((u[4][j][i]-0.5*du[4][j][i])
   -0.5*((u[1][j][i]-0.5*du[1][j][i])*(u[1][j][i]-0.5*du[1][j][i])+(u[2][j][i]-0.5*du[2][j][i])*(u[2][j][i]-0.5*du[2][j][i])+(u[3][j][i]-0.5*du[3][j][i])*(u[3][j][i]-0.5*du[3][j][i]))/(u[0][j][i]-0.5*du[0][j][i]))
   +0.5*(2.-kappa)*((u[5][j][i]-0.5*du[5][j][i])*(u[5][j][i]-0.5*du[5][j][i])+(u[6][j][i]-0.5*du[6][j][i])*(u[6][j][i]-0.5*du[6][j][i])+(u[7][j][i]-0.5*du[7][j][i])*(u[7][j][i]-0.5*du[7][j][i]))
   // end of Pi
   - (u[6][j][i]-0.5*du[6][j][i])*(u[6][j][i]-0.5*du[6][j][i]);    
  }
// f[3] ... momentum eq. in Z --------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[3][j][i] = (u[2][j][i]-0.5*du[2][j][i])*(u[3][j][i]-0.5*du[3][j][i])/(u[0][j][i]-0.5*du[0][j][i])
   - (u[6][j][i]-0.5*du[6][j][i])*(u[7][j][i]-0.5*du[7][j][i]);    
  }
// f[4] ... energy eq. ---------------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[4][j][i] = ((u[2][j][i]-0.5*du[2][j][i])*
   // e+Pi
   ((u[4][j][i]-0.5*du[4][j][i])
   +(kappa-1.)*((u[4][j][i]-0.5*du[4][j][i])
   -0.5*((u[1][j][i]-0.5*du[1][j][i])*(u[1][j][i]-0.5*du[1][j][i])+(u[2][j][i]-0.5*du[2][j][i])*(u[2][j][i]-0.5*du[2][j][i])+(u[3][j][i]-0.5*du[3][j][i])*(u[3][j][i]-0.5*du[3][j][i]))/(u[0][j][i]-0.5*du[0][j][i]))
   +0.5*(2.-kappa)*((u[5][j][i]-0.5*du[5][j][i])*(u[5][j][i]-0.5*du[5][j][i])+(u[6][j][i]-0.5*du[6][j][i])*(u[6][j][i]-0.5*du[6][j][i])+(u[7][j][i]-0.5*du[7][j][i])*(u[7][j][i]-0.5*du[7][j][i])))
   // end of e+Pi
   -(u[6][j][i]-0.5*du[6][j][i])*((u[5][j][i]-0.5*du[5][j][i])*(u[1][j][i]-0.5*du[1][j][i])+(u[6][j][i]-0.5*du[6][j][i])*(u[2][j][i]-0.5*du[2][j][i])+(u[7][j][i]-0.5*du[7][j][i])*(u[3][j][i]-0.5*du[3][j][i])))
   /(u[0][j][i]-0.5*du[0][j][i]);    
  }
// f[5] ... B-field eq. in X ---------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[5][j][i] = ((u[2][j][i]-0.5*du[2][j][i])*(u[5][j][i]-0.5*du[5][j][i]) 
   - (u[1][j][i]-0.5*du[1][j][i])*(u[6][j][i]-0.5*du[6][j][i]))
   /(u[0][j][i]-0.5*du[0][j][i]);    
  }
// f[6] ... B-field eq. in Y ---------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[6][j][i] = 0.;    
  }
// f[7] ... B-field eq. in Z ---------------------
  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i){
   f[7][j][i] = ((u[2][j][i]-0.5*du[2][j][i])*(u[7][j][i]-0.5*du[7][j][i]) 
   - (u[3][j][i]-0.5*du[3][j][i])*(u[6][j][i]-0.5*du[6][j][i]))
   /(u[0][j][i]-0.5*du[0][j][i]);     
  }
}

void calcX_up (double **up[]){
  for (int k=0;k<8;++k) for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i)up[k][j][i] = u[k][j][i];

  for (int k=0;k<8;++k) for (int j=2;j<nx+2;++j) for (int i=2;i<ny+2;++i){
    up[k][j][i] = u[k][j][i] - 0.5*dt/dx* (fLn[k][j][i] - fRn[k][j][i]);    
  }

}

void calcY_up (double **up[]){
  for (int k=0;k<8;++k) for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i)up[k][j][i] = u[k][j][i];
  
  for (int k=0;k<8;++k) for (int j=2;j<nx+2;++j) for (int i=2;i<ny+2;++i){
    up[k][j][i] = u[k][j][i] - 0.5*dt/dy* (fLn[k][j][i] - fRn[k][j][i]);
  }

}

void calcX_cmax (double *cmax[], double **u[]){
  double p[nx+4][ny+4],a[nx+4][ny+4],c=0.;
  for (int j=1;j<nx+2;++j) for (int i=1;i<ny+2;++i){
    p[j][i] = (kappa-1.)*(
    0.5*(u[4][j][i]+0.5*du[4][j][i]+u[4][j+1][i]-0.5*du[4][j+1][i])
    -0.5/(0.5*((u[0][j][i]+0.5*du[0][j][i])+(u[0][j+1][i]-0.5*du[0][j+1][i])))*
    (0.5*((u[1][j][i]+0.5*du[1][j][i])+(u[1][j+1][i]-0.5*du[1][j+1][i]))*0.5*((u[1][j][i]+0.5*du[1][j][i])+(u[1][j+1][i]-0.5*du[1][j+1][i])) 
    + 0.5*((u[2][j][i]+0.5*du[2][j][i])+(u[2][j+1][i]-0.5*du[2][j+1][i]))*0.5*((u[2][j][i]+0.5*du[2][j][i])+(u[2][j+1][i]-0.5*du[2][j+1][i]))
    + 0.5*((u[3][j][i]+0.5*du[3][j][i])+(u[3][j+1][i]-0.5*du[3][j+1][i]))*0.5*((u[3][j][i]+0.5*du[3][j][i])+(u[3][j+1][i]-0.5*du[3][j+1][i])))
    -0.5*
    ((0.5*(u[5][j][i]+0.5*du[5][j][i]+u[5][j+1][i]-0.5*du[5][j+1][i])) * (0.5*(u[5][j][i]+0.5*du[5][j][i]+u[5][j][i+1]-0.5*du[5][j][i+1]))
    +(0.5*(u[6][j][i]+0.5*du[6][j][i]+u[6][j+1][i]-0.5*du[6][j+1][i])) * (0.5*(u[6][j][i]+0.5*du[6][j][i]+u[6][j][i+1]-0.5*du[6][j][i+1]))
    +(0.5*(u[7][j][i]+0.5*du[7][j][i]+u[7][j+1][i]-0.5*du[7][j+1][i])) * (0.5*(u[7][j][i]+0.5*du[7][j][i]+u[7][j][i+1]-0.5*du[7][j][i+1])))
    );
//------ NEGATIVE PRESSURE CORRECTION --------------------------------------------------   
    if (p[j][i]<plim) p[j][i]=pset;
 
  }
  for (int j=1;j<nx+2;++j) for (int i=1;i<ny+2;++i){
    a[j][i] = (kappa*p[j][i] + 
    (0.5*(u[5][j][i]+0.5*du[5][j][i]+u[5][j+1][i]-0.5*du[5][j+1][i]))*(0.5*(u[5][j][i]+0.5*du[5][j][i]+u[5][j+1][i]-0.5*du[5][j+1][i])) + 
    (0.5*(u[6][j][i]+0.5*du[6][j][i]+u[6][j+1][i]-0.5*du[6][j+1][i]))*(0.5*(u[6][j][i]+0.5*du[6][j][i]+u[6][j+1][i]-0.5*du[6][j+1][i])) + 
    (0.5*(u[7][j][i]+0.5*du[7][j][i]+u[7][j+1][i]-0.5*du[7][j+1][i]))*(0.5*(u[7][j][i]+0.5*du[7][j][i]+u[7][j+1][i]-0.5*du[7][j+1][i])))
    /(0.5*(u[0][j][i]+0.5*du[0][j][i]+u[0][j+1][i]-0.5*du[0][j+1][i]));    
  }
  for (int j=1;j<nx+2;++j) for (int i=1;i<ny+2;++i){
    cmax[j][i] = abs(0.5*((u[1][j][i]+0.5*du[1][j][i])+(u[1][j+1][i]-0.5*du[1][j+1][i]))/(0.5*((u[0][j][i]+0.5*du[0][j][i])+(u[0][j+1][i]-0.5*du[0][j+1][i])))) + 
    1./sqrt(2.) * sqrt(a[j][i] + sqrt(a[j][i]*a[j][i] - 
    4.*kappa*p[j][i]*(0.5*(u[5][j][i]+0.5*du[5][j][i]+u[5][j+1][i]-0.5*du[5][j+1][i]))*(0.5*(u[5][j][i]+0.5*du[5][j][i]+u[5][j+1][i]-0.5*du[5][j+1][i]))
    /((0.5*(u[0][j][i]+0.5*du[0][j][i]+u[0][j+1][i]-0.5*du[0][j+1][i]))*(0.5*(u[0][j][i]+0.5*du[0][j][i]+u[0][j+1][i]-0.5*du[0][j+1][i])))));

    c=(c<cmax[j][i]) ? cmax[j][i]:c;
  }
  
if(cmethod==1)   for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i) cmax[j][i]=c;
}

void calcY_cmax (double *cmax[], double **u[]){
  double p[nx+4][ny+4],a[nx+4][ny+4],c=0.;
  for (int j=1;j<nx+2;++j) for (int i=1;i<ny+2;++i){

    p[j][i] = (kappa-1.)*(
    0.5*(u[4][j][i]+0.5*du[4][j][i]+u[4][j][i+1]-0.5*du[4][j][i+1])
    -0.5*(0.5*(u[0][j][i]+0.5*du[0][j][i]+u[0][j][i+1]-0.5*du[0][j][i+1]))*
    ((0.5*((u[1][j][i]+0.5*du[1][j][i])/(u[0][j][i]+0.5*du[0][j][i])+(u[1][j][i+1]-0.5*du[1][j][i+1])/(u[0][j][i+1]-0.5*du[0][j][i+1]))) * (0.5*((u[1][j][i]+0.5*du[1][j][i])/(u[0][j][i]+0.5*du[0][j][i])+(u[1][j][i+1]-0.5*du[1][j][i+1])/(u[0][j][i+1]-0.5*du[0][j][i+1])))
    + (0.5*((u[2][j][i]+0.5*du[2][j][i])/(u[0][j][i]+0.5*du[0][j][i])+(u[2][j][i+1]-0.5*du[2][j][i+1])/(u[0][j][i+1]-0.5*du[0][j][i+1]))) * (0.5*((u[2][j][i]+0.5*du[2][j][i])/(u[0][j][i]+0.5*du[0][j][i])+(u[2][j][i+1]-0.5*du[2][j][i+1])/(u[0][j][i+1]-0.5*du[0][j][i+1])))
    + (0.5*((u[3][j][i]+0.5*du[3][j][i])/(u[0][j][i]+0.5*du[0][j][i])+(u[3][j][i+1]-0.5*du[3][j][i+1])/(u[0][j][i+1]-0.5*du[0][j][i+1]))) * (0.5*((u[3][j][i]+0.5*du[3][j][i])/(u[0][j][i]+0.5*du[0][j][i])+(u[3][j][i+1]-0.5*du[3][j][i+1])/(u[0][j][i+1]-0.5*du[0][j][i+1])))
    )
    -0.5*
    ((0.5*(u[5][j][i]+0.5*du[5][j][i]+u[5][j][i+1]-0.5*du[5][j][i+1])) * (0.5*(u[5][j][i]+0.5*du[5][j][i]+u[5][j][i+1]-0.5*du[5][j][i+1])) 
    +(0.5*(u[6][j][i]+0.5*du[6][j][i]+u[6][j][i+1]-0.5*du[6][j][i+1])) * (0.5*(u[6][j][i]+0.5*du[6][j][i]+u[6][j][i+1]-0.5*du[6][j][i+1]))
    +(0.5*(u[7][j][i]+0.5*du[7][j][i]+u[7][j][i+1]-0.5*du[7][j][i+1])) * (0.5*(u[7][j][i]+0.5*du[7][j][i]+u[7][j][i+1]-0.5*du[7][j][i+1])))
    );  
//------ NEGATIVE PRESSURE CORRECTION --------------------------------------------------   
    if (p[j][i]<plim) p[j][i]=pset;

  }
  for (int j=1;j<nx+2;++j) for (int i=1;i<ny+2;++i){
    a[j][i] = (kappa*p[j][i] + 
    (0.5*(u[5][j][i]+0.5*du[5][j][i]+u[5][j][i+1]-0.5*du[5][j][i+1]))*(0.5*(u[5][j][i]+0.5*du[5][j][i]+u[5][j][i+1]-0.5*du[5][j][i+1])) + 
    (0.5*(u[6][j][i]+0.5*du[6][j][i]+u[6][j][i+1]-0.5*du[6][j][i+1]))*(0.5*(u[6][j][i]+0.5*du[6][j][i]+u[6][j][i+1]-0.5*du[6][j][i+1])) + 
    (0.5*(u[7][j][i]+0.5*du[7][j][i]+u[7][j][i+1]-0.5*du[7][j][i+1]))*(0.5*(u[7][j][i]+0.5*du[7][j][i]+u[7][j][i+1]-0.5*du[7][j][i+1])))
    /(0.5*(u[0][j][i]+0.5*du[0][j][i]+u[0][j][i+1]-0.5*du[0][j][i+1]));    
  }
  for (int j=1;j<nx+2;++j) for (int i=1;i<ny+2;++i){
    cmax[j][i] = abs(0.5*((u[2][j][i]+0.5*du[2][j][i])+(u[2][j][i+1]-0.5*du[2][j][i+1]))/(0.5*((u[0][j][i]+0.5*du[0][j][i])+(u[0][j][i+1]-0.5*du[0][j][i+1])))) + 
    1./sqrt(2.) * sqrt(a[j][i] + sqrt(a[j][i]*a[j][i] - 
    4.*kappa*p[j][i]*(0.5*(u[6][j][i]+0.5*du[6][j][i]+u[6][j][i+1]-0.5*du[6][j][i+1]))*(0.5*(u[6][j][i]+0.5*du[6][j][i]+u[6][j][i+1]-0.5*du[6][j][i+1]))
    /((0.5*(u[0][j][i]+0.5*du[0][j][i]+u[0][j][i+1]-0.5*du[0][j][i+1]))*(0.5*(u[0][j][i]+0.5*du[0][j][i]+u[0][j][i+1]-0.5*du[0][j][i+1])))));

    c=(c<cmax[j][i]) ? cmax[j][i]:c;
  }
  
if(cmethod==1)  for (int j=1;j<nx+3;++j) for (int i=1;i<ny+3;++i) cmax[j][i]=c;

}

void cmaxX (double *cmax[], double **u[]){
  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i)
  cmax[j][i]=abs(u[1][j][i]/u[0][j][i])+1./sqrt(2.) 
  * sqrt(
  (kappa*w[3][j][i] //(kappa-1.)*(u[4][j][i]-0.5*(u[1][j][i]*u[1][j][i]+u[2][j][i]*u[2][j][i]+u[3][j][i]*u[3][j][i])/u[0][j][i]-0.5*(u[5][j][i]*u[5][j][i]+u[6][j][i]*u[6][j][i]+u[7][j][i]*u[7][j][i]))
  +(u[5][j][i]*u[5][j][i]+u[6][j][i]*u[6][j][i]+u[7][j][i]*u[7][j][i]))/u[0][j][i]
  +sqrt( (kappa*w[3][j][i] //(kappa-1.)*(u[4][j][i]-0.5*(u[1][j][i]*u[1][j][i]+u[2][j][i]*u[2][j][i]+u[3][j][i]*u[3][j][i])/u[0][j][i]-0.5*(u[5][j][i]*u[5][j][i]+u[6][j][i]*u[6][j][i]+u[7][j][i]*u[7][j][i]))
  +(u[5][j][i]*u[5][j][i]+u[6][j][i]*u[6][j][i]+u[7][j][i]*u[7][j][i]))/u[0][j][i] * (kappa*w[3][j][i] //(kappa-1.)*(u[4][j][i]-0.5*(u[1][j][i]*u[1][j][i]+u[2][j][i]*u[2][j][i]+u[3][j][i]*u[3][j][i])/u[0][j][i]-0.5*(u[5][j][i]*u[5][j][i]+u[6][j][i]*u[6][j][i]+u[7][j][i]*u[7][j][i]))
  +(u[5][j][i]*u[5][j][i]+u[6][j][i]*u[6][j][i]+u[7][j][i]*u[7][j][i]))/u[0][j][i]
  -4*
  (kappa*w[3][j][i] //(kappa-1.)*(u[4][j][i]-0.5*(u[1][j][i]*u[1][j][i]+u[2][j][i]*u[2][j][i]+u[3][j][i]*u[3][j][i])/u[0][j][i]-0.5*(u[5][j][i]*u[5][j][i]+u[6][j][i]*u[6][j][i]+u[7][j][i]*u[7][j][i]))
  *(u[5][j][i]*u[5][j][i]))/(u[0][j][i]*u[0][j][i])
  ));
  
}

void cmaxY (double *cmax[], double **u[]){
  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i)
  cmax[j][i]=abs(u[2][j][i]/u[0][j][i])+1./sqrt(2.) 
  * sqrt(
  (kappa*w[3][j][i] //(kappa-1.)*(u[4][j][i]-0.5*(u[1][j][i]*u[1][j][i]+u[2][j][i]*u[2][j][i]+u[3][j][i]*u[3][j][i])/u[0][j][i]-0.5*(u[5][j][i]*u[5][j][i]+u[6][j][i]*u[6][j][i]+u[7][j][i]*u[7][j][i]))
  +(u[5][j][i]*u[5][j][i]+u[6][j][i]*u[6][j][i]+u[7][j][i]*u[7][j][i]))/u[0][j][i]
  +sqrt( (kappa*w[3][j][i] //(kappa-1.)*(u[4][j][i]-0.5*(u[1][j][i]*u[1][j][i]+u[2][j][i]*u[2][j][i]+u[3][j][i]*u[3][j][i])/u[0][j][i]-0.5*(u[5][j][i]*u[5][j][i]+u[6][j][i]*u[6][j][i]+u[7][j][i]*u[7][j][i]))
  +(u[5][j][i]*u[5][j][i]+u[6][j][i]*u[6][j][i]+u[7][j][i]*u[7][j][i]))/u[0][j][i] * (kappa*w[3][j][i] //(kappa-1.)*(u[4][j][i]-0.5*(u[1][j][i]*u[1][j][i]+u[2][j][i]*u[2][j][i]+u[3][j][i]*u[3][j][i])/u[0][j][i]-0.5*(u[5][j][i]*u[5][j][i]+u[6][j][i]*u[6][j][i]+u[7][j][i]*u[7][j][i]))
  +(u[5][j][i]*u[5][j][i]+u[6][j][i]*u[6][j][i]+u[7][j][i]*u[7][j][i]))/u[0][j][i]
  -4*
  (kappa*w[3][j][i] //(kappa-1.)*(u[4][j][i]-0.5*(u[1][j][i]*u[1][j][i]+u[2][j][i]*u[2][j][i]+u[3][j][i]*u[3][j][i])/u[0][j][i]-0.5*(u[5][j][i]*u[5][j][i]+u[6][j][i]*u[6][j][i]+u[7][j][i]*u[7][j][i]))
  *(u[6][j][i]*u[6][j][i]))/(u[0][j][i]*u[0][j][i])
  ));
  
}

void clean_fluxCT (double *bx[], double *by[]){
  for (int j=2;j<nx+2;++j) for (int i=2;i<ny+2;++i){
  bx[j][i]=bx[j][i]
//   bold[0][j][i] - dt/dy/16.*
  
  ;
  by[j][i]=by[j][i]
//   bold[1][j][i] - dt/dx/16.
  ;
  }
}

void clean_fluxCD (double *bx[], double *by[]){
  for (int j=2;j<nx+2;++j) for (int i=2;i<ny+2;++i){
  bx[j][i]=bx[j][i]
//   bold[0][j][i] - dt/dy/8.*
//   
  ;
  by[j][i]=by[j][i]
//   bold[1][j][i] + dt/dx/8.*
//   
  ;
  }  
}

void clean_fieldCD (double *bx[], double *by[]){
  double bx_neu[nx+4][ny+4],by_neu[nx+4][ny+4];
  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i){bx_neu[j][i]=bx[j][i];by_neu[j][i]=by[j][i];}
  for (int j=2;j<nx+2;++j) for (int i=2;i<ny+2;++i){
  bx_neu[j][i]=
  old[5][j][i] + 0.25*dt/dy* (
  (old[8][j][i+1]*old[6][j][i+1]-old[9][j][i+1]*old[5][j][i+1]+w[0][j][i+1]*u[6][j][i+1]-w[1][j][i+1]*u[5][j][i+1]) - 
  (old[8][j][i-1]*old[6][j][i-1]-old[9][j][i-1]*old[5][j][i-1]+w[0][j][i-1]*u[6][j][i-1]-w[1][j][i-1]*u[5][j][i-1]))
  ;
  by_neu[j][i]=
  old[6][j][i] - 0.25*dt/dx* (
  (old[8][j+1][i]*old[6][j+1][i]-old[9][j+1][i]*old[5][j+1][i]+w[0][j+1][i]*u[6][j+1][i]-w[1][j+1][i]*u[5][j+1][i]) -
  (old[8][j-1][i]*old[6][j-1][i]-old[9][j-1][i]*old[5][j-1][i]+w[0][j-1][i]*u[6][j-1][i]-w[1][j-1][i]*u[5][j-1][i])) 
  ;
  }
  for (int j=2;j<nx+2;++j) for (int i=2;i<ny+2;++i){bx[j][i]=bx_neu[j][i];by[j][i]=by_neu[j][i];}

}

void clean_fieldCT (double *bx[], double *by[]){
  double pp[4][nx+4][ny+4],pm[4][nx+4][ny+4],mp[4][nx+4][ny+4],mm[4][nx+4][ny+4],bx_neu[nx+4][ny+4],by_neu[nx+4][ny+4];
  for (int j=0;j<nx+4;++j) for (int i=0;i<ny+4;++i){bx_neu[j][i]=bx[j][i];by_neu[j][i]=by[j][i];}
  for(int j=2;j<nx+2;++j) for (int i=2;i<ny+2;++i){
   pp[0][j][i]= 0.125*(old[8][j][i]+old[8][j+1][i]+old[8][j][i+1]+old[8][j+1][i+1]+w[0][j][i]+w[0][j+1][i]+w[0][j][i+1]+w[0][j+1][i+1]);
   pp[1][j][i]= 0.125*(old[9][j][i]+old[9][j+1][i]+old[9][j][i+1]+old[9][j+1][i+1]+w[1][j][i]+w[1][j+1][i]+w[1][j][i+1]+w[1][j+1][i+1]);
   pp[2][j][i]= 0.125*(old[5][j][i]+old[5][j+1][i]+old[5][j][i+1]+old[5][j+1][i+1]+u[5][j][i]+u[5][j+1][i]+u[5][j][i+1]+u[5][j+1][i+1]);
   pp[3][j][i]= 0.125*(old[6][j][i]+old[6][j+1][i]+old[6][j][i+1]+old[6][j+1][i+1]+u[6][j][i]+u[6][j+1][i]+u[6][j][i+1]+u[6][j+1][i+1]);

   pm[0][j][i]= 0.125*(old[8][j][i]+old[8][j+1][i]+old[8][j][i-1]+old[8][j+1][i-1]+w[0][j][i]+w[0][j+1][i]+w[0][j][i-1]+w[0][j+1][i-1]);
   pm[1][j][i]= 0.125*(old[9][j][i]+old[9][j+1][i]+old[9][j][i-1]+old[9][j+1][i-1]+w[1][j][i]+w[1][j+1][i]+w[1][j][i-1]+w[1][j+1][i-1]);
   pm[2][j][i]= 0.125*(old[5][j][i]+old[5][j+1][i]+old[5][j][i-1]+old[5][j+1][i-1]+u[5][j][i]+u[5][j+1][i]+u[5][j][i-1]+u[5][j+1][i-1]);
   pm[3][j][i]= 0.125*(old[6][j][i]+old[6][j+1][i]+old[6][j][i-1]+old[6][j+1][i-1]+u[6][j][i]+u[6][j+1][i]+u[6][j][i-1]+u[6][j+1][i-1]);

   mp[0][j][i]= 0.125*(old[8][j][i]+old[8][j-1][i]+old[8][j][i+1]+old[8][j-1][i+1]+w[0][j][i]+w[0][j-1][i]+w[0][j][i+1]+w[0][j-1][i+1]);
   mp[1][j][i]= 0.125*(old[9][j][i]+old[9][j-1][i]+old[9][j][i+1]+old[9][j-1][i+1]+w[1][j][i]+w[1][j-1][i]+w[1][j][i+1]+w[1][j-1][i+1]);
   mp[2][j][i]= 0.125*(old[5][j][i]+old[5][j-1][i]+old[5][j][i+1]+old[5][j-1][i+1]+u[5][j][i]+u[5][j-1][i]+u[5][j][i+1]+u[5][j-1][i+1]);
   mp[3][j][i]= 0.125*(old[6][j][i]+old[6][j-1][i]+old[6][j][i+1]+old[6][j-1][i+1]+u[6][j][i]+u[6][j-1][i]+u[6][j][i+1]+u[6][j-1][i+1]);

   mm[0][j][i]= 0.125*(old[8][j][i]+old[8][j-1][i]+old[8][j][i-1]+old[8][j-1][i-1]+w[0][j][i]+w[0][j-1][i]+w[0][j][i-1]+w[0][j-1][i-1]);
   mm[1][j][i]= 0.125*(old[9][j][i]+old[9][j-1][i]+old[9][j][i-1]+old[9][j-1][i-1]+w[1][j][i]+w[1][j-1][i]+w[1][j][i-1]+w[1][j-1][i-1]);
   mm[2][j][i]= 0.125*(old[5][j][i]+old[5][j-1][i]+old[5][j][i-1]+old[5][j-1][i-1]+u[5][j][i]+u[5][j-1][i]+u[5][j][i-1]+u[5][j-1][i-1]);
   mm[3][j][i]= 0.125*(old[6][j][i]+old[6][j-1][i]+old[6][j][i-1]+old[6][j-1][i-1]+u[6][j][i]+u[6][j-1][i]+u[6][j][i-1]+u[6][j-1][i-1]);
  }
    
  for (int j=2;j<nx+2;++j) for (int i=2;i<ny+2;++i){
  bx_neu[j][i]=
  old[5][j][i] + 0.5*dt/dy* (pp[0][j][i]*pp[3][j][i]-pp[1][j][i]*pp[2][j][i] + mp[0][j][i]*mp[3][j][i]-mp[1][j][i]*mp[2][j][i] - pm[0][j][i]*pm[3][j][i]-pm[1][j][i]*pm[2][j][i] - mm[0][j][i]*mm[3][j][i]-mm[1][j][i]*mm[2][j][i]);

  by_neu[j][i]=
  old[6][j][i] - 0.5*dt/dx* (pp[0][j][i]*pp[3][j][i]-pp[1][j][i]*pp[2][j][i] - mp[0][j][i]*mp[3][j][i]-mp[1][j][i]*mp[2][j][i] + pm[0][j][i]*pm[3][j][i]-pm[1][j][i]*pm[2][j][i] - mm[0][j][i]*mm[3][j][i]-mm[1][j][i]*mm[2][j][i]);
  }
  for (int j=2;j<nx+2;++j) for (int i=2;i<ny+2;++i){bx[j][i]=bx_neu[j][i];by[j][i]=by_neu[j][i];}


}

void divergenceB (double *divB[]){  
  for (int j=2;j<nx+2;++j) for (int i=2;i<ny+2;++i)
    divB[j][i]=0.5*((u[5][j+1][i]-u[5][j-1][i])/dx+(u[6][j][i+1]-u[6][j][i-1])/dx);

}

