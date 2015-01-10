/*
 *  LLGlangevin3D.cpp
 *  
 *
 *  Ondrej Hovorka (o.hovorka@soton.ac.uk).
 *  (3. October 2013)
 * 
 *  This is a basic code which integrates the set of coupled stochastic 
 *  landau-lifshitz-gilbert (LLG) equation with noise, describing a system 
 *  of N coupled spins i = 1,...,N. At the moment the spin system sits on the 3D
 *  square lattice and relates the classical Heisenberg model with nearest neighbor
 *  in teraction. The normalization of the LLG euqation for i-th spin in the system 
 *  is as follows:
 *
 *  ds_i/dt = -gamma/(1+alpha^2)/mu_s * [ s_i x (H+h) + alpha * s_i x (s_i x (H+h)) ]
 *
 *  s_i, H, h are vectors. 
 * 
 *  s_i - dimensionless magnetization vector of i-th spin
 *  H - effective field in the units of energy [J]
 *  h - stochastic field in the units of energy [J]
 *  gama - gyromagnetic factor in [1/(Ts)]
 *  mu_s - magnetic moment in [J/T], e.g. Bohr magneton
 *  alpha - dimensionless damping factor
 *
 *  The effective field H is the sum of external fiels, exchange field, and anisotropy
 *  field, following the energy function:
 *
 *  E = -\sum_<ij> J_ij*s_i*s_j - \sum_i K_i.s_i - \sum_i H.s_i
 *  
 *  where K is the uniaxial anisotropy vector and J_ij is the interaction coupling
 *  between the nearest neighbor spins.
 *  
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <time.h>
#include <stdlib.h>

#include "randomc.h"
#include "stocc.h"

using namespace std;

void normalise_vector(vector<double> &v)
{
    /* 
       Normalise a vector. 
    */
    double norm;
    norm = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    v[0] = v[0]/norm;
    v[1] = v[1]/norm;
    v[2] = v[2]/norm;
}

vector<double> heff_eval_cub(vector<double> &sx, vector<double> &sy, vector<double> &sz,
                             vector<double> &s, vector<double> & h,
                             vector<double> &param, long &lj, long & num_neigh)
{
  /*
   * Calculate the effective field on a particle with cubic anisotropy K
   * in a local field B and acted on by an exhange field of strength J.
   */
  double J = param[4];
  double K1 = param[5];
  double K2 = param[6];
  double B = param[7];

  vector<double> h_int(3, 0.0);
  for(int l = 0; l < num_neigh; l++) //the total exchange field
  {
    h_int[0] += J*sx[lj+l];
    h_int[1] += J*sy[lj+l];
    h_int[2] += J*sz[lj+l];
  }

  vector<double> h_eff(3,0.0); //cubic anistoropy contribution defined here
  h_eff[0] = -2.0*s[0]*(K1*s[1]*s[1] + K1*s[2]*s[2] + K2*s[1]*s[1]*s[2]*s[2]);
  h_eff[1] = -2.0*s[1]*(K1*s[0]*s[0] + K1*s[2]*s[2] + K2*s[0]*s[0]*s[2]*s[2]);
  h_eff[2] = -2.0*s[2]*(K1*s[0]*s[0] + K1*s[1]*s[1] + K2*s[0]*s[0]*s[1]*s[1]);
  for(long l=0; l<2; l++) //overall sum for the effective field
  {
    h_eff[l] += h_int[l] + B*h[l];
  }

  return h_eff;
}

vector<double> heff_eval_uni(vector<double> &sx, vector<double> &sy, vector<double> &sz,
                         vector<double> &s, vector<double> &n, vector<double> &h, 
                         vector<double> &param, long &lj, long &num_neigh)
{
    /*
       Calculate the effective field on a particle with uniaxial anistoropy K 
       in a local field B and acted on by an exchange field of strength J. 
    */
    double J = param[4];
    double K1 = param[5];
    // ......param[6] is not used here
    double B = param[7];

    vector<double> h_int(3, 0.0);
    for(int l=0; l<num_neigh; l++) // the total exchange field 
    {
        h_int[0] += J*sx[lj+l];
        h_int[1] += J*sy[lj+l];
        h_int[2] += J*sz[lj+l];
    }

    double dotp=0.0; //dot product
    for(long l=0; l<3; l++) // the anisotropy field contribution
    {
        dotp += s[l]*n[l];
    }

    vector<double> h_eff(3, 0.0);
    for(long l=0; l<3; l++) // overall sum for the effective field 
    {
        h_eff[l] = h_int[l] + 2.0*K1*dotp*n[l] + B*h[l];
    }

    return h_eff;
}

vector<double> acoef(vector<double> &s, vector<double> &h_eff, 
                     vector<double> &param)
{
    /*
        'a' coefficients required by the Heun Scheme. 
    */
    double c1 = param[0];
    double c2 = param[1];
    double alpha = param[2];
    double dt = param[3];
  
    vector<double> a(3);

    a[0] = -c1*(s[1]*h_eff[2]-s[2]*h_eff[1])*dt
           -alpha*c1*(s[0]*(s[1]*h_eff[1]+s[2]*h_eff[2]))*dt
           +alpha*c1*(h_eff[0]*(s[1]*s[1]+s[2]*s[2]))*dt;
             
    a[1] = -c1*(s[2]*h_eff[0]-s[0]*h_eff[2])*dt
           -alpha*c1*(s[1]*(s[0]*h_eff[0]+s[2]*h_eff[2]))*dt
           +alpha*c1*(h_eff[1]*(s[0]*s[0]+s[2]*s[2]))*dt;

    a[2] = -c1*(s[0]*h_eff[1]-s[1]*h_eff[0])*dt
           -alpha*c1*(s[2]*(s[0]*h_eff[0]+s[1]*h_eff[1]))*dt
           +alpha*c1*(h_eff[2]*(s[0]*s[0]+s[1]*s[1]))*dt;

    return a;
}

vector<double> bcoef(vector<double> &s, vector<double> &h_eff, 
                     vector <double> &wi, vector<double> &param)
{
    /*
        'b' coefficients required by the Heun Scheme. 
    */
    double c1 = param[0];
    double c2 = param[1];
    double alpha = param[2];
    double dt = param[3];

    vector<double> b(3);
  
    dt = sqrt(dt);

    b[0] = alpha*c2*(s[1]*s[1]+s[2]*s[2])*dt*wi[0]
           +c2*(s[2]-alpha*s[0]*s[1])*dt*wi[1]
           -c2*(s[1]+alpha*s[0]*s[2])*dt*wi[2];
             
    b[1] = -c2*(s[2]+alpha*s[0]*s[1])*dt*wi[0]
           +alpha*c2*(s[0]*s[0]+s[2]*s[2])*dt*wi[1]
           +c2*(s[0]-alpha*s[1]*s[2])*dt*wi[2];

    b[2] = c2*(s[1]-alpha*s[1]*s[2])*dt*wi[0]
           -c2*(s[0]+alpha*s[1]*s[2])*dt*wi[1]
           +alpha*c2*(s[0]*s[0]+s[1]*s[1])*dt*wi[2];

    return b;
}

vector<long> lattice3D_periodic(long &num_part_side, long &num_neigh)
{
    /* 
       Generate 3D lattice with periodic boundary conditions. Return the 
       neighbor list.  
    */

    long t;
    long num_part = num_part_side*num_part_side*num_part_side;

    vector<long> neigh_list(num_part*num_neigh, 0);

    for(long i=0; i<num_part_side; i++)
    {
        for(long j=0; j<num_part_side; j++)
        {
            for(long k=0; k<num_part_side; k++)
            { 
                /*
                t = k*num_part_side*num_part_side+
                    j*num_part_side+i;

                is equivalent to building t from s[i][j][k], i.e. i is the fastest
                varying dimension in the t-vector and k is the most slowly varying 
                dimension. 

                - specifying the order of "for" loops changes the order in which the
                neighbor list is being built but does not change the spin ordering 
                along the lattice. To change the spin ordering we could specify, e.g.:
            
                t = i*num_part_side*num_part_side+
                    j*num_part_side+k;

                when the dimension i would be most slowly varyin within the array 
                indexed by t.
                */ 


                t = k*num_part_side*num_part_side+
                    j*num_part_side+i;

                neigh_list[num_neigh*t+0] = k*num_part_side*num_part_side+
                                   j*num_part_side+(i-1);

                neigh_list[num_neigh*t+1] = k*num_part_side*num_part_side+
                                   j*num_part_side+(i+1);

                neigh_list[num_neigh*t+2] = k*num_part_side*num_part_side+
                                   (j-1)*num_part_side+i;

                neigh_list[num_neigh*t+3] = k*num_part_side*num_part_side+
                                   (j+1)*num_part_side+i;

                neigh_list[num_neigh*t+4] = (k-1)*num_part_side*num_part_side+
                                   j*num_part_side+i;

                neigh_list[num_neigh*t+5] = (k+1)*num_part_side*num_part_side+
                                   j*num_part_side+i;
            }
        }
    }

    return neigh_list;
}


/*****************************************************
***                 Main Program                   ***
******************************************************/

int main(int argc, char* argv[])
{
    //write out a value however many of these timesteps
    long n_steps = 1e9;
    int data_frequency = (int)(n_steps/1e5);
    double dt = 1.0e-14;//atof(argv[4]);//1.0e-12; //[s] integration timestep

    // Basic model parameters and their units
    double mu0 = 4.0*3.14e-7;   // T/(A/m) [magnetic permeability]
    double kB = 1.38e-23;       // J/K [Boltzmann constant]
    double mu_s;                // J/T
    double mu_b = 9.274e-24;    // J/T [magnetic moment, Bohr magneton]
    mu_s = 66483.0*mu_b;
    double gamma = 1.76e11;     // 1/(Ts) [gyromagnetic ratio]
    
    double J = 0.0;             // 1.0e-21 J [coupling strength]
    J = J/mu_s;                 // T [renormalized coupling strength]
    
    double K1 = 3.0e4;//1.0e4; //1.53  // J/m^3 [uniaxial anisotropy constant]
    double K2 = .280e4; //.28
    double V = 3.05e-24;         // m^3 [particle volume]
    K1 = K1*V/mu_s;               // T [renormalized anisotropy]
    K2 = K2*V/mu_s;

    double Happ = 0.0;//2*K*V/(mu0*mu_s);        // A/m [external field]
    double B;                   // T [external field]
    B = mu0*Happ;
    
    double alpha = atof(argv[2]); //0.1;//0.1;         // LLG damping, dimensionless
    
    double sign = +1.0;         // to switch between the precession orientation

    double T = 300.0;             // K [temperature]
    double c, c1, c2; 
    c = 1.0/(1.0+alpha*alpha);
    c1 = c;
    c2 = c*sqrt(2.0*kB*T*alpha/mu_s);
 
    double D; 
    D = alpha*kB*T*mu_s/gamma;  // noise strength (square)
    
    // normalise the time (to computer friendly range)
    dt = dt*gamma;
    
    // basic info to set up 3D square lattice
    long num_neigh = 6;
    long num_part_side = 1;
    long num_part = num_part_side*num_part_side*num_part_side;

    // anisotropy vector - the same for all particles
    vector<double> n(3);
    n[0] = 1.0;
    n[1] = 0.0;
    n[2] = 0.0;
    normalise_vector(n);

    // field vector - the same for all particles
    vector<double> h(3);
    h[0] = 1.0;
    h[1] = 0.0;
    h[2] = 0.0;
    normalise_vector(h);

    // set the effective field vector 
    vector<double> h_eff(3);
    
    // spins - orient in the same direction
    vector<double> s(3, 0.0);
    vector<double> s_predict(3, 0.0);
    vector<double> sx(num_part, 1.0);
    vector<double> sy(num_part, 0.0);
    vector<double> sz(num_part, 0.0);
    vector<double> sx_predict(num_part);
    vector<double> sy_predict(num_part);
    vector<double> sz_predict(num_part);
    vector<double> mx(1);//(n_steps);
    vector<double> my(1);//(n_steps);
    vector<double> mz(1);//(n_steps);

    vector<double> a(3);
    vector<double> b(3);
    vector<double> a_predict(3);
    vector<double> b_predict(3);
    vector<double> ax(num_part);
    vector<double> ay(num_part);
    vector<double> az(num_part);
    vector<double> bx(num_part);
    vector<double> by(num_part);
    vector<double> bz(num_part);

    // random generator
    int seed = 1; //atoi(argv[2]);
    StochasticLib1 sto(seed);
    vector<double> wi(3);

    vector<double> wix(num_part);
    vector<double> wiy(num_part);
    vector<double> wiz(num_part);
   
    ofstream outfile0, outfile1;
    outfile0.open(argv[1]);//("magnetization_vs_time_3D.dat");

    /*********************************************************************/

    // set up spin lattice
    vector<long> neigh_list(num_part*num_neigh, 0);
    neigh_list = lattice3D_periodic(num_part_side, num_neigh);

    long lj;
    vector<double> param(7);
    param[0] = c1;
    param[1] = c2;
    param[2] = alpha;
    param[3] = dt;
    param[4] = J;
    param[5] = K1;
    param[6] = K2;
    param[7] = B;

    for(long i=0; i<n_steps; i++)
    {
        
        // keep track of time steps
        //if(i%1000 == 0.0) cout << "time: " << ((double) i)*dt/gamma << "\t" 
        //    << "tmax: " << ((double) n_steps)*dt/gamma << endl;

        // Solve the magnetization dynamics using the Heun integration scheme
        // (corrector-predictor type method). Requires renormalising the spin
        // magnitude in every step. 

        // Prediction step (Euler method)
        for(long j=0; j<num_part; j++)
        {
            wix[j] = sto.Normal(0.0, 1.0); // gauss random numbers
            wiy[j] = sto.Normal(0.0, 1.0);
            wiz[j] = sto.Normal(0.0, 1.0);

            wi[0] = wix[j];
            wi[1] = wiy[j];
            wi[2] = wiz[j];

            s[0] = sx[j];
            s[1] = sy[j];
            s[2] = sz[j];
            
            lj = neigh_list[num_neigh*j]; // pos. of j in the neigh-list
            h_eff = heff_eval_uni(sx,sy,sz,s,n,h,param,lj,num_neigh);
            //h_eff = heff_eval_cub(sx, sy, sz, s, h, param, lj, num_neigh);

            a = acoef(s, h_eff, param);
            b = bcoef(s, h_eff, wi, param);
            for(long l=0; l<3; l++)
            {
                s_predict[l] = s[l]+a[l]+b[l];
            }
            normalise_vector(s_predict); // renormalise the spin magnitude

            sx_predict[j] = s_predict[0];
            sy_predict[j] = s_predict[1];
            sz_predict[j] = s_predict[2];

            ax[j] = a[0];
            ay[j] = a[1];
            az[j] = a[2];

            bx[j] = b[0];
            by[j] = b[1];
            bz[j] = b[2];
        }

        // Correction step
        for(long j=0; j<num_part; j++)
        {
            wi[0] = wix[j];
            wi[1] = wiy[j];
            wi[2] = wiz[j];

            s[0] = sx[j];
            s[1] = sy[j];
            s[2] = sz[j];

            s_predict[0] = sx_predict[j];
            s_predict[1] = sy_predict[j];
            s_predict[2] = sz_predict[j];

            a[0] = ax[j];
            a[1] = ay[j];
            a[2] = az[j];

            b[0] = bx[j];
            b[1] = by[j];
            b[2] = bz[j];

            // correction
            h_eff = heff_eval_uni(sx,sy,sz,s_predict,n,h,param,lj,num_neigh);
            //h_eff = heff_eval_cub(sx, sy, sz, s_predict, h, param, lj, num_neigh);
            a_predict = acoef(s_predict, h_eff, param);
            b_predict = bcoef(s_predict, h_eff, wi, param);
            for(long l=0; l<3; l++)
            {
                s[l] += (a[l]+a_predict[l])/2.0 + (b[l]+b_predict[l])/2.0;
            }
            normalise_vector(s);

            sx[j] = s[0];
            sy[j] = s[1];
            sz[j] = s[2];
        }

        // calculate x, y, z components of magnetization
        for(long j=0; j<num_part; j++)
        {
            mx[0] += sx[j]/((double) num_part);
            my[0] += sy[j]/((double) num_part);
            mz[0] += sz[j]/((double) num_part);

            if(not(i%data_frequency)){
              outfile0 << ((double) i)*dt/gamma
                       << "\t" << mx[0]
                       << "\t" << my[0]
                       << "\t" << mz[0]
                       << "\t" << wi[0]
                       << "\t" << wi[1]
                       << "\t" << wi[2]
                       << endl;
            }

            mx[0] = 0.0;
            my[0] = 0.0;
            mz[0] = 0.0;
        }

        #if 0
        if(i>9000) // Here we export energy to check the boltzmann distribution
        {
            for(long j=0; j<num_part; j++)
                outfile1 << acos(h[0]*sx[j]+h[1]*sy[j]+h[2]*sz[j]) << endl;
        }
        #endif

    }

    /*for(long i=0; i<n_steps; i++)
        outfile0 << ((double) i)*dt/gamma 
                 << "\t" << mx[i]
                 << "\t" << my[i]
                 << "\t" << mz[i]
                 << "\t" << sqrt(mx[i]*mx[i]+my[i]*my[i]+mz[i]*mz[i])
                 << endl; */

     
    outfile0.close();

    return 0;
}
