#include <iostream> // input output
#include <cmath> // (for sqrt)
#include <vector>
#include <tuple>
#include "readParameters.hpp"
#include "GetPot.hpp"
#include "gnuplot-iostream.hpp"// interface with gnuplot
/*!
  @file main.cpp
  @brief Temperature distribution in a 1D bar.

  @detail
    We solve  \f$ -T^{\prime\prime}(x)+act*(T(x)-T_e)=0, 0<x<L \f$ with 
    boundary conditions \f$ T(0)=To; T^\prime(L)=0\f$
    
    **************************************************
    Linear finite elements
    Resolution by Thomas algorithm.
    **************************************************
    
    Example adapted by Luca Formaggia from  a code found in 
    "Simulation numerique an C++" di I. Danaila, F. Hecht e
    O. Pironneau.
*/
//! helper function
void printHelp()
{
  std::cout<<"USAGE: main [-h] [-v] -p parameterFile (default: parameters.pot)"<<std::endl;
  std::cout<<"-h this help"<<std::endl;
  std::cout<<"-v verbose output"<<std::endl;
}

//! main program
int main(int argc, char** argv)
{
  using namespace std; // avoid std::
  int status(0); // final program status
  GetPot   cl(argc, argv);
  if( cl.search(2, "-h", "--help") )
    {
      printHelp();
      return 0;
    }
  // check if we want verbosity
  bool verbose=cl.search(1,"-v");
  // Get file with parameter values
  string filename = cl.follow("parameters.pot","-p");
  cout<<"Reading parameters from "<<filename<<std::endl;
  // read parameters
  const parameters param=readParameters(filename,verbose);
  // Transfer parameters to local variables
  // I use references to save memory (not really an issue here, it is just
  // to show a possible  use of references)
  // Here I use auto (remember that you need const and & if you want constant references)
  const auto& L= param.L;  // Bar length
  const auto& a1=param.a1; // First longitudinal dimension
  const auto& a2=param.a2; //  Second longitudinal dimension
  const auto& To=param.To; // Dirichlet condition
  const auto& Te=param.Te; // External temperature (Centigrades)
  const auto& k=param.k;  // Thermal conductivity
  const auto& hc=param.hc; // Convection coefficient
  const auto& M=param.M; // Number of grid elements
  const auto& result_name=param.result_name; // name of the output file
  
  //! Precomputed coefficient for adimensional form of equation
  const auto act=2.*(a1+a2)*hc*L*L/(k*a1*a2);

  // mesh size
  const auto h=1./M;
  
  // Solution vector
  std::vector<double> theta(M+1,0);
  
// matrix
   std::vector <double> c(M-1,-1);
   std::vector <double> b(M-1,-1);
   std::vector <double> a(M,1);
   for(int i=0;i<M-1;i++)
       a[i]= 2+ h*h*act;
 // termine noto
    std::vector <double> g(M,0);
    g[0]=(To-Te)/Te;
     
 // Thomas algorithm 
   
 // fattorizzazione
    for(int i=0;i<M-1;i++){
	b[i]=b[i]/a[i];
        a[i+1]=a[i+1]-b[i]*c[i];
     }
     
  // soluzione
     std::vector<double> y(M);
     y[0]=g[0];
     for (int i=1; i<M;i++)
          y[i]=g[i]-b[i-1]*y[i-1];

     theta[M]=y[M-1]/a[M-1];
     for(int i=1;i<M;i++)
       theta[M-i]=(y[M-i-1]-c[M-i-1]*theta[M-i+1])/a[M-i-1];
     theta[0]=(To-Te)/Te;
 // Analitic solution
    
    vector<double> thetaa(M+1,0);
     for(int m=0;m <= M;m++)
       thetaa[m]=Te+(To-Te)*cosh(sqrt(act)*(1-m*h))/cosh(sqrt(act));

     // writing results with format
     // x_i u_h(x_i) u(x_i) and lauch gnuplot 

     Gnuplot gp;
     std::vector<double> coor(M+1);
     std::vector<double> sol(M+1);
     std::vector<double> exact(M+1);

     cout<<"Result file:"<<result_name<<endl;
     ofstream f(result_name);
     for(int m = 0; m<= M; m++)
       {
	 // \t writes a tab 
         f<<m*h*L<<"\t"<<Te*(1.+theta[m])<<"\t"<<thetaa[m]<<endl;
	 // An example of use of tie and tuples!
         
	 std::tie(coor[m],sol[m],exact[m])=
	   std::make_tuple(m*h*L,Te*(1.+theta[m]),thetaa[m]);

       }
     // Using temporary files (another nice use of tie)
     gp<<"plot"<<gp.file1d(std::tie(coor,sol))<<
       "w lp title 'uh',"<< gp.file1d(std::tie(coor,exact))<<
       "w l title 'uex'"<<std::endl;
     f.close();
     return status;
}
