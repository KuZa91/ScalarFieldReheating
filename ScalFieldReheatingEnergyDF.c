//===================================================//
// Scalar Fielf Reheating                            //
//===================================================//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <random>
#include <fftw3.h>

#define REAL 0
#define IMAG 1

using namespace std;

double N_Norm_Squared(double *a, int64_t n); // Function to estimate the N dimensional norm squared of a vector a
double NRG(double ****Field, double ****Mom, double a, double lambda, double m, int64_t N, int64_t n_ord); // Function for estimating the energy of the field configuration 
double* AVG_Phi(double ****Field, double a,  int64_t N, int64_t n_ord); // Function that return the average value of the components of the field
double Omega_k(int i, int j, int k, int N, double a, double m); // Function for estimating the energy omega of the free theory in Fourier Space
fftw_complex* Diagonal(int N, fftw_complex *field);
double* F_k(int N, fftw_complex *diag);
double* derF_k(int N, double d_eta, fftw_complex *diag_old, fftw_complex *diag_now);


int main ()
{
	//***************************************
    // Defining the variables of the system *
    //***************************************

	int64_t opt =1, N, N_order = 4, N_print = 6, N_samples, step, i, j, k, l, app, app_2, samples, seed = 0;
    int dim[4];
    double lambda, m, Lenght, a, d_eta, eta, eta_max;
    int *print_val;
    bool logic = true, done;
    double ****Phi, ****pi;
    char fname[] = "FieldNP00eta0000.dat"; //modify fname[7], fname[8], fname[12], fname[13], fname[14], fname[15] to change file name

    //****************************************************
    // Initializing the seed and the normal distribution *
    //****************************************************

    mt19937 genrnd(seed);

    normal_distribution<> Gauss{0,1};

    uniform_real_distribution<double> Unif(0.0,1.0); // Generate a random number between 0 and 1, will be used to initialize randomically the spins to +-1
    
    //**********************************************************
    // Reading and setting the values of the working variables *
    //**********************************************************

    
    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;
    cout<< "This program will simulate the Reheating due to a O(N) scalar field"<< endl;
    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;
    
    do 
    {
        cout<< "Press 0 for a simulation with the standard parameters or 1 to input each parameter"<< endl;
        cin >> opt ;
        if((opt == 1)||(opt == 0))
        {
            logic = false;
        }
    }
    while (logic);

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do 
    {
        if(opt ==1)
        {
            cout<< "Please Insert the Lattice steps of the grid (>0 sugg 20)"<< endl;
            cin >> N ;
        }
        else
        {
            N = 20;
            cout<< "Lattice steps of the grid initialized to : "<<N<< endl;
        }
    }
    while (N <= 0);

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

     do 
    {
        if(opt ==1)
        {
            cout<< "Please Insert the Lattice size of the grid (>0. sugg 1.)"<< endl;
            cin >> Lenght ;
        }
        else
        {
            Lenght = 1.;
            cout<< "Lattice lenght initialized to : "<<Lenght<< endl;
        }
    }
    while (Lenght <= 0);

    a = Lenght/N;

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do 
    {
        if(opt ==1)
        {
            cout<< "Please Insert the number of samples of the field to generate during the simulation (>0 max 100)"<< endl;
            cin >> N_samples;
        }
        else
        {
            N_samples = 1;
            cout<< "Number of requested samples initialized to : "<<N_samples<< endl;
        }
    }
    while (N_samples <= 0 || N_samples > 100);

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do 
    {
        if(opt ==1)
        {
   		   cout<< "Please Insert the value of the coupling constant lambda (>0 sugg 0.1)"<< endl;
   		   cin >> lambda;
        }
        else
        {
            lambda = 0.1;
            cout<< "Coupling constant lambda initialized to : "<<lambda<< endl;
        }
    }
    while (lambda <= 0.); 

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do 
    {
        if(opt ==1)
        {
           cout<< "Please Insert the value of the mass of the field (>0 sugg 1.)"<< endl;
           cin >> m;
        }
        else
        {
            m = 1.;
            cout<< "Mass of the field initialized to : "<<m<< endl;
        }
    }
    while (m <= 0.); 

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do 
    {
        if(opt ==1)
        {
            cout<< "Please Insert the step delta_eta for evolution of the simulation in minkowski time (>0. sugg 0.01)"<< endl;
            cin >> d_eta;
        }
        else
        {
            d_eta = 0.01;
            cout<< "The evolution step delta_eta is initialized to : "<<d_eta<< endl;
        }
    }
    while (d_eta <= 0.); 

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    do 
    {
        if(opt ==1)
        {
            cout<< "Please Insert the final time for the evolution of the simulation in minkowski time (>0. max 10000.)"<< endl;
            cin >> eta_max;
        }
        else
        {
            eta_max = 1000.;
            cout<< "The final time for the evolution of the simulation is initialized to : "<<eta_max<< endl;
        }
    }
    while (eta_max <= 0. || eta_max > 10000.); 

    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;

    //*********************************************************************************************************************************
    // Initializing the Correlation vectors, we need one for the previous and one for the next step to estimate the second derivative *
    //*********************************************************************************************************************************

    double *F_0, *F_1, *F_2, *F_3;
    double *dd_corr_0, *dd_corr_1, *dd_corr_2, *dd_corr_3;
    double *n_0, *n_1, *n_2, *n_3;
    double *phi_avg;

    phi_avg = new double[N_order];

    F_0 = new double[(N/2 + 1)];
    F_1 = new double[(N/2 + 1)];
    F_2 = new double[(N/2 + 1)];
    F_3 = new double[(N/2 + 1)];

    dd_corr_0 = new double[(N/2 + 1)];
    dd_corr_1 = new double[(N/2 + 1)];
    dd_corr_2 = new double[(N/2 + 1)];
    dd_corr_3 = new double[(N/2 + 1)];

    n_0 = new double[(N/2 + 1)];
    n_1 = new double[(N/2 + 1)];
    n_2 = new double[(N/2 + 1)];
    n_3 = new double[(N/2 + 1)];

    //*****************************************
    // Initializing the starting lattice grid *
    //*****************************************

    Phi = new double***[N];
    pi = new double***[N];

    for(i = 0; i<N; i++)
    {
        Phi[i] = new double**[N];
        pi[i] = new double**[N];

        for(j=0; j<N; j ++)
        {
            Phi[i][j] = new double*[N];
            pi[i][j] = new double*[N];

            for(k = 0; k<N; k++)
            {
                Phi[i][j][k] = new double[N_order];
                pi[i][j][k] = new double[N_order];

                for(l = 0; l < N_order; l++)
                {
                    //******************************************************
                    // Initialization as a gaussian for the starting field *
                    //******************************************************

                    Phi[i][j][k][l] = Gauss(genrnd);

                    //**********************************************
                    // Initialization to 0 for the initial momenta *
                    //**********************************************

                    pi[i][j][k][l] = 0.;

                }
            }

        }
    }

    //**************************************************************************************************************************
    // Opening the file in which the data will be stored and the array that will choose the saving positions of the simulations*
    //**************************************************************************************************************************

    ofstream outputfile;

    outputfile.open("Energy.dat", ofstream::out /*| ofstream::app*/); // use "out" only for writing new file and "app" to append

    print_val = new int[N_print];

    print_val[0] = int(eta_max/(100.*d_eta)); // saving simulation values at 1% of simulation
    print_val[1] = int(eta_max/(20.*d_eta)); // saving simulation values at 5% of simulation
    print_val[2] = int(eta_max/(10.*d_eta)); // saving simulation values at 10% of simulation
    print_val[3] = int(eta_max/(5.*d_eta)); // saving simulation values at 20% of simulation
    print_val[4] = int(eta_max/(2.*d_eta)); // saving simulation values at 50% of simulation
    print_val[5] = int(eta_max/d_eta); // saving simulation values at the end of simulation



    //******************************************************************
    // Allocating the memory and defining variables for the field fftw *
    //******************************************************************

    fftw_complex *field_fft_0, *field_fft_1, *field_fft_2, *field_fft_3;
    fftw_complex *result_0, *result_1, *result_2, *result_3;
    fftw_complex *diag_old_0, *diag_old_1, *diag_old_2, *diag_old_3;
    fftw_complex *diag_now_0, *diag_now_1, *diag_now_2, *diag_now_3;

    field_fft_0 = (fftw_complex*) fftw_malloc( N * N * N * sizeof(fftw_complex));
    field_fft_1 = (fftw_complex*) fftw_malloc( N * N * N * sizeof(fftw_complex));
    field_fft_2 = (fftw_complex*) fftw_malloc( N * N * N * sizeof(fftw_complex));
    field_fft_3 = (fftw_complex*) fftw_malloc( N * N * N * sizeof(fftw_complex));

    result_0 = (fftw_complex*) fftw_malloc( N * N * N * sizeof(fftw_complex)); 
    result_1 = (fftw_complex*) fftw_malloc( N * N * N * sizeof(fftw_complex)); 
    result_2 = (fftw_complex*) fftw_malloc( N * N * N * sizeof(fftw_complex)); 
    result_3 = (fftw_complex*) fftw_malloc( N * N * N * sizeof(fftw_complex)); 

    diag_old_0 = (fftw_complex*) fftw_malloc( N * sizeof(fftw_complex));
    diag_old_1 = (fftw_complex*) fftw_malloc( N * sizeof(fftw_complex));
    diag_old_2 = (fftw_complex*) fftw_malloc( N * sizeof(fftw_complex));
    diag_old_3 = (fftw_complex*) fftw_malloc( N * sizeof(fftw_complex));

    diag_now_0 = (fftw_complex*) fftw_malloc( N * sizeof(fftw_complex));
    diag_now_1 = (fftw_complex*) fftw_malloc( N * sizeof(fftw_complex));
    diag_now_2 = (fftw_complex*) fftw_malloc( N * sizeof(fftw_complex));
    diag_now_3 = (fftw_complex*) fftw_malloc( N * sizeof(fftw_complex));


    //******************************************************************
    // Starting the proper body of the simulation *
    //******************************************************************

    eta = 0.;

    cout<< "The simulation may now start !"<< endl;

    for (samples = 0; samples < N_samples; samples ++)
    {

        //*********************************************
        // Initializing the samples field and momenta *
        //*********************************************

        app = 0;
        step = 0;

        for(i = 0; i<N; i++)
        {
    
            for(j=0; j<N; j ++)
            {
        
                for(k = 0; k<N; k++)
                {
                    
                    for(l = 0; l < N_order; l++)
                    {   
                        //******************************************************
                        // Initialization as a gaussian for the starting field *
                        //******************************************************

                        Phi[i][j][k][l] = Gauss(genrnd);                        

                        //**********************************************
                        // Initialization to 0 for the initial momenta *
                        //**********************************************

                        pi[i][j][k][l] = 0.;


                        //************************************************************************************************************
                        // Initializing the fftw field vector in row-major order, Phi_{i,j,k,l} -> Arr[ l + N_ord*(k + N*(j + N*i))] *
                        //************************************************************************************************************

                        switch (l) 
                        {
                            case 0:
                                    field_fft_0[app][REAL] = Phi[i][j][k][l];
                                    field_fft_0[app][IMAG] = 0.;  
                                    break;

                            case 1:
                                    field_fft_1[app][REAL] = Phi[i][j][k][l];
                                    field_fft_1[app][IMAG] = 0.;  
                                    break;

                            case 2:
                                    field_fft_2[app][REAL] = Phi[i][j][k][l];
                                    field_fft_2[app][IMAG] = 0.;  
                                    break;

                            case 3:
                                    field_fft_3[app][REAL] = Phi[i][j][k][l];
                                    field_fft_3[app][IMAG] = 0.;  
                                    break;            
                        }
                    }
                    app ++;
                }
            }
        }


        //****************************************************************************************************************
        // Estimating the FFTW of the field in order to estimate the initial condition of it by implementing free theory *
        //****************************************************************************************************************

        
        for(l = 0; l < N_order; l++)
        {
            switch (l) 
                        {
                            case 0: {
                                    fftw_plan FFT0 = fftw_plan_dft_3d(N, N, N, field_fft_0,
                                                     result_0, FFTW_FORWARD, FFTW_ESTIMATE);
                                    fftw_execute(FFT0);
                                    fftw_destroy_plan(FFT0);
                                    break;
                                    }

                            case 1: {
                                    fftw_plan FFT1 = fftw_plan_dft_3d(N, N, N, field_fft_1,
                                                     result_1, FFTW_FORWARD, FFTW_ESTIMATE);
                                    fftw_execute(FFT1);
                                    fftw_destroy_plan(FFT1);
                                    break;
                                    }

                            case 2: {
                                    fftw_plan FFT2 = fftw_plan_dft_3d(N, N, N, field_fft_2,
                                                     result_2, FFTW_FORWARD, FFTW_ESTIMATE);
                                    fftw_execute(FFT2);
                                    fftw_destroy_plan(FFT2);
                                    break;
                                    }

                            case 3: {
                                    fftw_plan FFT3 = fftw_plan_dft_3d(N, N, N, field_fft_3,
                                                     result_3, FFTW_FORWARD, FFTW_ESTIMATE);
                                    fftw_execute(FFT3);
                                    fftw_destroy_plan(FFT3);
                                    break;  
                                    }         
                        }   
        }

        //*****************************************************************************
        // Implementing free theory initial conditions in the void case, where n_p = 0*
        //*****************************************************************************

        app = 0;

        for(i = 0; i<N; i++)
        {
    
            for(j=0; j<N; j ++)
            {
        
                for(k = 0; k<N; k++)
                {
                        
                    
                        result_0[app][REAL] *= 1./sqrt(2.*Omega_k(i, j, k, N, a, m));
                        result_0[app][IMAG] *= 1./sqrt(2.*Omega_k(i, j, k, N, a, m));

                        result_1[app][REAL] *= 1./sqrt(2.*Omega_k(i, j, k, N, a, m));
                        result_1[app][IMAG] *= 1./sqrt(2.*Omega_k(i, j, k, N, a, m));

                        result_2[app][REAL] *= 1./sqrt(2.*Omega_k(i, j, k, N, a, m));
                        result_2[app][IMAG] *= 1./sqrt(2.*Omega_k(i, j, k, N, a, m));

                        result_3[app][REAL] *= 1./sqrt(2.*Omega_k(i, j, k, N, a, m));
                        result_3[app][IMAG] *= 1./sqrt(2.*Omega_k(i, j, k, N, a, m));

                        app ++;
                }
            }
        }  

        //****************************************************
        // Going back from Fourier space to Coordinate space *
        //****************************************************

        for(l = 0; l < N_order; l++)
        {
            switch (l) 
                        {
                            case 0: {

                                    fftw_plan IFFT0 =  fftw_plan_dft_3d(N, N, N, result_0,
                                                       field_fft_0, FFTW_BACKWARD, FFTW_ESTIMATE);
                                    fftw_execute(IFFT0);
                                    fftw_destroy_plan(IFFT0);
                                    break;
                                    }

                            case 1: {
                                    fftw_plan IFFT1 = fftw_plan_dft_3d(N, N, N, result_1,
                                                      field_fft_1, FFTW_BACKWARD, FFTW_ESTIMATE);
                                    fftw_execute(IFFT1);
                                    fftw_destroy_plan(IFFT1);
                                    break;
                                    }

                            case 2: {
                                    fftw_plan IFFT2 = fftw_plan_dft_3d(N, N, N, result_2,
                                                      field_fft_2, FFTW_BACKWARD, FFTW_ESTIMATE);
                                    fftw_execute(IFFT2);
                                    fftw_destroy_plan(IFFT2);
                                    break;
                                    }

                            case 3: {
                                    fftw_plan IFFT3 = fftw_plan_dft_3d(N, N, N, result_3,
                                                      field_fft_3, FFTW_BACKWARD, FFTW_ESTIMATE);
                                    fftw_execute(IFFT3);
                                    fftw_destroy_plan(IFFT3);
                                    break;
                                    }        
                        }   
        }

        app = 0;

        for(i = 0; i<N; i++)
        {
    
            for(j=0; j<N; j ++)
            {
        
                for(k = 0; k<N; k++)
                {
                    
                    for(l = 0; l < N_order; l++)
                    {   
                        
                        switch (l) 
                        {
                            case 0:
                                {    
                                    Phi[i][j][k][l] = (1./lambda) + (field_fft_0[app][REAL]/pow(N,3.)); //We need the N^3 term because the inverse fftw is not normalized, we add the expectation value on the 0 component                           
                                    field_fft_0[app][IMAG] = 0.;
                                    break;
                                }

                            case 1:
                                {
                                    Phi[i][j][k][l] = field_fft_1[app][REAL]/pow(N,3.);
                                    field_fft_1[app][IMAG] = 0.;
                                    break;
                                }

                            case 2:
                                {
                                    Phi[i][j][k][l] = field_fft_2[app][REAL]/pow(N,3.);
                                    field_fft_2[app][IMAG] = 0.;
                                    break;
                                }

                            case 3:
                                {
                                    Phi[i][j][k][l] = field_fft_3[app][REAL]/pow(N,3.);
                                    field_fft_3[app][IMAG] = 0.;
                                    break;
                                }       
                        }
                        app ++;

                    }
                }
            }
        }        
        
        //***************************************************************************************
        // Estimating the step at d_eta/2 in order to implement the Leap Frog Symplectic solver *
        //***************************************************************************************

        if(eta == 0.)
        {
            for(i = 0; i<N; i++)
            {
        
                for(j=0; j<N; j ++)
                {
                

                    for(k = 0; k<N; k++)
                    {
                    

                        for(l = 0; l < N_order; l++)
                        {   
                            //**********************************************************************************
                            // Evolution for half the step of the momenta in order to run the leap frog solver *
                            //**********************************************************************************

                            pi[i][j][k][l] += (d_eta/(2.*a*a))*((Phi[(i + 1)%N][j][k][l] - 2.*Phi[i][j][k][l] + Phi[(i -1 + N)%N][j][k][l]) + 
                                            (Phi[i][(j +1)%N][k][l] - 2.*Phi[i][j][k][l] + Phi[i][(j -1 + N)%N][k][l]) + 
                                            (Phi[i][j][(k +1)%N][l] - 2.*Phi[i][j][k][l] + Phi[i][j][(k -1 + N)%N][l]) - 
                                            (a*a*lambda/6.)*N_Norm_Squared(Phi[i][j][k], N_order)*Phi[i][j][k][l] - pow(m*a,2.)*Phi[i][j][k][l]) ;

                        }
                    }

                }
            }
        }

        //****************************************************
        // Time evolution of the field from eta 0 to eta max *
        //****************************************************

        do
        {
            for(i = 0; i<N; i++)
            {
        
                for(j=0; j<N; j ++)
                {
                
                    for(k = 0; k<N; k++)
                    {
                    
                        for(l = 0; l < N_order; l++)
                        {   
                            Phi[i][j][k][l] += d_eta*pi[i][j][k][l];
                        }
                    }

                }
            }

            for(i = 0; i<N; i++)
            {
        
                for(j=0; j<N; j ++)
                {
                
                    for(k = 0; k<N; k++)
                    {
                    
                        for(l = 0; l < N_order; l++)
                        {
                            pi[i][j][k][l] += (d_eta/(a*a))*((Phi[(i + 1)%N][j][k][l] - 2.*Phi[i][j][k][l] + Phi[(i -1 + N)%N][j][k][l]) + 
                                              (Phi[i][(j +1)%N][k][l] - 2.*Phi[i][j][k][l] + Phi[i][(j -1 + N)%N][k][l]) + 
                                              (Phi[i][j][(k +1)%N][l] - 2.*Phi[i][j][k][l] + Phi[i][j][(k -1 + N)%N][l]) - 
                                              (a*a*lambda/6.)*N_Norm_Squared(Phi[i][j][k], N_order)*Phi[i][j][k][l]
                                              - pow(m*a,2.)*Phi[i][j][k][l]) ;

                        }
                    }
                }
            }   
            

            eta += d_eta;
            step ++;

            phi_avg = AVG_Phi(Phi, a, N, N_order);

            if( outputfile.is_open() )
            {   
                
                outputfile<< setprecision(20)  // precision of following floating point output
                << setfill(' ')      // character used to fill the column
                << left              // left alignment  -- remember all comes to left of affected object here.
                << scientific
                << setw(30) << eta << setw(30) << NRG(Phi, pi, a, lambda, m, N, N_order) << setw(30)
                << phi_avg[0] << setw(30) << phi_avg[1] << setw(30) << phi_avg[2] << setw(30) << phi_avg[3] 
                << endl;
                

            } 
            else 
            {
                cerr << "Cannot write to file ! " << endl; //cerr use no buffer, direct message regardless of stack turn
            }


            //********************************************************************************************
            // Checking if the value of eta is in the list of values where we need to estimate the np(k) *
            //********************************************************************************************

            for(app_2 = 0; app_2 < N_print; app_2 ++)
            {
                break ;
                if((step >= print_val[app_2] - 1 )&&(step <= print_val[app_2]))
                {
                    done = false;
                    if(step < print_val[app_2])
                    {
                        app = 0;
                        for(i = 0; i<N; i++)
                        {
    
                            for(j=0; j<N; j ++)
                            {
        
                                for(k = 0; k<N; k++)
                                {
                    
                                    for(l = 0; l < N_order; l++)
                                    {   
                                        
                                        //************************************************************************************************************
                                        // Initializing the fftw field vector in row-major order, Phi_{i,j,k,l} -> Arr[ l + N_ord*(k + N*(j + N*i))] *
                                        //************************************************************************************************************

                                        switch (l) 
                                        {
                                            case 0:
                                                {
                                                    field_fft_0[app][REAL] = Phi[i][j][k][l];
                                                    field_fft_0[app][IMAG] = 0.;  
                                                    break;
                                                }

                                            case 1:
                                                {
                                                    field_fft_1[app][REAL] = Phi[i][j][k][l];
                                                    field_fft_1[app][IMAG] = 0.;  
                                                    break;
                                                }

                                            case 2:
                                                {
                                                    field_fft_2[app][REAL] = Phi[i][j][k][l];
                                                    field_fft_2[app][IMAG] = 0.;  
                                                    break;
                                                }

                                            case 3:
                                                {
                                                    field_fft_3[app][REAL] = Phi[i][j][k][l];
                                                    field_fft_3[app][IMAG] = 0.;  
                                                    break;
                                                }            
                                        }
                                    }
                                    app ++;
                                }
                            }
                        }

                        for(l = 0; l < N_order; l++)
                        {
                            switch (l) 
                            {
                                case 0: {
                                            fftw_plan FFT0 = fftw_plan_dft_3d(N, N, N, field_fft_0,
                                                     result_0, FFTW_FORWARD, FFTW_ESTIMATE);
                                            fftw_execute(FFT0);
                                            diag_old_0 = Diagonal(N, result_0);
                                            fftw_destroy_plan(FFT0);
                                            break;
                                        }

                                case 1: {
                                            fftw_plan FFT1 = fftw_plan_dft_3d(N, N, N, field_fft_1,
                                                     result_1, FFTW_FORWARD, FFTW_ESTIMATE);
                                            fftw_execute(FFT1);
                                            diag_old_1 = Diagonal(N, result_1);
                                            fftw_destroy_plan(FFT1);
                                            break;
                                        }

                                case 2: {
                                            fftw_plan FFT2 = fftw_plan_dft_3d(N, N, N, field_fft_2,
                                                     result_2, FFTW_FORWARD, FFTW_ESTIMATE);
                                            fftw_execute(FFT2);
                                            diag_old_2 = Diagonal(N, result_2);
                                            fftw_destroy_plan(FFT2);
                                            break;
                                        }

                                case 3: {
                                            fftw_plan FFT3 = fftw_plan_dft_3d(N, N, N, field_fft_3,
                                                     result_3, FFTW_FORWARD, FFTW_ESTIMATE);
                                            fftw_execute(FFT3);
                                            diag_old_3 = Diagonal(N, result_3);
                                            fftw_destroy_plan(FFT3);
                                            break;  
                                        }         
                            }   
                        }
                    }
                    else
                    {
                        app = 0;
                        for(i = 0; i<N; i++)
                        {
            
                            for(j=0; j<N; j ++)
                            {
                
                                for(k = 0; k<N; k++)
                                {
                            
                                    for(l = 0; l < N_order; l++)
                                    {   
                                                
                                        //************************************************************************************************************
                                        // Initializing the fftw field vector in row-major order, Phi_{i,j,k,l} -> Arr[ l + N_ord*(k + N*(j + N*i))] *
                                        //************************************************************************************************************

                                        switch (l) 
                                        {
                                            case 0:
                                                    {
                                                        field_fft_0[app][REAL] = Phi[i][j][k][l];
                                                        field_fft_0[app][IMAG] = 0.;  
                                                        break;
                                                    }

                                            case 1:
                                                    {
                                                        field_fft_1[app][REAL] = Phi[i][j][k][l];
                                                        field_fft_1[app][IMAG] = 0.;  
                                                        break;
                                                    }

                                            case 2:
                                                    {
                                                        field_fft_2[app][REAL] = Phi[i][j][k][l];
                                                        field_fft_2[app][IMAG] = 0.;  
                                                        break;
                                                    }

                                            case 3:
                                                    {
                                                        field_fft_3[app][REAL] = Phi[i][j][k][l];
                                                        field_fft_3[app][IMAG] = 0.;  
                                                        break;
                                                    }            
                                            
                                        }
                                        app ++;
                                    }
                                }
                            }
                        }

                        for(l = 0; l < N_order; l++)
                        {
                            switch (l) 
                            {
                                    case 0: {
                                                fftw_plan FFT0 = fftw_plan_dft_3d(N, N, N, field_fft_0,
                                                        result_0, FFTW_FORWARD, FFTW_ESTIMATE);
                                                fftw_execute(FFT0);
                                                diag_now_0 = Diagonal(N, result_0);
                                                F_0 = F_k(N, diag_now_0);
                                                dd_corr_0 = derF_k(N, d_eta, diag_old_0, diag_now_0);
                                                for(app = 0; app < N/2 + 1; app ++)
                                                {
                                                    n_0[app] = F_0[app]*sqrt(abs(dd_corr_0[app]/F_0[app])) - 0.5;
                                                }
                                                fftw_destroy_plan(FFT0);
                                                break;
                                            }

                                    case 1: {
                                                fftw_plan FFT1 = fftw_plan_dft_3d(N, N, N, field_fft_1,
                                                        result_1, FFTW_FORWARD, FFTW_ESTIMATE);
                                                fftw_execute(FFT1);
                                                diag_now_1 = Diagonal(N, result_1);
                                                F_1 = F_k(N, diag_now_1);
                                                dd_corr_1 = derF_k(N, d_eta, diag_old_1, diag_now_1);
                                                for(app = 0; app < N/2 + 1; app ++)
                                                {
                                                    n_1[app] = F_1[app]*sqrt(abs(dd_corr_1[app]/F_1[app])) - 0.5;
                                                }
                                                fftw_destroy_plan(FFT1);
                                                break;
                                            }

                                    case 2: {
                                                fftw_plan FFT2 = fftw_plan_dft_3d(N, N, N, field_fft_2,
                                                            result_2, FFTW_FORWARD, FFTW_ESTIMATE);
                                                fftw_execute(FFT2);
                                                diag_now_2 = Diagonal(N, result_2);
                                                F_2 = F_k(N, diag_now_2);
                                                dd_corr_2 = derF_k(N, d_eta, diag_old_2, diag_now_2);
                                                for(app = 0; app < N/2 + 1; app ++)
                                                {
                                                    n_2[app] = F_2[app]*sqrt(abs(dd_corr_2[app]/F_2[app])) - 0.5;
                                                }
                                                fftw_destroy_plan(FFT2);
                                                break;
                                            }

                                    case 3: {
                                                fftw_plan FFT3 = fftw_plan_dft_3d(N, N, N, field_fft_3,
                                                        result_3, FFTW_FORWARD, FFTW_ESTIMATE);
                                                fftw_execute(FFT3);
                                                diag_now_3 = Diagonal(N, result_3);
                                                F_3 = F_k(N, diag_now_3);
                                                dd_corr_3 = derF_k(N, d_eta, diag_old_3, diag_now_3);
                                                for(app = 0; app < N/2 + 1; app ++)
                                                {
                                                    n_3[app] = F_3[app]*sqrt(abs(dd_corr_3[app]/F_3[app])) - 0.5;
                                                }
                                                fftw_destroy_plan(FFT3);
                                                break;
                                            }  
                            }         
                        }  
                        done = true;                     
                        
                    }

                    //*****************************************************************************************************************************************
                    // Save the values of occupation numbers to file, the name of the output file will be automatically setted by the parameters of the cycle *
                    //*****************************************************************************************************************************************

                    if(done)
                    {
                        fname[7] = int(samples/10) + '0';
                        fname[8] = int(samples%10) + '0';

                        if(int(print_val[i]*d_eta) == 10000)
                        {
                            fname[12] = '9';
                            fname[13] = '9';
                            fname[14] = '9';
                            fname[15] = '9';
                        }
                        else
                        {
                            fname[12] = int((print_val[app_2]*d_eta)/1000.) + '0';
                            fname[13] = (int(print_val[app_2]*d_eta)%1000)/100 + '0';
                            fname[14] = (int(print_val[app_2]*d_eta)%100)/10 + '0';
                            fname[15] = int((print_val[app_2]*d_eta))%10 + '0';
                        }

                        outputfile.open(fname, ofstream::out /*| ofstream::app*/); // use "out" only for writing new file and "app" to append

                        if( outputfile.is_open() )
                        {   
                            for(app = 0; app < N/2 + 1; app ++)
                            {  
                                outputfile<< setprecision(20)  // precision of following floating point output
                                << setfill(' ')      // character used to fill the column
                                << left              // left alignment  -- remember all comes to left of affected object here.
                                << scientific
                                << setw(30) << ((M_PI *sqrt(3))/(a))*(1. - double(abs(app - N/2))/(N/2)) << setw(30) << n_0[app] << setw(30)
                                << n_1[app] << setw(30) << n_2[app] << setw(30) << n_3[app] << setw(30) << (1./3.)*(n_1[app] + n_2[app] + n_3[app]) 
                                << endl;
                            }

                        } 
                        else 
                        {
                            cerr << "Cannot write to file ! " << endl; //cerr use no buffer, direct message regardless of stack turn
                        }
                    }
                    outputfile.close();
                    
                }
            }

            //************************************************************************
            // Print the energy for the sample at a certain eta during the evolution *
            //************************************************************************

            if(step%print_val[2] == 0)
            {
                cout<<setprecision(10)  // precision of following floating point output
                << setfill(' ')      // character used to fill the column
                << left              // left alignment  -- remember all comes to left of affected object here.
                << scientific
                <<"We are at eta "<< eta <<" for sample N."<< samples<<" and the energy of the configuration is "
                << NRG(Phi, pi, a, lambda, m, N, N_order) << " !"<<endl;
            }
            
        }while(eta < eta_max);

        eta = 0.;
        outputfile.close();

    }
    cout<< "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"<< endl;
}

//*******************************************************
// Definition of the functions needed in the simulation *
//*******************************************************

double N_Norm_Squared(double *a, int64_t n) 
{
    //*************************************************************
    // Return the N dimensional norm squared of the  given vector *
    //*************************************************************

    int64_t i;

    double result = 0.;
    
    for (i = 0; i < n; i++) 
    {
        result += a[i]*a[i];
    }
    return result;
}

double NRG(double ****Field, double ****Mom, double a, double lambda, double m, int64_t N, int64_t n_ord)
{
    //**************************************************************************************************************************************************
    // This function will estimate the energy of the field configuration, in order to check that it would remain almost constant during time evolution *
    //**************************************************************************************************************************************************

    int64_t i,j,k,l;

    double result = 0.;

    for(i = 0; i<N; i++)
    {
        
        for(j=0; j<N; j ++)
        {
                

            for(k = 0; k<N; k++)
            {
                result += 0.5*N_Norm_Squared(Mom[i][j][k], n_ord)  + (lambda/24.)*pow(N_Norm_Squared(Field[i][j][k], n_ord),2.) + 
                m*m*N_Norm_Squared(Field[i][j][k], n_ord)/2.;

                for(l = 0; l < n_ord; l++)
                {   

                    result += 0.5*(pow((Field[(i + 1)%N][j][k][l] - Field[i][j][k][l])/(a),2.) + pow((Field[i][(j +1)%N][k][l] - Field[i][j][k][l])/a,2.)
                              + pow((Field[i][j][(k + 1)%N][l] - Field[i][j][k][l])/(a),2.));

                }
            }

        }
    }

    return pow(a,3.)*result;

}

double* AVG_Phi(double ****Field, double a,  int64_t N, int64_t n_ord)
{
    //*************************************************************************************************
    // This function return the average value of all the field components given a field configuration *
    //*************************************************************************************************

    int64_t i,j,k,l;
    double *result;

    result = new double[n_ord];

    for(l = 0; l < n_ord; l++)
    {
        result[l] = 0.;
    }

    for(i = 0; i<N; i++)
    {
        
        for(j=0; j<N; j ++)
        {
                

            for(k = 0; k<N; k++)
            {

                for(l = 0; l < n_ord; l++)
                {   

                    result[l] += pow(a,3.)*Field[i][j][k][l];

                }
            }

        }
    }

    return result;

}


double Omega_k(int i, int j, int k, int N, double a, double m)
{
    //**************************************************************************************************************************
    // This function will estimate the energy omega in fourier space for a free theory using the formula omega = sqrt(k² + m²) *
    //**************************************************************************************************************************

    double result = 0.;
    result += pow((M_PI/(a))*(1 - double(abs(i - N/2))/(N/2)),2.) + pow((M_PI/(a))*(1 - double(abs(j - N/2))/(N/2)),2.) + 
    pow((M_PI/(a))*(1 - double(abs(k - N/2))/(N/2)),2.) + m*m;
    return sqrt(result);
}

double* F_k(int N, fftw_complex *diag)
{
    //************************************************************************************************************
    // This function will return the 2 point function F_k in Fourier space given the diagonal of a 3d fftw field *
    //************************************************************************************************************

    int count;
    double *Corr_func;

    Corr_func = new double[N/2 + 1];

    for(count = 0; count < N/2 + 1; count ++)
    {
        Corr_func[count] = diag[count][REAL]*diag[(N-count)%N][REAL] - diag[count][IMAG]*diag[(N-count)%N][IMAG]; // <field(k) field(-k)> = Real(field(k)*field(-k))
    }
    return Corr_func;
}

double* derF_k(int N, double d_eta, fftw_complex *diag_old, fftw_complex *diag_now)
{
    //*****************************************************************************************************************
    // This function will return the derivative of the F functions given the diagonal at the actual and previous step *
    //*****************************************************************************************************************

    int count;
    double *der;

    der = new double[N/2 + 1];

    for(count = 0; count < N/2 + 1; count ++)
    {
        der[count] = (1./(d_eta*d_eta))*(diag_now[count][REAL] - diag_old[count][REAL])*(diag_now[(N-count)%N][REAL] - diag_old[(N-count)%N][REAL]);
    }
    return der;
}

fftw_complex* Diagonal(int N, fftw_complex *field)
{
    //*******************************************************************************
    // This function will return the diagonal of a 3d fftw vector in row major order*
    //*******************************************************************************

    int count;
    fftw_complex *diag;
    diag = (fftw_complex*) fftw_malloc( N* sizeof(fftw_complex));

    for(count = 0; count < N; count ++)
    {
        diag[count][REAL] = (1./pow(N,3./2.))*field[count + N*(count + N*count)][REAL]; //in row major order Mat[i][j][k] - > Arr[ k + N*(j + N*i)]
        diag[count][IMAG] = (1./pow(N,3./2.))*field[count + N*(count + N*count)][IMAG]; //in row major order Mat[i][j][k] - > Arr[ k + N*(j + N*i)]
    }
    return diag;
}