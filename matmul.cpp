#include <iostream>
#include <stdlib.h>     // srand, rand
#include <time.h>       // time 
#include <chrono>
#include <ctime>
#include <math.h>       // sqrt 
#include <algorithm>    // min
using namespace std;

double dRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int matrix_multiply(int n, double* sequential, double* parallel)
{	
    int i, j, k;
    std::chrono::time_point<std::chrono::system_clock>  startTime, endTime;

    // Initialize random seed
    srand (time(NULL));

    // Initialize arrays
    double** A = new double*[n];
    for(int i = 0; i < n; ++i)
        A[i] = new double[n];
    
    double** B = new double*[n];
    for(int i = 0; i < n; ++i)
        B[i] = new double[n];

    double** C = new double*[n];
    for(int i = 0; i < n; ++i)
        C[i] = new double[n];

	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			A[i][j] =  dRand(0, 5);
            B[i][j] =  dRand(0, 5);
		}
	}

    //****************** Sequential *******************

    // Initializing elements of matrix C to 0.
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			C[i][j] = 0;
		}
	}

	// Multiplying matrix firstMatrix and secondMatrix and storing in array mult.
    startTime = std::chrono::system_clock::now();

	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			for(k=0; k<n; ++k)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}

    endTime = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = endTime-startTime; 
    *sequential = elapsed_seconds.count();

    //******************* Parallel *******************************

    // Initializing elements of matrix mult to 0.
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			C[i][j] = 0;
		}
	}

	// Multiplying matrix firstMatrix and secondMatrix and storing in array mult.
    startTime = std::chrono::system_clock::now();

    #pragma omp parallel for 
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			for(int k=0; k<n; ++k)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}

    endTime = std::chrono::system_clock::now();

    elapsed_seconds = endTime-startTime;
    *parallel = elapsed_seconds.count() ; 

    // clean memory
    for(int i = 0; i < n; ++i) {
        delete [] A[i];
        delete [] B[i];
        delete [] C[i];
    }
    delete [] A;
    delete [] B;
    delete [] C;

	return 0;
}

int getMinimumSamples(double* arr,int  samples){
    double mean, smean;
    double total =0;
    double stotal = 0;
    for (int i = 0; i < samples; i ++){
        total += arr[i];
    }
    mean = total/samples;
    for (int i = 0; i < samples; i ++){
        stotal +=  pow ((arr[i] - mean), 2.0);
    }
    smean = stotal/(samples-1);
    double z = sqrt(smean);
    
    double n = pow((100*1.96*z/(5*mean)), 2.0);
    // cout << "mean " << mean << " total " << total << " smean " << smean << " stotal " << stotal << " z " << z << " n " << n << endl;
    if(n<1){
        return 1;
    }
    return  (int)(ceil(n));
}

int main(){

    int testRuns = 20;
    int samples;
    double parallelTime = 0.0;
    double sequentialTime = 0.0;
    double totalSequenceTime;
    double totalParallelTime;
    double averageSequentialTime;
    double averageParallelTime;

    double* parallel = &parallelTime;
    double* sequential = &sequentialTime;

    
    for (int n =200; n<=2000;n+=200){

        // Test run number of 'testRuns' samples to determine minimum samples for ±5% accuracy and 95% confidence level
        double* samplesSequencial = new double[testRuns];
        double* samplesParallel = new double[testRuns];
        for(int j = 0; j<testRuns; j++){
            matrix_multiply(n, sequential, parallel);
            samplesSequencial[j] = *sequential;
            samplesParallel[j] = *parallel;
        }

        // Get minimum sample size for accuracy of ±5% and 95% confidence level for both sequential, parallel algos
        int samplesS  = getMinimumSamples(samplesSequencial, testRuns);
        int samplesP  = getMinimumSamples(samplesParallel, testRuns);

        samples = max(samplesS, samplesP); // Set minimum sample size as the maximum of above two
        cout << "n: "<< n<< " Minimum samples [±5% 95%] :" << samples << endl;

        // Do for number of 'samples' samples
        totalSequenceTime = 0.0;
        totalParallelTime = 0.0;
        for(int j = 0; j<samples; j++){
            cout << "iteration " << j << " of " << n << "\r";
            matrix_multiply(n, sequential, parallel);
            totalSequenceTime += *sequential;
            totalParallelTime += *parallel;
        }
        averageSequentialTime = totalSequenceTime/samples;
        averageParallelTime = totalParallelTime/samples;

        cout << "Average Sequential: " << averageSequentialTime << "\t";
        cout << "Average Parallel: " << averageParallelTime << "\n\n";
    }
}