 
#pragma once 
#include <vector>
//CPLUS_INCLUDE_PATH =~/projects 
//export CPLUS_INCLUDE_PATH
#include <math.h> 
using namespace std;

 

//LINEAR GAUSS MARKOV model
//or Hagan's numeraire

//a struct to initiata volatilty and mean reversion
struct lgmonefactor
{
    vector<double> volatility;
    double meanreversionparameter;
    //ctor
    //move-no need
    //
};

//functor to calculate mean reversion integral
struct lgmmeanreversion //int_0^t exp(-\int_0^s param du) ds
{
    double param;
    lgmmeanreversion(double p): param(p){}
    double operator() (const double& t) const //time
    {
        return (1-exp(-t*param))/param;
    }
};

//should take a model (volatility and mean rversion param)
//and some other data (strike etc.)
//and produce a swaption price
struct lgmswaptionprice
{
    double operator() (double vol, double strike, size_t numOfPeriods, 
                double meanreversionparam,
                double breakEvenRate,
                double startTime, //swap start
               const vector<double>& tau, //period year fractions
               const vector<double>& paymTimes, 
                 const vector<double>& discVector)
    {
        lgmmeanreversion H(meanreversionparam); 
        double aux = 0;
        //auto e = discVector.begin(); e != discVector.end(); ++e
        for (size_t i = 0; i<numOfPeriods; ++i)
        {
            double h = H(paymTimes[i]) - H(paymTimes[0]);
            if (i == 0) 
                aux += -discVector[i] * normalCdf(breakEvenRate/sqrt(vol));
            else
            {                
                aux += tau[i]*discVector[i] * strike
                        * normalCdf((breakEvenRate + vol* h)/sqrt(vol));
            }
            
            if  (i == numOfPeriods-1)
            {
                aux += discVector[i] 
                        * normalCdf((breakEvenRate + vol* h)/sqrt(vol));
            } 
        }
        return aux;
    }
};

class lgmpde //only explicit method based discrete version
{
    vector<double> payAtEnd;
    vector<double,double> upperLowerBoundary; //this should be zero for lgm

    public
};
