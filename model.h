 
#pragma once 
#include <vector>
//CPLUS_INCLUDE_PATH =~/projects 
//export CPLUS_INCLUDE_PATH
#include <math.h> 
#include "nonlinsolver/methods.h"
#include "marketdata/mdata.h"
using namespace std;




//LINEAR GAUSS MARKOV model
//or Hagan's numeraire

//a struct to initiata volatilty and mean reversion
struct lgmonefactor
{
    vector<double> vol; // volatility;
    double meanrev; // mean reversion parameter;
    lgmmeanreversion H;
    lgmonefactor (double expiry = 0.0, double m = 0.0):
        vol(int(expiry/0.25)), meanrev(m) , H(m)   {} 
        //give expiry in years, 
        //assume vol vector will be provided
        //for each starting three month period

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

struct swaption
{
    double swaptionExpiry; //in years = swap start time
    double swapTenor; //in years
    double strike;
    size_t numOfPeriods;
    vector<unsigned short> paymTimes;
}

//should take a model (volatility and mean rversion param)
//and some other data (strike etc.)
//and produce a swaption price
//THIS WILL BE A FUNCTOR in volatility
struct lgmswaptionprice
{
    swaption& s;
    timezero& z; 
    lgmonefactor& model; //model->meanrev
    lgmswapprice price;


    lgmswaptionprice() {}

    double operator() (double vol)  
    {             
        //there is composition here
        //so please see if it is better to write separate functions
        //and compose them, ala functional programming
        double breakEvenRate = bisectionMethod(price,-50,50,0.0000001);
         
        double aux = 0;
        //auto e = discVector.begin(); e != discVector.end(); ++e
        for (size_t i = 0; i<s.numOfPeriods; ++i)
        {
            double h = model.H(s.paymTimes[i]) - model.H(s.paymTimes[0]);
            if (i == 0) 
                aux += -z.discVector[i] * normalCdf(breakEvenRate/sqrt(vol));
            else
            {                
                aux += tau[i]*z.discVector[i] * s.strike
                        * normalCdf((breakEvenRate + vol* h)/sqrt(vol));
            }
            
            if  (i == s.numOfPeriods-1)
            {
                aux += z.discVector[i] 
                        * normalCdf((breakEvenRate + vol* h)/sqrt(vol));
            } 
        }
        return aux;
    }
};

/*
class lgmpde //only explicit method based discrete version
{
    vector<double> payAtEnd;
    vector<std::tuple<double,double>> upperLowerBoundary; //this should be zero for lgm
 
};
*/

//using formula 5.7a from Hagan's "Evaluating and hedging ..."
struct lgmswapprice
{
    lgmonefactor& model;  
    timezero& z;
    double vol;
    double fixrate;
    double swapstarttime;
     
    double operator()(double markov) const
    {
        size_t n = paymTimes.size();
        double value = 0.0;
        for (size_t i; i<n;++n)
        {
            double hdiff = model.H(paymTimes[i]) - model.H(swapstarttime);

            value += tau[i] * z.discVector[i]*exp(k())
        }

        return value*fixrate;
    }

};
