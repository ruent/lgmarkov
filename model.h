 
#pragma once 
#include <vector>
//CPLUS_INCLUDE_PATH =~/projects 
//export CPLUS_INCLUDE_PATH
#include <math.h> 
#include "nonlinsolver/methods.h"
#include "marketdata/mdata.h"
using namespace std;


const int mdsep = 1/12; //market data is for monthly increments

//LINEAR GAUSS MARKOV model
//or Hagan's numeraire

struct lgmonefactor
{
    //a struct to initiata volatilty and mean reversion
     
    // let the vol vector include values for all possible 1 months
    // then you can reduce this to the available swt of swaption expiries
    //this can make vol lookup easier
    vector<double> vol; 
    double meanrev; // mean reversion parameter;
    lgmmeanreversion H;
    lgmonefactor (double expiry = 0.0, double m = 0.0):
        vol(int(expiry/0.25)), meanrev(m) , H(m)   {} 
        //give expiry in years, 
        //assume vol vector will be provided
        //for each starting three month period

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
    double expiry; //in years = swap start time
    //this too is a multiple of 1months as of time 0
    double tenor; //in years
    double strike;
    size_t numOfPeriods;
    vector<unsigned int> paymTimes; 
    //paymTimes are as in marketdata convention:
    //multiples of 1 months as of time 0
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


    lgmswaptionprice(swaption _s, timezero _z, lgmonefactor _m): 
            s(_s), z(_z), m(_m)
      {
          price()
      }

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


struct lgmswapprice
{
    //using formula 5.7a from Hagan's "Evaluating and hedging ..."
    //swap is from a swaption
    //no need  for an independent swap, I hope!  
    lgmonefactor& m;  //mean reversion function H and vol data
    timezero& z; //cached discount factors
    swaption& s;
    
    //this is not really the swap price
    //but the one divided by exp(-H_0 y*-0.5 H_0^2 v) 
    double operator()(double x) const //x is the X ~N(0,\alpha_t)
    {         
        double value = 0.0;
        size_t j = 1;
        size_t J = s.paymTimes.size(); //
        double hdiff, hsqrdiff;
        for(auto const& i: s.paymTimes) 
        {
            if (j == J) break;
            z.disc[i];
            hdiff = m.H(i * mdsep) - m.H(s.expiry * mdsep);
            hsqrdiff = hdiff *(m.H(i *mdsep) + m.H(s.expiry * mdsep));
            value += z.tau[i] * z.disc[i] * exp(-hdiff * x - 0.5* hsqrdiff *m.vol[s.expiry]);
            ++j;
        }
        //VOL here for the expiry, so is fixed
        hdiff = m.H(s.paymTimes[J-1] * mdsep) - m.H(s.expiry * mdsep);
        hsqrdiff = hdiff *(m.H(s.paymTimes[J-1] *mdsep) + m.H(s.expiry * mdsep));
        value = s.strike * value;
        value += z.disc[s.paymTimes[J-1]]* exp(-hdiff * x - 0.5* hsqrdiff *m.vol[s.expiry]);
        value += -z.disc[s.paymTimes[0]] ;
        return value
    }

};
