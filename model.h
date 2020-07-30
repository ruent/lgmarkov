 
#pragma once 
#include <vector>
//CPLUS_INCLUDE_PATH =~/projects 
//export CPLUS_INCLUDE_PATH
#include <math.h> 
#include "solvers/nonlinearsolvers.h"
using namespace std;




//LINEAR GAUSS MARKOV model
//or Hagan's numeraire

//since I cannot get this data natively
//I'll cache it here in this struct
//I'll assume data for each +1month, relative to
//price date
//calling functions can calculate the time multiple and retrive the
//stored numbers
struct timezero
{
    const vector<double>& tau; //period year fractions
    const vector<double>& endTimes;  
    const vector<double>& discVector; //a discount value for each endtime
};


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

struct basicswaption
{
    double swaptionExpiry; //in years = swap start time
    double swapTenor; //in years
    double strike;
    size_t numOfPeriods;
}

//should take a model (volatility and mean rversion param)
//and some other data (strike etc.)
//and produce a swaption price
//THIS WILL BE A FUNCTOR in volatility
struct lgmswaptionprice
{
    basicswaption& s;
    timezero& tZero; 
    lgmonefactor& model; //model->meanrev
    lgmswapprice swapprice;


    lgmswaptionprice() {}

    double operator() (double vol)  
    {             
        //there is composition here
        //so please see if it is better to write separate functions
        //and compose them, ala functional programming
        double breakEvenRate = bisectionMethod(swapprice,-50,50,0.0000001);
         
        double aux = 0;
        //auto e = discVector.begin(); e != discVector.end(); ++e
        for (size_t i = 0; i<s.numOfPeriods; ++i)
        {
            double h = model.H(s.paymTimes[i]) - model.H(paymTimes[0]);
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
    lgmonefactor* model;  
    timezero& tZero;
    double vol;
    double fixrate;
    double swapstarttime;
     
    double operator()(double markov) const
    {
        size_t n = paymTimes.size();
        double value = 0.0;
        for (size_t i; i<n;++n)
        {
            double hdiff = model->H(paymTimes[i]) - model->H(swapstarttime);

            value += tau[i] * discVector[i]*exp(k())
        }

        return value*fixrate;
    }

};
