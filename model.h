 
#pragma once 
#include <vector>
#include <map>
#include <math.h>
 
#include "../nonlinsolver/methods.h"
#include "../marketdata/mdata.h"
#include "../stats/statistics.h"
#include "../elements/dates.h"

const double mdsep = 1.0/12.0; //market data is for monthly increments

//LINEAR GAUSS MARKOV model
//or Hagan's numeraire

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

//a compiler handled version (see Gottschling Ch 5)
//"let the compiler compute"
//constexpr double lgmmeanrev(double& p, double& t )
//{
//    return (1-exp(-t*p))/p;
//}

//a parallel to the timezero data
struct lgmnumeraire  
{
    //variance is the "integral \alpha^2(t) dt" term  
    lgmmeanreversion& H;
    lgmnumeraire(lgmmeanreversion& _H): H(_H){}
    double operator()(const double& disc, const double& variance,
               const double& t, const double& x) const
    {
        double h = H(t);
        return (1/disc)*exp(h*x +0.5* h*h*variance);
    }
};

struct variancemap
{
    std::map<size_t, double, std::less<size_t>> variance;
    std::map<size_t, double, std::less<size_t>> volatility; 
    //the integral of vol^2 is variance
    //i.e. the alpha term in the documents
    variancemap(std::map<size_t, double, std::less<size_t>> _variance):
                variance(_variance){
        double val_prev;
        double time_prev;    
        for (auto k = variance.begin(), e = variance.end(); k != e; ++k){
            if(k== variance.begin()){
                volatility.insert({k->first,sqrt(k->second/k->first)});
            }
            else{
                volatility.insert({k->first,sqrt((k->second-val_prev)/(k->first- time_prev))});
            }
            val_prev = k->second;
            time_prev = k->first;
        }
	} 

    //variancemap(const variancemap& rhs):
    //    variancemap(rhs.variance){
    //}





    double operator()(size_t lookup) const
    {
        for (auto k = variance.begin(), e = variance.end(); k != e; ++k){
			if (lookup< k->first){
                return k->second;
            }
		}
        return 0.0;
    }           
};

struct lgmonefactor
{
    //a struct to initiata volatility and mean reversion     
    // let the variance vector be paired with integer excel dates
    //so (t,v) means v[0,t] where $t$ is the integer expiry date
    
    variancemap vmap;
    lgmmeanreversion H;
    //lgmonefactor (double m, std::map<size_t, double, std::less<size_t>> & variance):
    //    H(m), vmap(variance)   {} 
    lgmonefactor (double m, variancemap& _variancemap):
        H(m), vmap(_variancemap)  {} 
    //IT LOOKS like default copy-constructer for varaincemap works
};

struct swaption
{
    unsigned int expiry; //a date, as in excel 
    double strike;
    double swaptenor;     
    vector<unsigned int> swapdates; //start and all the payment dates
    double mvalue; //market price (from volatility using Black's formula)
    std::vector<unsigned int> swappaymTimes; //fixed leg payment times
    swaption(double e, double s, vector<unsigned int> dates, double val):
        expiry(e), strike(s), swaptenor(e+dates.size()), swapdates(dates),
        mvalue(val) 
        {}
};

struct lgmswapprice
{
    //using formula 5.7a from Hagan's "Evaluating and hedging ..."
    //swap is from a swaption
    //no need  for an independent swap!  
    lgmonefactor& m;  //mean reversion function H and vol data
    timezero& z; //cached discount factors
    swaption& s;
    size_t pDate;
    DayCountCalculator dayc; //is this good memory wise? maybe by reference?

    lgmswapprice(swaption& _s, timezero& _z, lgmonefactor& _m, size_t _pDate, std::string dayCMethod):
    s(_s), z(_z), m(_m),  pDate(_pDate), dayc(dayCMethod){;}

    //this is not really the swap price
    //but the one divided by exp(-H_0 y*-0.5 H_0^2 v) 
    //this division is only for rrot finding purposes, no other intention
    double operator()(double x) const //x is the X ~N(0,\alpha_t)
    {      
        //cout<<z.tau[0] << z.disc[0]<<"\n";   
        double value = 0.0;
        double hdiff, hsqrdiff;
        //cout<<"before"<<"\n";
        double expiry_yf = dayc.yFrac(pDate,s.expiry);
        for(auto const& i: s.swapdates)
        {    
            double paymdate_yf = dayc.yFrac(pDate,i);
            hdiff = m.H(paymdate_yf) - m.H(expiry_yf);
            hsqrdiff = hdiff *(m.H(paymdate_yf) + m.H(expiry_yf));
            value += z.tau(i) * z.disc(i) * exp(-hdiff * x - 0.5* hsqrdiff *m.vmap(s.expiry));
        }
        
        value *= s.strike ;
        //VOL here for the expiry, so is fixed
        size_t J = s.swapdates.size(); //
        double lastpaym_yf = dayc.yFrac(pDate,s.swapdates[J-1]);
        hdiff = m.H(lastpaym_yf) - m.H(expiry_yf);
        hsqrdiff = hdiff *(m.H(lastpaym_yf) + m.H(expiry_yf));
        value += z.disc(s.swapdates[J-1])* exp(-hdiff * x - 0.5* hsqrdiff *m.vmap(s.expiry));
        //assuming as usual xpiry =swap start date
        value += -z.disc(s.expiry) ;
         
        return value;
    }
};

//should take a model (volatility and mean rversion param)
//and some other data (strike etc.)
//and produce a swaption price
//THIS WILL BE A FUNCTOR in volatility
struct lgmswaptionprice
{
    swaption& s;
    timezero& z; 
    lgmonefactor& m; //model->meanrev
    lgmswapprice p;    

    lgmswaptionprice(swaption& _s, timezero& _z, lgmonefactor& _m, size_t pDate, std::string dayCMethod ): 
            s(_s), z(_z), m(_m), p(_s,_z,_m, pDate, dayCMethod){;}

    double operator() (double variance)  const
    {             
        //there is composition here
        //so please see if it is better to write separate functions
        //and compose them, ala functional programming
        double breakEvenRate = bisectionMethod(p,-50,50,0.0000001);
        std::cout<<"\n"<<"(breakEvenRate " << breakEvenRate <<")"<<"\n"; 
        double aux = 0;
        double h;
        double expiry_yf = p.dayc.yFrac(p.pDate,s.expiry);
        for(auto const& i: s.swapdates) 
        {
            double paymdate_yf = p.dayc.yFrac(p.pDate,i);
            h = m.H(paymdate_yf) - m.H(expiry_yf);
            aux += z.tau(i)*z.disc(i) 
                        * normalCdf((breakEvenRate + variance* h)/sqrt(variance));
            
        }
        aux *= s.strike;
        //-D_0
        aux += -z.disc(s.expiry)  * normalCdf(breakEvenRate/sqrt(variance));
        //+D_end (the last h from above loop should be the the
        //right h for below)
        size_t J = s.swapdates.size(); //
        aux += z.disc(s.swapdates[J-1])  * normalCdf((breakEvenRate + variance* h)/sqrt(variance));
        return aux -s.mvalue;
    }
};

 

