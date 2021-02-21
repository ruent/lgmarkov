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
    double shift;
    lgmmeanreversion(double p, double s = 0.0): param(p), shift(s){}
    double operator() (const double& t) const //time
    {
        return shift+ (1-exp(-t*param))/param;
    }
};

struct meanReversionThreeStates //int_0^t exp(-\int_0^s param du) ds
{
    double paramA;
    double paramB;
    double paramC;
    double cutoffA; //paramA applies on the interval upto cutoffA
    double cutoffB;  
     
    double shift;
    meanReversionThreeStates(double pA, double pB, double pC, double cA, double cB, double s = 0.0): 
             paramA(pA),paramB(pB),paramC(pC), 
             cutoffA(cA),cutoffB(cB),shift(s){
                 if (paramA == 0.0 || paramB == 0.0 || paramC == 0.0){
                     std::cout<< "zero here is not a good idea!"<<"\n";   
                 }   
             }
    double operator() (const double& t) const //time
    {
        //e^{-k_1 psi_1 -k_2 psi_2 + k_3 (psi_1+psi_2)}
        double x;
        if(t <= cutoffA){
            x = (1-exp(-t*paramA))/paramA;
        }
        else if(t <= cutoffB){
            x = (1-exp(-cutoffA*paramA))/paramA;
            x += exp(-paramA * cutoffA+ paramB * cutoffA) 
                           *(exp(-paramB*cutoffA)-exp(-t*paramB)) / paramB;
        }
        else{
            x = (1-exp(-cutoffA*paramA))/paramA;
            x += exp(-paramA * cutoffA+ paramB * cutoffA) 
                           *(exp(-paramB*cutoffA)-exp(-cutoffB*paramB)) / paramB;
            x += exp(-paramA * cutoffA - paramB * (cutoffB-cutoffA) + paramC*cutoffB) 
                        *(exp(-paramC*cutoffB)-exp(-t*paramC))/paramC;        
        }      
        return shift+ x;
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
    int initdate; //should be an excel date, it is used explicitly in vol
    //presence of initdate is implicit in variance as it will be coming from
    //a calculation involving this date; still it is a good idea to put it
    //somewhere in there
    //FIX AN ANNUAL VOL GAUGE! ALL THROUGHOUT!
    //THIS REQUIRES EXCEL DATE differences to be turned into year fractions
    std::map<int, double, std::less<int>> variance; 
    //variance is time frac times annual vol times annual vol
    std::map<int, double, std::less<int>> volatility; //annualized vol
    //the integral of vol^2 is variance
    //i.e. the integral of alpha^2 term in the documents
    variancemap(std::map<int, double, std::less<int>> _variance, int _initdate):
                variance(_variance), initdate(_initdate){
                double val_prev;
                double time_prev;    
                for (auto k = variance.begin(), e = variance.end(); k != e; ++k){
                    if(k== variance.begin()){
                        volatility.insert({k->first,sqrt((k->second)*365/(k->first- initdate))});
                    }
                    else{
                        volatility.insert({k->first,sqrt((k->second-val_prev)*365/(k->first- time_prev))});
                    }
                    val_prev = k->second;
                    time_prev = k->first;
                }
	}     

    //so, think of this as a pw-flat curve, that is right-continuous
    //the first time map contains an excel date higher than the lookup value
    //corresponding number is your return
    double operator()(size_t lookup) const {
        bool needForExtrapolation = true;
        for (auto k = volatility.begin(), e = volatility.end(); k != e; ++k)
		{
			if (lookup< k->first){
                return k->second;
                needForExtrapolation = false;
            }
		}
        if (needForExtrapolation){
            //cout<< "Extrpolating from the last value: " << variance.end()->first ;
            //cout<<  " " <<variance.end()->second <<"\n";
            return volatility.end()->second;
        }
    }
    // THIS PART IS NEEDED to be an integral of volatility IF YOU EVER NEED TOTAL VARIANCE 
    // IN BETWEEN swap dates; for now, due to urgency it is just plain lookup
    double total_variance(size_t lookup) const {
        bool needForExtrapolation = true;
        double total_variance = 0.0;
        for (auto k = variance.begin(), e = variance.end(); k != e; ++k){
			if (lookup< k->first){
                return k->second;
                needForExtrapolation = false;
            }
		}
        if (needForExtrapolation){
            return variance.end()->second;
        }
    }     
};

struct lgmonefactor
{
    //a struct to initiata volatility and mean reversion     
    // let the variance vector be paired with integer excel dates
    //so (t,v) means v[0,t] where $t$ is the integer expiry date    
    variancemap vmap;    
    meanReversionThreeStates H; //initiated with meanreversion parameter m, and a shift s
    lgmonefactor(double A,double B,double C, size_t c, size_t d, 
              double s, std::map<int, double, std::less<int>> variance, int initdate):
        H(A,B,C, c,d), vmap(variance, initdate)   {}     
};

struct swaption
{
    unsigned int expiry; 
    double strike;
    double swaptenor;     
    vector<unsigned int> swapdates; //start and all the payment dates
    //payment dates should be fixed leg dates
    //the reason is float leg, assuming a single-currency framework ...
    double mvalue; //market price (from volatility using Black's formula)
    std::vector<unsigned int> swappaymTimes; //fixed leg payment times
    bool isReceiver;
    swaption(double e, double s, vector<unsigned int> dates, double val, bool isReceiver_):
        expiry(e), strike(s), swaptenor(e+dates.size()), swapdates(dates),
        mvalue(val), isReceiver(isReceiver_)
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
    DayCountCalculator dayc; //is this good memory wise? maybe by reference? maybe a miniscule differece only!

    lgmswapprice(swaption& _s, timezero& _z, lgmonefactor& _m, size_t _pDate, std::string dayCMethod):
    s(_s), z(_z), m(_m),  pDate(_pDate), dayc(dayCMethod){;}

    //this is not really the swap price, but swap price multiplied discount at expiry
    double operator()(double x) const {//x is the X ~N(0,\alpha_t)
        double value = 0.0;
        double hdiff;
        //cout<<"before"<<"\n";
        double expiry_yf = dayc.yFrac(pDate,s.expiry);
        double some_term = 0.0;
        for(auto const& i: s.swapdates){               
            double paymdate_yf = dayc.yFrac(pDate,i);
            hdiff = m.H(paymdate_yf) - m.H(expiry_yf);
            value += z.tau(i) * z.disc(i) * exp(-hdiff * x - 0.5* hdiff*hdiff *m.vmap.total_variance(s.expiry));
            some_term += z.tau(i) * z.disc(i) *s.strike;
        }        
        value *= s.strike;
        //VOL here for the expiry, so is fixed
        size_t J = s.swapdates.size(); //
        double lastpaym_yf = dayc.yFrac(pDate,s.swapdates[J-1]);
        hdiff = m.H(lastpaym_yf) - m.H(expiry_yf);
        double p_last = z.disc(s.swapdates[J-1])* exp(-hdiff * x - 0.5* hdiff *hdiff *m.vmap.total_variance(s.expiry));
        value += p_last;
        //assuming as usual xpiry =swap start date
        value += -z.disc(s.expiry) ;  
        some_term += -z.disc(s.expiry) +p_last;
        std::cout << "swap price ... approximate some_term: " << some_term<<"\n";
        if (s.isReceiver){
            return value;
        }
        else return -value;
    }

    //operator() above is swapprice times discount at expiry of the swaption (swap start date).
    //rebased swap price is to be used in the pde solver
    //rebasing like this eliminates the need for the numeraire but please make
    //sure you rebase it explicitly in future versions...
    double rebased_swap_price(double x) const{
        double value{};
        double h;
        double expiry_yf = dayc.yFrac(pDate,s.expiry);
        double value_to_compare{};
        for(auto const& i: s.swapdates){               
            double paymdate_yf = dayc.yFrac(pDate,i);
            h = m.H(paymdate_yf) ;
            value += z.tau(i) * z.disc(i) * exp(-h * x - 0.5* h*h *m.vmap.total_variance(s.expiry));
            value_to_compare += z.tau(i) * z.disc(i); 
        }        
        value *= s.strike;
        value_to_compare *= s.strike;
        //VOL here for the expiry, so is fixed
        size_t J = s.swapdates.size(); //
        double lastpaym_yf = dayc.yFrac(pDate,s.swapdates[J-1]);
        h = m.H(lastpaym_yf) ;
        double p_last = z.disc(s.swapdates[J-1])* exp(-h * x - 0.5* h * h *m.vmap.total_variance(s.expiry));
        value += p_last;
        value_to_compare += z.disc(s.swapdates[J-1]);
        h = m.H(expiry_yf);
        value += -z.disc(s.expiry) * exp(-h * x - 0.5* h * h *m.vmap.total_variance(s.expiry));
        value_to_compare += -z.disc(s.expiry) ;
        std::cout << " rebased value: " << value << "plain value: " <<value_to_compare <<"\n";
        if (s.isReceiver){
            return value;
        }
        else return -value;
        
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
    size_t pDate;
    std::string dayCMethod;
    //found variance and the underlying swap breakeven state value can be saved as 
    //memeber variables, but this requires to remove the const from operator()
    //not going to do that for now
    double fixed_point_breakeven;
    double fixed_point_variance;

    lgmswaptionprice(swaption& _s, timezero& _z, lgmonefactor& _m, size_t _pDate, std::string _dayCMethod): 
            s(_s), z(_z), m(_m), pDate(_pDate), dayCMethod(_dayCMethod){;}

    double operator() (double variance){    
        //the 7 below means variance will prevail for expiry plus a week
        std::map<int, double, std::less<int>> _variance{{s.expiry+7, variance}};
        lgmonefactor _m(m.H.paramA,m.H.paramB,m.H.paramC,m.H.cutoffA,m.H.cutoffB, 
                            m.H.shift, _variance, pDate); //model->meanrev  
        lgmswapprice p(s,z,_m, pDate, dayCMethod);
        //there is composition here
        //so please see if it is better to write separate functions
        //and compose them, ala functional programming
        std::cout<<" " <<"\n"; 
        std::cout<<"Underlying variance: " << variance << "do brute force bisection"<<"\n";         
        double k = 0.0;
        for(size_t i=0; i<999; ++i){
            k = i*0.1;
            if (p(-k) *p(k)<0){
                std::cout<<"found something: " << k<<"\n";
                break;
            } 
        } 
        double breakEvenRate = bisectionMethodNonConst(p,-k,k,0.0000001);   
        fixed_point_breakeven = breakEvenRate; 
        fixed_point_variance = variance;       
        std::cout<<"breakEvenRate: " << breakEvenRate <<"\n";
        double aux = 0;
        double h;
        double expiry_yf = p.dayc.yFrac(p.pDate,s.expiry);
        double some_term = 0.0;
        for(auto const& i: s.swapdates){
            double paymdate_yf = p.dayc.yFrac(p.pDate,i);
            h = m.H(paymdate_yf) - m.H(expiry_yf);
            aux += z.tau(i)*z.disc(i) 
                        * normalCdf((breakEvenRate + variance* h)/sqrt(variance));   
            some_term += z.tau(i)*z.disc(i) *s.strike;
        }        
        aux *= s.strike;
        std::cout<<"intermediate steps: sums " << aux <<"\n"; 
        //-D_0
        aux += -z.disc(s.expiry)  * normalCdf(breakEvenRate/sqrt(variance));
        std::cout<<"intermediate steps: D_0 term: " << -z.disc(s.expiry)  * normalCdf(breakEvenRate/sqrt(variance)) <<"\n"; 
        //+D_end (the last h from above loop should be the the
        //right h for below)
        size_t J = s.swapdates.size(); //
        aux += z.disc(s.swapdates[J-1])  * normalCdf((breakEvenRate + variance* h)/sqrt(variance));
        some_term += -z.disc(s.expiry)+ z.disc(s.swapdates[J-1])  ;
        std::cout<<"swaption price: h: " << h << " variance: "<< variance  <<"\n"; 
        std::cout<<"swaption price: D_n term: " << z.disc(s.swapdates[J-1])  * normalCdf((breakEvenRate + variance* h)/sqrt(variance)) <<"\n"; 
        std::cout<<"swaption price aux: "<< aux << "target value: " << s.mvalue <<"\n"; 
        std::cout<< "swaption price D_n -D_0 +sum : ... " << some_term<<"\n"; 
        if (!s.isReceiver){             
            for(auto const& i: s.swapdates){
                aux += - -z.tau(i)*z.disc(i) *s.strike;            
            } 
            aux += z.disc(s.expiry) -z.disc(s.swapdates[J-1]);
        }
        return aux -s.mvalue;
    }
};