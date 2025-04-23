// Hg Pup Production Model V1b (Does not include Stage 4 + Skew-Normal Birth Curve)
// Last modified 2024-02-14 by EKJ

#include <TMB.hpp>                                // Links in the TMB libraries
#include <cmath>
#include <iostream>


// integrate skew normal
// https://kaskr.github.io/adcomp/namespaceromberg.html
template<class Type>
struct int_dsn {
  // Parameters in integrand (things you're not integrating over)
  Type xi, omega, alphabday;
  // Evaluate integrand
  Type operator()(Type x){
    return dsn((x - xi)/omega, alphabday, false)/omega;
  }
};

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  using namespace Eigen;
  using namespace atomic;
  using namespace Rmath;

  DATA_VECTOR(nwhite);                       // number of white pups counted
  DATA_VECTOR(nmoult);                       // number of moulted pups counted
  DATA_IVECTOR(obsdays);                      // days on which counts happened
  
  DATA_VECTOR(pmoult);                        // probability of ending moult in days post-start
  DATA_VECTOR(pleave);                        // probability of leaving at each age
  DATA_VECTOR(bday);                          // this is just a vector 1:ndays

  DATA_SCALAR(pcwhite);                       // probability of being counted as white
  DATA_SCALAR(pcmoult);                       // probability of being counted as moulted
  DATA_SCALAR(powhite);                       // probability of being counted if white
  DATA_SCALAR(pomoult);                       // probability of being counted if moulted
  DATA_INTEGER(ndays);                        // days of the season
  //DATA_INTEGER(age);                          // maximum age of pups
  
  PARAMETER(N);                               // total number born
  PARAMETER(mubday);                          // mean birthday
  PARAMETER(sdbday);                          // sd birthday
  PARAMETER(alphabday);                       // skew parameter for birth curve
  
  vector<Type> born(ndays);                   // unobserved true # born on each day
  born.setZero();
  vector<Type> totPups(ndays);                // unobserved true # on beach
  totPups.setZero();
  vector<Type> totWhite(ndays);               // unobserved true # whitecoats
  totWhite.setZero();
  vector<Type> totMoult(ndays);               // unobserved true # moulted
  totMoult.setZero();
  vector<Type> totBorn(ndays);                // unobserved true # cumulatively born on each day
  totBorn.setZero();
  vector<Type> pborn(ndays);
  pborn.setZero();
  vector<Type> pborn1(ndays);
  pborn.setZero();
  vector<Type> Ew(ndays);               // expected # whitecoats obs
  Ew.setZero();
  vector<Type> Em(ndays);               // expected # moulted pups obs
  Em.setZero();
  
  vector<Type> Emm(ndays);               // expected # moulted pups obs
  Emm.setZero();
  vector<Type> Ewm(ndays);               // expected # moulted pups obs
  Ewm.setZero();
  
  vector<Type> wpupvec(obsdays.size());               // expected # obs on obsdays
  wpupvec.setZero();
  vector<Type> mpupvec(obsdays.size());               // expected # obs on moultdays
  mpupvec.setZero();
  
  vector<Type> sdwhitevec(ndays);               // vector to save sdwhite
  sdwhitevec.setZero();
  vector<Type> sdmoultvec(ndays);               // vector to save sdmoult
  sdmoultvec.setZero();

  Type cumborn = 0.0;                         // cumulative # born (single val)
  Type wpups = 0;                                 // pups available to be counted (single val)
  Type mpups = 0;
  Type wobs = 0;                                   // number of pups counted (single val)
  Type mobs = 0;
  Type nll = 0.0;                             // Negative log-likelihood

  Type delta = 0.0;
  Type omega = 0.0;
  Type xi = 0.0;
  
  Type mumday = 23.0;
  Type sdmday = 5.0;
  Type mulday = 31.5;
  Type sdlday = 7.0;
  
  // convert to params for skew normal
  
  delta = alphabday/sqrt(1+pow(alphabday,2));
  omega = sqrt(pow(sdbday,2)/(1-(2*pow(delta,2))/3.14159));
  xi = mubday - omega*delta*sqrt(2/3.14159);
  
  // integrate birth curve for each day of the season
  int_dsn<Type> f={xi, omega, alphabday};
  for(int i = 0; i < ndays; i++){
    Type r = i;
    pborn(i) = romberg::integrate(f, r, r+1.0);
  }
  
  // standardize birth curve so that all pups are born
  
  for(int i = 1; i < ndays; i++){
    
    pborn1(i) = pborn(i)/sum(pborn);
    
  }
  
  // Birth process model
  for(int i = 0; i < ndays; i++){
    
    born(i) = pborn1(i) * N;                // deterministic birth process
    cumborn += born(i);                     // add to pups already born
    totBorn(i) = cumborn;                   // save cumulative # born by this day

    
  } // end process model
  
  // moult and leaving process model
  
  for(int d = 0; d < ndays; d++){
    
    for(int i = 0; i < d+1; i++){
      
      Type age = d-i;
      
      totWhite(d) += born(i) * (1 - pnorm(age, mumday, sdmday)) * (1 - pnorm(age, mulday, sdlday));
      
      totMoult(d) += born(i) * pnorm(age, mumday, sdmday) * (1 - pnorm(age, mulday, sdlday));
      
    }
    
    totPups(d) = totWhite(d) + totMoult(d);
    
    Ew(d) = totWhite(d)*powhite*pcwhite + totMoult(d)*pomoult*(1-pcmoult);
    Em(d) = totMoult(d)*pomoult*pcmoult + totWhite(d)*powhite*(1-pcwhite);

  }
  

  // Observation model
  
  int x = obsdays.size();

  for(int i=0; i<x; i++){

  // Calculate # pups observed/classified as white
  wpups = (totWhite(obsdays(i))*powhite*pcwhite) +              // s1-3 pups can't be misIDed
          (totMoult(obsdays(i))*pomoult*(1-pcmoult));         // Moulted may be IDed as S4
  
  // Calculate binomial variance
  // NOTE the sqrt(pow(n)) args are to make sure variance isn't negative
  Type vwhite = sqrt(pow(totWhite(obsdays(i)), 2))*powhite*(1-powhite) + // var w cc
                sqrt(pow(totMoult(obsdays(i)), 2))*pomoult*(1-pcmoult)*(1-(pomoult*(1-pcmoult))); // var m ic
    
  Type sdwhite = sqrt(vwhite);
  //Type sdwhite = 0.1;
  Type wobs = nwhite(i);                          // true # obs = data
  
// calculate likelihood, add to nll
  //nll -= log(pnorm(wobs+0.5, wpups, sdwhite) - pnorm(wobs-0.5, wpups, sdwhite));
  nll -= dnorm(wobs, wpups, sdwhite, true);
     
  // Calculate # pups observed/classified as moulted
  mpups = (totMoult(obsdays(i))*pomoult*pcmoult) + (totWhite(obsdays(i))*powhite*(1-pcwhite));
    
  // Calculate binomial variance
  // NOTE the sqrt(pow(n)) args are to make sure variance isn't negative
  Type vmoult = sqrt(pow(totMoult(obsdays(i)),2))*pomoult*pcmoult*(1-(pomoult*pcmoult)) + // var m cc
                sqrt(pow(totWhite(obsdays(i)), 2))*powhite*(1-pcwhite)*(1-(powhite*(1-pcwhite))); // var S4 ic
     
  Type sdmoult = sqrt(vmoult);
  //Type sdmoult = 0.1;
  Type mobs = nmoult(i);                          // true # obs = data
   
  //nll -= log(pnorm(mobs+0.5, mpups, sdmoult) - pnorm(mobs-0.5, mpups, sdmoult));
  nll -= dnorm(mobs, mpups, sdmoult, true);

  wpupvec(i) = wpups;
  mpupvec(i) = mpups;
  sdwhitevec(i) = sdwhite;
  sdmoultvec(i) = sdmoult;
     
 }
  

  REPORT(pborn1); 
  REPORT(born);
  REPORT(totBorn);
  REPORT(totPups);
  REPORT(totWhite);
  REPORT(totMoult);
  REPORT(wpupvec);
  REPORT(mpupvec);
  REPORT(sdwhitevec);
  REPORT(sdmoultvec);
  REPORT(Ew);
  REPORT(Em);
  REPORT(Emm);
  REPORT(Ewm);
  
  ADREPORT(mubday);
  ADREPORT(sdbday);
  ADREPORT(alphabday);
  ADREPORT(N);
  
  return nll;                                 // negative log-likelihood is returned

  } 
