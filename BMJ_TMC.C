#include "TMath.h"
#include "TArrayD.h"
#include "BMJ_TMC.h"
#include <iostream>
#include <fstream>

Double_t tmin(Double_t q2, Double_t x){
 x_2=TMath::Power(x,2);
 eps=2*x*M/TMath::Sqrt(q2);
 eps_2=TMath::Power(eps,2);
 t_min=-q2*(2*(1-x)*(1-TMath::Sqrt(1+eps_2))+eps_2)/(4*x*(1-x)+eps_2);
 return t_min;
}

Double_t F1p(Double_t q2) { return ((-1.)*q2/4./TMath::Power(0.938271998,2.)*GMp(q2)+GEp(q2))/(1.-q2/4./TMath::Power(0.938271998,2.)); }
Double_t F2p(Double_t q2) { return ( GMp(q2)-GEp(q2))/(1.-q2/4./TMath::Power(0.938271998,2.)); }
/*Double_t F1p(Double_t q2) { return (4.*TMath::Power(0.938272,2)-2.79285*q2)/(TMath::Power(1-1.4084507042253522*q2,2)*(4.*TMath::Power(0.938272,2)-q2)); }
  Double_t F2p(Double_t q2) { return (7.1714*TMath::Power(0.938272,2))/(TMath::Power(1-1.4084507042253522*q2,2)*(4.*TMath::Power(0.938272,2)-q2)); }*/

Double_t GMp(Double_t q2) { return KellyM(q2); }
Double_t GEp(Double_t q2) { return KellyE(q2); }


Double_t KellyE(Double_t q2)
  /* JJ Kelly PRC 70, 068202 (2004)
   */
{
  Int_t ia=1, ib=3;
  Double_t a[1]={-0.24};
  Double_t b[3]={10.98, 12.82, 21.97};
  Double_t Mp = 0.938;
  Double_t tau = -q2/(4.*pow(Mp,2));
  Double_t GEKn = 1.0;
  Double_t GEKd = 1.0;
  for (int ja=0; ja<ia; ja++){
    GEKn+=a[ja]*pow(tau,ja+1);
  }
  for (int jb=0; jb<ib; jb++){
    GEKd+=b[jb]*pow(tau,jb+1);
  }
  return GEKn/GEKd;
}

Double_t KellyM(Double_t q2)
  /* JJ Kelly PRC 70, 068202 (2004)
     Magnetic Form Factor fit
     Returned value is ratio to dipole*mu_p
  */
{
  Int_t ia=1, ib=3;
  Double_t a[1]={0.12};
  Double_t b[3]={10.97, 18.86, 6.55};Double_t Mp = 0.938;
  Double_t tau = -q2/(4.*pow(Mp,2));
  Double_t GMKn = 1.0;
  Double_t GMKd = 1.0;
  for (int ja=0; ja<ia; ja++)
    {
      GMKn+=a[ja]*pow(tau,ja+1);
    }
  for (int jb=0; jb<ib; jb++)
    {
      GMKd+=b[jb]*pow(tau,jb+1);
    }
  return 2.79285*GMKn/GMKd;
}


Double_t BH(Double_t q2, Double_t x, Double_t t, Double_t kz, Double_t phig){
  //From BMK, BH unpolarized target case ( 
  yy=q2/2./M/x/kz;
  //phig=TMath::Pi()-phig; // conversion phi_trento a phi_BMK
  TArrayD coeff(3);
  F1_2=TMath::Power(F1p(t),2);
  F2_2=TMath::Power(F2p(t),2);
  F_sum_2=TMath::Power(F2p(t)+F1p(t),2);
  x_2=TMath::Power(x,2);
  eps=2*x*M/TMath::Sqrt(q2);
  eps_2=TMath::Power(eps,2);
  J=(1-yy-yy*eps_2/2.)*(1+t/q2)-(1-x)*(2-yy)*t/q2;
  
  t_min=-q2*(2*(1-x)*(1-TMath::Sqrt(1+eps_2))+eps_2)/(4*x*(1-x)+eps_2);
  
  K_2=-t*(1-x)*(1-yy-TMath::Power(yy*eps/2.,2))*(1-t_min/t)*(TMath::Sqrt(1+eps_2)+(4*x*(1-x)+eps_2)/(4*(1-x))*(t-t_min)/q2)/q2;
  
  K=TMath::Sqrt(K_2);

  P1=-(J+2*K*TMath::Cos((TMath::Pi()-phig)))/(yy*(1+eps_2)); 
  P2=1+t/q2+(J+2*K*TMath::Cos((TMath::Pi()-phig)))/(yy*(1+eps_2));
  prefactor_BH=hbar_c*TMath::Power(alpha,3)/(TMath::Pi()*8*q2*TMath::Sqrt(1+eps_2))/x/yy/TMath::Power(1+eps_2,2)/t/P1/P2;
  prefactor_BH=prefactor_BH*Jacobian;
  coeff[0]=
    8*K_2*((2+3*eps_2)*q2*(F1_2-t*F2_2/(4*M_2))/t+2*x_2*F_sum_2)
    +TMath::Power(2-yy,2)*((eps_2+2)*(4*x_2*M_2/t*TMath::Power(1+t/q2,2)+4*(1-x)*(1+x*t/q2))*(F1_2-t/4./M_2*F2_2)+4*x_2*(x+(1-x+eps_2/2.)*TMath::Power(1-t/q2,2)-x*(1-2*x)*TMath::Power(t/q2,2))*F_sum_2)
    +8*(1+eps_2)*(1-yy-TMath::Power(eps*yy/2.,2))*(2*eps_2*(1-t/4./M_2)*(F1_2-t/4./M_2*F2_2)-x_2*TMath::Power(1-t/q2,2)*F_sum_2);

  coeff[1]=8*K*(2-yy)*((4*x_2*M_2/t-2*x-eps_2)*(F1_2-t/4./M_2*F2_2)+2*x_2*(1-(1-2*x)*t/q2)*F_sum_2);

  coeff[2]=8*x_2*K_2*(4*M_2/t*(F1_2-t/4./M_2*F2_2)+2*F_sum_2);
  return prefactor_BH*(coeff[0]+TMath::Cos((TMath::Pi()-phig))*coeff[1]+TMath::Cos(2*(TMath::Pi()-phig))*coeff[2]);
}

TArrayD BH_DVCS_interference(Double_t q2, Double_t x, Double_t t, Double_t kz, Double_t phig, Int_t helstate, Int_t part){
  yy=q2/2./M/x/kz;
  //phig=TMath::Pi()-phig; // conversion phi_trento to phi_BMK
  TArrayD coeff_CFF(2);
  x_2=TMath::Power(x,2);
  eps=2*x*M/TMath::Sqrt(q2);
  eps_2=TMath::Power(eps,2);
  t_min=-q2*(2*(1-x)*(1-TMath::Sqrt(1+eps_2))+eps_2)/(4*x*(1-x)+eps_2); 
  K_2=-t*(1-x)*(1-yy-TMath::Power(yy*eps/2.,2))*(1-t_min/t)*(TMath::Sqrt(1+eps_2)+(4*x*(1-x)+eps_2)/(4*(1-x))*(t-t_min)/q2)/q2;
  K=TMath::Sqrt(K_2); 
  J=(1-yy-yy*eps_2/2.)*(1+t/q2)-(1-x)*(2-yy)*t/q2;
  P1=-(J+2*K*TMath::Cos((TMath::Pi()-phig)))/(yy*(1+eps_2)); 
  P2=1+t/q2+(J+2*K*TMath::Cos((TMath::Pi()-phig)))/(yy*(1+eps_2));

  prefactor_I=hbar_c*TMath::Power(alpha,3)*x*yy/(TMath::Pi()*8*q2*TMath::Sqrt(1+eps_2))/x/TMath::Power(yy,3)/t/P1/P2;
  prefactor_I=Jacobian*prefactor_I;
  Kp_2=q2*K_2/(1-yy-TMath::Power(eps*yy/2.,2));
  Kp=TMath::Sqrt(Kp_2);
  //Double_t p_eff=TMath::Sqrt(2)/(2-x)*Kp/TMath::Sqrt(q2); 
  Double_t p_eff=TMath::Sqrt(2*(1-yy-TMath::Power(yy*eps/2.,2)))/TMath::Power(1+eps_2,2.);
  coeff_CFF[0]=0; coeff_CFF[1]=0;
  ///////////////// Part 1 Helicity conserved ++
  //Following coefficient are just for H
  if (helstate==0){
    if (part==0){
      coeff_CFF[0]=-4*(2-yy)*(1+TMath::Sqrt(1+eps_2))/TMath::Power(1+eps_2,2)*(Kp_2/q2*TMath::Power(2-yy,2)/TMath::Sqrt(1+eps_2)+t/q2*(1-yy-TMath::Power(eps*yy/2.,2))*(2-x)*(1+(2*x*(2-x+(TMath::Sqrt(1+eps_2)-1)/2.+eps_2/2./x)*t/q2+eps_2)/(2-x)/(1+TMath::Sqrt(1+eps_2)))); //n=0
     
      coeff_CFF[0]+=(-16*K*(1-yy-TMath::Power(eps*yy/2.,2))/TMath::Power(1+eps_2,2.5)*(x*t/q2*(1+(1-x)*(TMath::Sqrt(eps_2+1)-1)/2./x+eps_2/4./x)-3*eps_2/4.)-4*K*(2-2*yy+yy*yy+eps_2*yy*yy/2.)*(1+TMath::Sqrt(1+eps_2)-eps_2)/TMath::Power(1+eps_2,2.5)*(1-(1-3*x)*t/q2+x*t/q2*(1-TMath::Sqrt(1+eps_2)+3*eps_2)/(1+TMath::Sqrt(1+eps_2)-eps_2)))*TMath::Cos((TMath::Pi()-phig)); 

      coeff_CFF[0]+=(8*(2-yy)*(1-yy-TMath::Power(eps*yy/2.,2))/TMath::Power(1+eps_2,2)*(Kp_2/q2*2*eps_2/(1+eps_2+TMath::Sqrt(1+eps_2))+x*t*(t-t_min)/q2/q2*(1-x-(TMath::Sqrt(1+eps_2)-1)/2.+eps_2/2./x)))*TMath::Cos(2*(TMath::Pi()-phig)); 

      coeff_CFF[0]+=(-8*K*(1-yy-TMath::Power(eps*yy/2.,2))*(TMath::Sqrt(1+eps_2)-1)/TMath::Power(1+eps_2,2.5)*(t/q2*(1-x)+(TMath::Sqrt(1+eps_2)-1)/2.*(1+t/q2)))*TMath::Cos(3*(TMath::Pi()-phig)); 

      coeff_CFF[1]=(8*K*(2-yy)*yy/(1+eps_2)*(1+(1-x+(TMath::Sqrt(1+eps_2)-1)/2.)/(1+eps_2)*(t-t_min)/q2))*TMath::Sin((TMath::Pi()-phig));
      coeff_CFF[1]+=-(4*yy*(1-yy-TMath::Power(eps*yy/2.,2))/TMath::Power(1+eps_2,1.5)*(1+TMath::Sqrt(1+eps_2)-2*x)*(t-t_min)/q2*((eps_2-x*(TMath::Sqrt(1+eps_2)-1))/(1+TMath::Sqrt(eps_2+1)-2*x)-(2*x+eps_2)/(2*TMath::Sqrt(1+eps_2))*(t-t_min)/q2))*TMath::Sin(2*(TMath::Pi()-phig));   
   
    }

 
  
   //Here I have implemented two CFF according to hotfix which give delta C when Q goes to infinity
  //C^V
    if (part==1){
      coeff_CFF[0]=8*(2-yy)*x*t/TMath::Power(1+eps_2,2)/q2*(TMath::Power(2-yy,2)*Kp_2/q2/TMath::Sqrt(1+eps_2)+(1-yy-yy*yy*eps_2/4.)*(1+TMath::Sqrt(1+eps_2))*(1+t/q2)*(1+(TMath::Sqrt(1+eps_2)-1+2*x)*t/q2/(1+TMath::Sqrt(1+eps_2)))/2.); //n=0
      coeff_CFF[0]+=16*K*x*t/q2/TMath::Sqrt(TMath::Power(1+eps_2,5))*(TMath::Power(2-yy,2)*(1-(1-2*x)*t/q2)+(1-yy-yy*yy*eps_2/4.)*(1-2*x+TMath::Sqrt(1+eps_2))*(t-t_min)/q2/2.)*TMath::Cos((TMath::Pi()-phig)); //n=1
      coeff_CFF[0]+=8*(2-yy)*(1-yy-yy*yy*eps_2/4.)*x*t/q2/TMath::Power(1+eps_2,2)*(4*Kp_2/q2/TMath::Sqrt(1+eps_2)+(1+TMath::Sqrt(1+eps_2)-2*x)*(1+t/q2)*(t-t_min)/q2/2.)*TMath::Cos(2*(TMath::Pi()-phig)); //n=2
      coeff_CFF[0]+=-8*K*(1-yy-yy*yy*eps_2/4.)*x*t/q2/TMath::Sqrt(TMath::Power(1+eps_2,5))*(TMath::Sqrt(1+eps_2)-1+(1+TMath::Sqrt(1+eps_2)-2*x)*t/q2)*TMath::Cos(3*(TMath::Pi()-phig)); //n=3
      coeff_CFF[1]=-TMath::Sin((TMath::Pi()-phig))*8*K*(2-yy)*yy*x*t/TMath::Power(1+eps_2,2)/q2*(TMath::Sqrt(1+eps_2)-1+t*(1-2*x+TMath::Sqrt(1+eps_2))/q2);
      coeff_CFF[1]+=-TMath::Sin(2*(TMath::Pi()-phig))*4*(1-yy-yy*yy*eps_2/4.)*yy*x*t/q2/TMath::Power(1+eps_2,2)*(1-(1-2*x)*t/q2)*(TMath::Sqrt(1+eps_2)-1+t*(1-2*x+TMath::Sqrt(1+eps_2))/q2);
    }
    if (part==2){
  //C^A
      coeff_CFF[0]=8*(2-yy)*t/q2/TMath::Power(1+eps_2,2)*(TMath::Power(2-yy,2)*Kp_2*(1+TMath::Sqrt(1+eps_2)-2*x)/q2/TMath::Sqrt(1+eps_2)/2.+(1-yy-yy*yy*eps_2/4.)*((1+TMath::Sqrt(1+eps_2))*(1+TMath::Sqrt(1+eps_2)-x+t*(TMath::Sqrt(1+eps_2)-1+x*(3+TMath::Sqrt(1+eps_2)-2*x)/(1+TMath::Sqrt(1+eps_2)))/q2)/2.-2*Kp_2/q2)); //n=0
      coeff_CFF[0]+=-16*K*t/q2/TMath::Power(1+eps_2,2)*((1-yy-yy*yy*eps_2/4.)*(1-(1-2*x)*t/q2+(4*x*(1-x)+eps_2)*(t-t_min)/(4.*TMath::Sqrt(1+eps_2))/q2)-TMath::Power(2-yy,2)*(1-x/2.+(1-t/q2)*(1+TMath::Sqrt(1+eps_2)-2*x)/4.+(t-t_min)*(4*x*(1-x)+eps_2)/q2/2./TMath::Sqrt(1+eps_2)))*TMath::Cos((TMath::Pi()-phig)); //n=1
      coeff_CFF[0]+=4*(2-yy)*(1-yy-yy*yy*eps_2/4.)*t/q2/TMath::Power(1+eps_2,2)*(4*Kp_2*(1-2*x)/q2/TMath::Sqrt(1+eps_2)-(3-TMath::Sqrt(1+eps_2)-2*x+eps_2/x)*x*(t-t_min)/q2)*TMath::Cos(2*(TMath::Pi()-phig)); //n=2
      coeff_CFF[0]+=16*K*(1-yy-yy*yy*eps_2/4.)*t*(t-t_min)/q2/q2/TMath::Power(TMath::Sqrt(1+eps_2),5)*(x*(1-x)+eps_2/4.)*TMath::Cos(3*(TMath::Pi()-phig)); //n=3 
      coeff_CFF[1]=TMath::Sin((TMath::Pi()-phig))*8*K*(2-yy)*yy*t/(1+eps_2)/q2*(1-(1-2*x)*(t-t_min)*(1-2*x+TMath::Sqrt(1+eps_2))/q2/2./(1+eps_2));
      coeff_CFF[1]+=-TMath::Sin(2*(TMath::Pi()-phig))*8*(1-yy-yy*yy*eps_2/4.)*yy*t*(t-t_min)/q2/q2/TMath::Power(1+eps_2,2)*(1-x/2.+3*eps_2/4.)*(1+TMath::Sqrt(1+eps_2)-2*x)*(1+t*(4*x*(1-x)+eps_2)/q2/(4-2*x+3*eps_2));
    }

  }
  if (helstate==1){
    //Following ones are for 0+
    //coeff_CFF (retrait de peff car on regarde directement F0+ et non Feff
    if (part==0){
      coeff_CFF[0]=p_eff*12*K*(2-yy)*(eps_2+(2-6*x-eps_2)*t/3./q2)/TMath::Power(1+eps_2,0.5); //n=0
      coeff_CFF[0]+=p_eff*8*(TMath::Power(2-yy,2)*(t-t_min)/q2*(1-x+((1-x)*x+eps_2/4.)*(t-t_min)/TMath::Sqrt(1+eps_2)/q2)+(1-yy-TMath::Power(yy*eps/2.,2))/TMath::Sqrt(1+eps_2)*(1-(1-2*x)*t/q2)*(eps_2-2*(1+eps_2/2./x)*x*t/q2))*TMath::Cos((TMath::Pi()-phig)); //n=1 cos
      coeff_CFF[0]+=-p_eff*8*K*(2-yy)*(1+eps_2/2.)/TMath::Power(1+eps_2,0.5)*(1+(1+eps_2/2./x)*x*t/(1+eps_2/2.)/q2)*TMath::Cos(2*(TMath::Pi()-phig)); //n=2 cos 
      coeff_CFF[1]=p_eff*8*(2-yy)*yy*Kp_2/q2*TMath::Sin((TMath::Pi()-phig)); //n=1 sin      
      coeff_CFF[1]+=p_eff*8*K*yy*(1+eps_2/2.)*(1+(1+eps_2/2./x)*x*t/(1+eps_2/2.)/q2)*TMath::Sin(2*(TMath::Pi()-phig)); //n=2 sin
    }

    if (part==1){
    //coeff_CFF V
      coeff_CFF[0]=p_eff*24*K*(2-yy)*x*t/q2/TMath::Power(1+eps_2,0.5)*(1-(1-2*x)*t/q2);//n=0
      coeff_CFF[0]+=TMath::Cos((TMath::Pi()-phig))*p_eff*16*x*t/q2/TMath::Power(1+eps_2,0.5)*(Kp_2*TMath::Power(2-yy,2)/q2+TMath::Power(1-(1-2*x)*t/q2,2)*(1-yy-yy*yy*eps_2/4.));//n=1 cos
      coeff_CFF[0]+=TMath::Cos(2*(TMath::Pi()-phig))*p_eff*8*K*(2-yy)*x*t/q2/TMath::Power(1+eps_2,0.5)*(1-(1-2*x)*t/q2);//n=2 cos
      coeff_CFF[1]=TMath::Sin((TMath::Pi()-phig))*p_eff*4*yy*(2-yy)*x*t/q2*(4*(1-x)*t/q2*(1+x*t/q2)+TMath::Power(1+t/q2,2)*eps_2);//n=1 sin
      coeff_CFF[1]+=-TMath::Sin(2*(TMath::Pi()-phig))*p_eff*8*K*yy*x*t/q2*(1-(1-2*x)*t/q2);//n=2 sin
    }
    if (part==2){
      //coeff_CFF A
      coeff_CFF[0]=p_eff*4*K*(2-yy)*t/q2/TMath::Power(1+eps_2,0.5)*(8-6*x+5*eps_2)*(1-t*(2-12*x*(1-x)-eps_2)/q2/(8-6*x+5*eps_2));//n=0
      coeff_CFF[0]+=TMath::Cos((TMath::Pi()-phig))*p_eff*8*t/q2/TMath::Power(1+eps_2,0.5)*(Kp_2*(1-2*x)*TMath::Power(2-yy,2)/q2+(1-(1-2*x)*t/q2)*(1-yy-yy*yy*eps_2/4.)*(4-2*x+3*eps_2+t/q2*(4*x*(1-x)+eps_2)));//n=1 cos
      coeff_CFF[0]+=TMath::Cos(2*(TMath::Pi()-phig))*p_eff*8*K*(2-yy)*t/q2*(1-x+(t-t_min)*(4*x*(1-x)+eps_2)/2./q2/TMath::Sqrt(1+eps_2));//n=2 cos
      coeff_CFF[1]=-TMath::Sin((TMath::Pi()-phig))*p_eff*8*yy*(2-yy)*(1-2*x)*t*Kp_2/q2/q2;//n=1 sin
      coeff_CFF[1]+=-TMath::Sin(2*(TMath::Pi()-phig))*p_eff*2*K*yy*t/q2*(4-4*x+2*eps_2+2*t*(4*x*(1-x)+eps_2)/q2);//n=2 sin
    }
  }

  if (helstate==2){
    //Following ones are for -+
    if (part==0){
      coeff_CFF[0]=8*(2-yy)/TMath::Power(1+eps_2,1.5)*(TMath::Power(2-yy,2)*(TMath::Sqrt(1+eps_2)-1)/(2*(1+eps_2))*Kp_2/q2+(1-yy-yy*yy*eps_2/4.)/TMath::Sqrt(1+eps_2)*(1-x-(TMath::Sqrt(1+eps_2)-1)/2.+eps_2/(2*x))*x*t*(t-t_min)/q2/q2);//n=0
      coeff_CFF[0]+=TMath::Cos((TMath::Pi()-phig))*8*K/TMath::Power(1+eps_2,1.5)*(TMath::Power(2-yy,2)*(2-TMath::Sqrt(1+eps_2))/(1+eps_2)*((1-t/q2)*(TMath::Sqrt(1+eps_2)-1+eps_2)/(2*(2-TMath::Sqrt(1+eps_2)))-x*t/q2)+2*(1-yy-yy*yy*eps_2/4.)/TMath::Sqrt(1+eps_2)*((1-TMath::Sqrt(1+eps_2)+eps_2/2.)/2./TMath::Sqrt(1+eps_2)+t/q2*(1-3*x/2.+(x+eps_2/2.)/(2*TMath::Sqrt(1+eps_2)))));//n=1 cos
      coeff_CFF[0]+=TMath::Cos(2*(TMath::Pi()-phig))*4*(2-yy)*(1-yy-yy*yy*eps_2/4.)*(1+TMath::Sqrt(1+eps_2))/TMath::Power(1+eps_2,2.5)*((2-3*x)*t/q2+x*t*t*(1-2*x+(2-2*x)/(1+TMath::Sqrt(1+eps_2)))/q2/q2+(1+t*((TMath::Sqrt(1+eps_2)+x+(1-x)*t/q2)/(1+TMath::Sqrt(1+eps_2)))/q2)*eps_2);// n=2 cos
      coeff_CFF[0]+=-TMath::Cos(3*(TMath::Pi()-phig))*8*K*(1-yy-yy*yy*eps_2/4.)*(1+TMath::Sqrt(1+eps_2)+eps_2/2.)/TMath::Power(1+eps_2,2.5)*(1+(1+TMath::Sqrt(1+eps_2)+eps_2/2./x)/(1+TMath::Sqrt(1+eps_2)+eps_2/2.)*x*t/q2);//n=3 cos
      
      coeff_CFF[1]=TMath::Sin((TMath::Pi()-phig))*4*K*(2-yy)*yy/TMath::Power(1+eps_2,2)*(1-TMath::Sqrt(1+eps_2)+2*eps_2-2*(1+(TMath::Sqrt(1+eps_2)-1)/(2*x))*x*t/q2);//n=1 sin
      coeff_CFF[1]+=TMath::Sin(2*(TMath::Pi()-phig))*2*yy*(1-yy-yy*yy*eps_2/4.)/TMath::Power(1+eps_2,2)*(1+TMath::Sqrt(1+eps_2))*(eps_2-2*(1+eps_2/2./x)*x*t/q2)*(1+t*(TMath::Sqrt(1+eps_2)-1+2*x)/(1+TMath::Sqrt(1+eps_2))/q2);//n=2 sin
     
    }
    //Coeff_CFF V
    if (part==1){
      coeff_CFF[0]=4*(2-yy)*x*t/q2/TMath::Power(1+eps_2,2.5)*(2*Kp_2/q2*(2-2*yy+yy*yy+yy*yy*eps_2/2.)-(1-(1-2*x)*t/q2)*(1-yy-yy*yy*eps_2/4.)*(TMath::Sqrt(1+eps_2)-1+t*(TMath::Sqrt(1+eps_2)+1-2*x)/q2));//n=0
      coeff_CFF[0]+=TMath::Cos((TMath::Pi()-phig))*8*K*x*t/q2/TMath::Power(1+eps_2,2.5)*(2*(1-(1-2*x)*t/q2)*(2-2*yy+yy*yy+yy*yy*eps_2/2.)+(1-yy-yy*yy*eps_2/4.)*(3-TMath::Sqrt(1+eps_2)-(3*(1-2*x)+TMath::Sqrt(1+eps_2))*t/q2));//n=1 cos (TMath::Pi()-phig)
      coeff_CFF[0]+=TMath::Cos(2*(TMath::Pi()-phig))*4*(2-yy)*(1-yy-yy*yy*eps_2/4.)*x*t/q2/TMath::Power(1+eps_2,2.5)*(4*Kp_2/q2+1+TMath::Sqrt(1+eps_2)+t*((1-2*x)*(1-2*x-TMath::Sqrt(1+eps_2))*t/q2-2+4*x+2*x*TMath::Sqrt(1+eps_2))/q2);//n=2 cos (TMath::Pi()-phig)
      coeff_CFF[0]+=TMath::Cos(3*(TMath::Pi()-phig))*8*K*(1-yy-yy*yy*eps_2/4.)*(1+TMath::Sqrt(1+eps_2))/TMath::Power(1+eps_2,2.5)*x*t/q2*(1-t*(1-2*x-TMath::Sqrt(1+eps_2))/q2/(1+TMath::Sqrt(1+eps_2)));//n=3 cos
      
      coeff_CFF[1]=TMath::Sin((TMath::Pi()-phig))*8*K*yy*(2-yy)*x*t/q2/TMath::Power(1+eps_2,2)*(1+TMath::Sqrt(1+eps_2))*(1-t*(1-2*x-TMath::Sqrt(1+eps_2))/q2/(1+TMath::Sqrt(1+eps_2)));//n=1 sin
      coeff_CFF[1]+=TMath::Sin(2*(TMath::Pi()-phig))*4*yy*(1-yy-yy*yy*eps_2/4.)*x*t/q2/TMath::Power(1+eps_2,2)*(1+TMath::Sqrt(1+eps_2))*(1-t*(1-2*x)/q2)*(1-t*(1-2*x-TMath::Sqrt(1+eps_2))/(1+TMath::Sqrt(1+eps_2))/q2);//n=2 si
    }
  //Coeff_CFF A
    if (part==2){
      coeff_CFF[0]=4*(2-yy)*t/q2/TMath::Power(1+eps_2,2)*((t-t_min)*(1-yy-yy*yy*eps_2/4.)*(2*x*x-eps_2-3*x+x*TMath::Sqrt(1+eps_2))/q2+Kp_2*(4-2*x*TMath::Power(2-yy,2)-4*yy+yy*yy-yy*yy*TMath::Power(1+eps_2,1.5))/q2/TMath::Sqrt(1+eps_2));//n=0
      coeff_CFF[0]+=TMath::Cos((TMath::Pi()-phig))*4*K*t/q2/TMath::Power(1+eps_2,2.5)*((2-2*yy+yy*yy+yy*yy*eps_2/2.)*(5-4*x+3*eps_2-TMath::Sqrt(1+eps_2)-t/q2*(1-eps_2-TMath::Sqrt(1+eps_2)-2*x*(4-4*x-TMath::Sqrt(1+eps_2))))+(1-yy-yy*yy*eps_2/4.)*(8+5*eps_2-6*x+2*x*TMath::Sqrt(1+eps_2)-t/q2*(2-eps_2+2*TMath::Sqrt(1+eps_2)-4*x*(3-3*x+TMath::Sqrt(1+eps_2)))));//n=1 cos
      coeff_CFF[0]+=TMath::Cos(2*(TMath::Pi()-phig))*16*(2-yy)*(1-yy-yy*yy*eps_2/4.)*t/TMath::Power(1+eps_2,1.5)/q2*(Kp_2*(1-2*x)/(1+eps_2)/q2-(1-x)/(4*x*(1-x)+eps_2)*(2*x*x-eps_2-3*x-x*TMath::Sqrt(1+eps_2))-(t-t_min)*(2*x*x-eps_2-3*x-x*TMath::Sqrt(1+eps_2))/(4*TMath::Sqrt(1+eps_2))/q2);//n=2 cos
      coeff_CFF[0]+=TMath::Cos(3*(TMath::Pi()-phig))*16*K*(1-yy-yy*yy*eps_2/4.)*t/q2/TMath::Power(1+eps_2,2)*(1-x+(t-t_min)*(4*x*(1-x)+eps_2)/q2/(4*TMath::Sqrt(1+eps_2)));//n=3 cos
      
      coeff_CFF[1]=TMath::Sin((TMath::Pi()-phig))*4*K*yy*(2-yy)*t/q2/TMath::Power(1+eps_2,2)*(3+2*eps_2+TMath::Sqrt(1+eps_2)-2*x-2*x*TMath::Sqrt(1+eps_2)-t*(1-2*x)*(1-2*x-TMath::Sqrt(1+eps_2))/q2);// n=1 sin
      coeff_CFF[1]+=TMath::Sin(2*(TMath::Pi()-phig))*2*yy*(1-yy-yy*yy*eps_2/4.)*t/q2/TMath::Power(1+eps_2,2)*(4-2*x+3*eps_2+t*(4*x*(1-x)+eps_2)/q2)*(1+TMath::Sqrt(1+eps_2)-t*(1-2*x-TMath::Sqrt(1+eps_2))/q2);// n=2 sin
      
    }
   }

  coeff_CFF[0]=coeff_CFF[0]*prefactor_I;
  coeff_CFF[1]=coeff_CFF[1]*prefactor_I;
  return coeff_CFF;
}


TArrayD Cprod(Int_t CFFa,Int_t CFFb,Int_t heltr,Int_t heltrf){
  TArrayD prod(2);
  prod[0]=CFF[heltr][2*CFFa]*CFF[heltrf][2*CFFb]+CFF[heltr][2*CFFa+1]*CFF[heltrf][2*CFFb+1];
  prod[1]=CFF[heltr][2*CFFa+1]*CFF[heltrf][2*CFFb]-CFF[heltr][2*CFFa]*CFF[heltrf][2*CFFb+1]; 
  //if (heltr==1&&heltrf!=1) cout<<CFF[heltr][2*CFFa+1]*CFF[heltrf][2*CFFb]<<" et "<<CFF[heltr][2*CFFa]*CFF[heltrf][2*CFFb+1]<<endl;
  return prod;
}

TArrayD Jacobian_Cprod(Int_t CFFa,Int_t CFFb,Int_t heltr,Int_t heltrf, Int_t JRow){
  TArrayD prod(2);
  prod[0]=0;
  prod[1]=0;
  Int_t JHel=JRow/8; //Find which twist
  Int_t JCFF=(JRow%8)/2; //Find H,E,Htilde or Etilde
  Int_t JReal=(JRow%8)%2; //Find Imaginary or real
  if(JHel==heltr&&JCFF==CFFa){
    if (JReal==0){//Real
      prod[0]=CFF[heltrf][2*CFFb];
      prod[1]=-CFF[heltrf][2*CFFb+1]; 
    } 
    if (JReal==1){//Imaginary
      prod[0]=CFF[heltrf][2*CFFb+1];
      prod[1]=CFF[heltrf][2*CFFb]; 
    }
  }
  if (JHel==heltrf&&JCFF==CFFb){
    if (JReal==0){
      prod[0]=CFF[heltr][2*CFFa];
      prod[1]=CFF[heltr][2*CFFa+1];
    }
    if (JReal==1){
      prod[0]=CFF[heltr][2*CFFa+1];
      prod[1]=-CFF[heltr][2*CFFa];
    }
  }
  //if (heltr==1&&heltrf!=1) cout<<CFF[heltr][2*CFFa+1]*CFF[heltrf][2*CFFb]<<" et "<<CFF[heltr][2*CFFa]*CFF[heltrf][2*CFFb+1]<<endl;
  return prod;
}

TArrayD CU_DVCS(Double_t q2, Double_t x, Double_t t,Int_t heltrans,Int_t heltransf){
  TArrayD ReIm(2);//0 is real part et 1 is imaginary part
  for (Int_t k=0;k<2;k++){
    ReIm[k]=4*(1-x)*(1+x*t/q2)/TMath::Power(2-x+x*t/q2,2)*(Cprod(0,0,heltrans,heltransf)[k]+Cprod(2,2,heltrans,heltransf)[k])+((2+t/q2)*eps_2)/TMath::Power(2-x+x*t/q2,2)*Cprod(2,2,heltrans,heltransf)[k]-t/4./M/M*Cprod(1,1,heltrans,heltransf)[k]-x*x/TMath::Power(2-x+x*t/q2,2)*(TMath::Power(1+t/q2,2)*(Cprod(0,1,heltrans,heltransf)[k]+Cprod(1,0,heltrans,heltransf)[k]+Cprod(1,1,heltrans,heltransf)[k])+Cprod(2,3,heltrans,heltransf)[k]+Cprod(3,2,heltrans,heltransf)[k]+t/4./M/M*Cprod(3,3,heltrans,heltransf)[k]);
  }
 
  return ReIm;

}

TArrayD Jacobian_CU_DVCS(Double_t q2, Double_t x, Double_t t,Int_t heltrans,Int_t heltransf, Int_t JRow){
  TArrayD ReIm(2);//0 is real part et 1 is imaginary part
  for (Int_t k=0;k<2;k++){
    ReIm[k]=4*(1-x)*(1+x*t/q2)/TMath::Power(2-x+x*t/q2,2)*(Jacobian_Cprod(0,0,heltrans,heltransf,JRow)[k]+Jacobian_Cprod(2,2,heltrans,heltransf,JRow)[k])+((2+t/q2)*eps_2)/TMath::Power(2-x+x*t/q2,2)*Jacobian_Cprod(2,2,heltrans,heltransf,JRow)[k]-t/4./M/M*Jacobian_Cprod(1,1,heltrans,heltransf,JRow)[k]-x*x/TMath::Power(2-x+x*t/q2,2)*(TMath::Power(1+t/q2,2)*(Jacobian_Cprod(0,1,heltrans,heltransf,JRow)[k]+Jacobian_Cprod(1,0,heltrans,heltransf,JRow)[k]+Jacobian_Cprod(1,1,heltrans,heltransf,JRow)[k])+Jacobian_Cprod(2,3,heltrans,heltransf,JRow)[k]+Jacobian_Cprod(3,2,heltrans,heltransf,JRow)[k]+t/4./M/M*Jacobian_Cprod(3,3,heltrans,heltransf,JRow)[k]);
  }
 
  return ReIm;

}

TArrayD CU_I(Double_t q2, Double_t x, Double_t t, Double_t kz, Double_t phi,Int_t heltrans,Int_t comp){
  TArrayD ReIm(2);//0 is real part et 1 is imaginary part
  ReIm[0]=0;ReIm[1]=0;
  TArrayD kinfac=BH_DVCS_interference(q2,x,t,kz,phi,heltrans,comp);
  if (comp==0){
    for (Int_t k=0;k<2;k++){
      ReIm[k]+=kinfac[k]*(F1p(t)*CFF[heltrans][k]-t/4./M/M*F2p(t)*CFF[heltrans][2+k]+x/(2-x+x*t/q2)*(F2p(t)+F1p(t))*CFF[heltrans][4+k]);
    }
  }
  //kinfac=BH_DVCS_interference(1);
  if (comp==1){
    for (Int_t k=0;k<2;k++){
      ReIm[k]+=kinfac[k]*x/(2-x+x*t/q2)*(F2p(t)+F1p(t))*(CFF[heltrans][k]+CFF[heltrans][k+2]);
    } 
  }
  //kinfac=BH_DVCS_interference(2);
  if (comp==2){
    for (Int_t k=0;k<2;k++){
      ReIm[k]+=kinfac[k]*x/(2-x+x*t/q2)*(F2p(t)+F1p(t))*CFF[heltrans][4+k];
    }
  }
  return ReIm;

}

TArrayD Jacobian_CU_I(Double_t q2, Double_t x, Double_t t, Double_t kz, Double_t phi,Int_t heltrans,Int_t comp, Int_t JRow){
  TArrayD ReIm(2);//0 is real part et 1 is imaginary part
  ReIm[0]=0;ReIm[1]=0;
  TArrayD kinfac=BH_DVCS_interference(q2,x,t,kz,phi,heltrans,comp);
  Int_t JHel=JRow/8; //Find which twist
  Int_t JCFF=(JRow%8)/2; //Find H,E,Htilde or Etilde
  Int_t JReal=(JRow%8)%2; //Find Imaginary or real
  if (JHel==heltrans){
    if (comp==0){
      if (JCFF==0) ReIm[JReal]+=kinfac[JReal]*F1p(t);
      if (JCFF==1) ReIm[JReal]+=kinfac[JReal]*(-t/4./M/M*F2p(t));
      if (JCFF==2) ReIm[JReal]+=kinfac[JReal]*x/(2-x+x*t/q2)*(F2p(t)+F1p(t));
    }
    //kinfac=BH_DVCS_interference(1);
    if (comp==1){
      if (JCFF==0||JCFF==1) ReIm[JReal]+=kinfac[JReal]*x/(2-x+x*t/q2)*(F2p(t)+F1p(t));
    }
    //kinfac=BH_DVCS_interference(2);
    if (comp==2){
      if (JCFF==2) ReIm[JReal]+=kinfac[JReal]*x/(2-x+x*t/q2)*(F2p(t)+F1p(t));
    }
  }
  return ReIm;

}


void XS(Double_t q2, Double_t x, Double_t t, Double_t kz, Double_t phi,Bool_t TMC){
  yy=q2/2./M/x/kz;
  //phi=TMath::Pi()-phi; // conversion phi_trento to phi_BMK
  //From BM hotfix, DVCS unpolarized target case
  TArrayD coeff(4);
  x_2=TMath::Power(x,2);
  eps=2*x*M/TMath::Sqrt(q2);
  eps_2=TMath::Power(eps,2);
  t_min=-q2*(2*(1-x)*(1-TMath::Sqrt(1+eps_2))+eps_2)/(4*x*(1-x)+eps_2);
  t_max=-q2*(2*(1-x)*(1+TMath::Sqrt(1+eps_2))+eps_2)/(4*x*(1-x)+eps_2);
  K_2=-t*(1-x)*(1-yy-TMath::Power(yy*eps/2.,2))*(1-t_min/t)*(TMath::Sqrt(1+eps_2)+(4*x*(1-x)+eps_2)/(4*(1-x))*(t-t_min)/q2)/q2;
  K=TMath::Sqrt(K_2);
  KT=TMath::Sqrt((1-x)*x+eps_2/4.)*TMath::Sqrt((t_min-t)*(t-t_max)/q2);
  polv=(1-yy-eps_2*yy*yy/4.)/(1-yy+yy*yy/2.+eps_2*yy*yy/4.);
  prefactor_DVCS=hbar_c*TMath::Power(alpha,3)*x/(TMath::Pi()*8*q2*TMath::Sqrt(1+eps_2))/q2/yy;// Prefactor DVCS
  Jacobian=1/(2*x*M*kz); 
  prefactor_DVCS=Jacobian*prefactor_DVCS;
 

  Double_t XO=TMath::Sqrt(2*q2)*KT/(q2+t)/TMath::Sqrt(1+eps_2);// gamma in BMP is epsilon in BMJ
  Double_t X=(q2-t+2*x*t)/TMath::Sqrt(1+eps_2)/(q2+t)-1;
 
  for (int i=0;i<8;i++){
    if (!TMC){
      // if no TMC we want to extract BMP CFFs
      /*CFF[0][i]=((1+TMath::Sqrt(1+eps_2))/(2.*TMath::Sqrt(1+eps_2))+((1-x)*x*x*(4*M_2-t)*(1+t/q2))/(q2*TMath::Sqrt(1+eps_2)*TMath::Power(2-x+x*t/q2,2)))*CFF_BMP[0][i]+(1-TMath::Sqrt(1+eps_2))/(2.*TMath::Sqrt(1+eps_2))*2.*KT*KT/M_2/TMath::Power(2-x+x*t/q2,2)*CFF_BMP[2][i]+4*x*x*KT*KT/q2/TMath::Sqrt(1+eps_2)/TMath::Power(2-x+x*t/q2,3)*CFF_BMP[1][i];
      CFF[1][i]=-TMath::Sqrt(2)*KT/TMath::Sqrt(q2)/(2-x+x*t/q2)/TMath::Sqrt(1+eps_2)*(2*x*(1+2*x*x*(4*M_2-t)/q2/(2-x+x*t/q2))*CFF_BMP[1][i]+x*(1+2*x*(4*M_2-t)/q2/(2-x+x*t/q2))*CFF_BMP[0][i]+(4.*x*x*M_2-t*(2*x-eps_2))/(2*M_2*(2-x+x*t/q2))*CFF_BMP[2][i]);
      CFF[2][i]=((1-TMath::Sqrt(1+eps_2))/(2.*TMath::Sqrt(1+eps_2))+((1-x)*x*x*(4*M_2-t)*(1+t/q2))/(q2*TMath::Sqrt(1+eps_2)*TMath::Power(2-x+x*t/q2,2)))*CFF_BMP[0][i]+(1+TMath::Sqrt(1+eps_2))/(2.*TMath::Sqrt(1+eps_2))*2.*KT*KT/M_2/TMath::Power(2-x+x*t/q2,2)*CFF_BMP[2][i]+4*x*x*KT*KT/q2/TMath::Sqrt(1+eps_2)/TMath::Power(2-x+x*t/q2,3)*CFF_BMP[1][i];*/
      CFF[0][i]=CFF_BMP[0][i];
      CFF[1][i]=CFF_BMP[1][i];
      CFF[2][i]=CFF_BMP[2][i];
    }
    if (TMC){
      CFF[0][i]=CFF_BMP[0][i]+X*(CFF_BMP[0][i]+CFF_BMP[2][i])/2.-XO*CFF_BMP[1][i];
      CFF[1][i]=-(1+X)*CFF_BMP[1][i]+(CFF_BMP[0][i]+CFF_BMP[2][i])*XO;
      CFF[2][i]=CFF_BMP[2][i]+X*(CFF_BMP[0][i]+CFF_BMP[2][i])/2.-XO*CFF_BMP[1][i];
      for (int j=0;j<24;j++){
	//All the formulaes have been derived wrt BMJ CFF. In case of TMC, we need to apply a jacobian since we change our parameters basis.
	if (j/8==0) {
	  JacobianTMC[j][j%8]=1+X/2.;
	  JacobianTMC[j][j%8+8]=XO;
	  JacobianTMC[j][j%8+16]=X/2.;
	}
	if (j/8==1) {
	  JacobianTMC[j][j%8]=-XO;
	  JacobianTMC[j][j%8+8]=-(1+X);
	  JacobianTMC[j][j%8+16]=-XO;
	}
	if (j/8==2) {
	  JacobianTMC[j][j%8]=X/2.;
	  JacobianTMC[j][j%8+8]=XO;
	  JacobianTMC[j][j%8+16]=1+X/2.;
	}
      }
    }
  }

  result=0;
  result_pol=0;
  DVCS=0;
  DVCS_pol=0;
  inter=0;
  inter_pol=0;
 
  // Compute BH part;
  BH_cont=BH(q2, x, t, kz, phi);
  result+=BH_cont;

  // Compute the DVCS2 part;
  DVCS+=prefactor_DVCS*2*(2-2*yy+TMath::Power(yy,2)+eps_2*TMath::Power(yy,2)/2.)/(1+eps_2)*(CU_DVCS(q2, x, t,0,0)[0]+CU_DVCS(q2, x, t,2,2)[0]+2*polv*CU_DVCS(q2, x, t,1,1)[0]);
  DVCS+=prefactor_DVCS*(2-yy)*(4*TMath::Sqrt(2*(1-yy-eps_2*yy*yy/4.)))*(CU_DVCS(q2, x, t,1,0)[0]+CU_DVCS(q2, x, t,1,2)[0])*TMath::Cos(TMath::Pi()-phi)/(1+eps_2);
  DVCS+=prefactor_DVCS*8*(1-yy-eps_2*yy*yy/4.)*CU_DVCS(q2, x, t,2,0)[0]*TMath::Cos(2*(TMath::Pi()-phi))/(1+eps_2);
  result+=DVCS;
 
  DVCS_pol+=prefactor_DVCS*(yy*TMath::Sqrt(1+eps_2))*(4*TMath::Sqrt(2*(1-yy-eps_2*yy*yy/4.)))*(CU_DVCS(q2, x, t,1,0)[1]+CU_DVCS(q2, x, t,1,2)[1])*TMath::Sin(TMath::Pi()-phi)/(1+eps_2);
  result_pol-=DVCS_pol;//minus because there is a minus in BMJ

  //Compute Interference part
  for (Int_t ii=0;ii<3;ii++){
    for (Int_t kk=0;kk<3;kk++){
      inter+=CU_I(q2, x, t, kz, phi,ii,kk)[0];
      inter_pol+=CU_I(q2, x, t, kz,phi,ii,kk)[1];
    }
  }
  result+=inter;
  result_pol+=inter_pol;
  

  ///////*****************************///////////////////////Evaluation of the error on interference and DVCS contribution
  if (Eval_error){ 
    for (int kk=0;kk<24;kk++){
      Jacobian_DVCS_BMJ[kk]=0.0;
      Jacobian_DVCS_pol_BMJ[kk]=0.0;
      Jacobian_DVCS_BMP[kk]=0.0;
      Jacobian_DVCS_pol_BMP[kk]=0.0;
      
      Jacobian_Inter_BMJ[kk]=0.0;
      Jacobian_Inter_pol_BMJ[kk]=0.0;
      Jacobian_Inter_BMP[kk]=0.0;
      Jacobian_Inter_pol_BMP[kk]=0.0;
    }
    err_DVCS=0.0;
    err_DVCS_pol=0.0;
    err_Inter=0.0;
    err_Inter_pol=0.0;
    for (int kk=0;kk<24;kk++){
      //Gradient for DVCS contribution
      Jacobian_DVCS_BMJ[kk]+=prefactor_DVCS*2*(2-2*yy+TMath::Power(yy,2)+eps_2*TMath::Power(yy,2)/2.)/(1+eps_2)*(Jacobian_CU_DVCS(q2, x, t,0,0,kk)[0]+Jacobian_CU_DVCS(q2, x, t,2,2,kk)[0]+2*polv*Jacobian_CU_DVCS(q2, x, t,1,1,kk)[0]); 
      Jacobian_DVCS_BMJ[kk]+=prefactor_DVCS*(2-yy)*(4*TMath::Sqrt(2*(1-yy-eps_2*yy*yy/4.)))*(Jacobian_CU_DVCS(q2, x, t,1,0,kk)[0]+Jacobian_CU_DVCS(q2, x, t,1,2,kk)[0])*TMath::Cos(TMath::Pi()-phi)/(1+eps_2);
      Jacobian_DVCS_BMJ[kk]+=prefactor_DVCS*8*(1-yy-eps_2*yy*yy/4.)*Jacobian_CU_DVCS(q2, x, t,2,0,kk)[0]*TMath::Cos(2*(TMath::Pi()-phi))/(1+eps_2);
    
      Jacobian_DVCS_pol_BMJ[kk]+=prefactor_DVCS*(yy*TMath::Sqrt(1+eps_2))*(4*TMath::Sqrt(2*(1-yy-eps_2*yy*yy/4.)))*(Jacobian_CU_DVCS(q2, x, t,1,0,kk)[1]+Jacobian_CU_DVCS(q2, x, t,1,2,kk)[1])*TMath::Sin(TMath::Pi()-phi)/(1+eps_2);

      //Gradient for Inter contribution
      for (Int_t ii=0;ii<3;ii++){
	for (Int_t jj=0;jj<3;jj++){
	  Jacobian_Inter_BMJ[kk]+=Jacobian_CU_I(q2, x, t, kz, phi,ii,jj,kk)[0];
	  Jacobian_Inter_pol_BMJ[kk]+=Jacobian_CU_I(q2, x, t, kz,phi,ii,jj,kk)[1];
	}
      }
    }
 
    for (int  kk=0;kk<24;kk++){
      if (TMC){
      	for (int aa=0;aa<24;aa++){
	  Jacobian_DVCS_BMP[kk]+=JacobianTMC[kk][aa]*Jacobian_DVCS_BMJ[aa];
	  Jacobian_DVCS_pol_BMP[kk]+=JacobianTMC[kk][aa]*Jacobian_DVCS_pol_BMJ[aa];
	  Jacobian_Inter_BMP[kk]+=JacobianTMC[kk][aa]*Jacobian_Inter_BMJ[aa];
	  Jacobian_Inter_pol_BMP[kk]+=JacobianTMC[kk][aa]*Jacobian_Inter_pol_BMJ[aa];
	}
      }
    
      if (!TMC){
	Jacobian_DVCS_BMP[kk]=Jacobian_DVCS_BMJ[kk];
	Jacobian_DVCS_pol_BMP[kk]=Jacobian_DVCS_pol_BMJ[kk];
	Jacobian_Inter_BMP[kk]=Jacobian_Inter_BMJ[kk];
	Jacobian_Inter_pol_BMP[kk]=Jacobian_Inter_pol_BMJ[kk];
      }
    }

    for (int kk=0;kk<24;kk++){
      for (int aa=0;aa<24;aa++){
	err_DVCS+=Jacobian_DVCS_BMP[kk]*covariance[kk][aa]*Jacobian_DVCS_BMP[aa]; //tJ*cov*J
	err_DVCS_pol+=Jacobian_DVCS_pol_BMP[kk]*covariance[kk][aa]*Jacobian_DVCS_pol_BMP[aa];
	err_Inter+=Jacobian_Inter_BMP[kk]*covariance[kk][aa]*Jacobian_Inter_BMP[aa]; //tJ*cov*J
	err_Inter_pol+=Jacobian_Inter_pol_BMP[kk]*covariance[kk][aa]*Jacobian_Inter_pol_BMP[aa];
      }
    }
    err_DVCS=TMath::Sqrt(err_DVCS);
    err_DVCS_pol=TMath::Sqrt(err_DVCS_pol);
    err_Inter=TMath::Sqrt(err_Inter);
    err_Inter_pol=TMath::Sqrt(err_Inter_pol);
  }

    return;
}

Double_t chi0(Double_t q2, Double_t x, Double_t t){
 
  //phi=TMath::Pi()-phi; // conversion phi_trento to phi_BMK
  //From BM hotfix, DVCS unpolarized target case
  TArrayD coeff(4);
  x_2=TMath::Power(x,2);
  eps=2*x*M/TMath::Sqrt(q2);
  eps_2=TMath::Power(eps,2);
  t_min=-q2*(2*(1-x)*(1-TMath::Sqrt(1+eps_2))+eps_2)/(4*x*(1-x)+eps_2);
  t_max=-q2*(2*(1-x)*(1+TMath::Sqrt(1+eps_2))+eps_2)/(4*x*(1-x)+eps_2);
  KT=TMath::Sqrt((1-x)*x+eps_2/4.)*TMath::Sqrt((t_min-t)*(t-t_max)/q2);
 

  Double_t XO=TMath::Sqrt(2*q2)*KT/(q2+t)/TMath::Sqrt(1+eps_2);// gamma in BMP is epsilon in BMJ

  return XO;
} 

Double_t chi(Double_t q2, Double_t x, Double_t t){
 
  //phi=TMath::Pi()-phi; // conversion phi_trento to phi_BMK
  //From BM hotfix, DVCS unpolarized target case
  TArrayD coeff(4);
  x_2=TMath::Power(x,2);
  eps=2*x*M/TMath::Sqrt(q2);
  eps_2=TMath::Power(eps,2);
  t_min=-q2*(2*(1-x)*(1-TMath::Sqrt(1+eps_2))+eps_2)/(4*x*(1-x)+eps_2);
  t_max=-q2*(2*(1-x)*(1+TMath::Sqrt(1+eps_2))+eps_2)/(4*x*(1-x)+eps_2);
  KT=TMath::Sqrt((1-x)*x+eps_2/4.)*TMath::Sqrt((t_min-t)*(t-t_max)/q2);
 
  Double_t X=(q2-t+2*x*t)/TMath::Sqrt(1+eps_2)/(q2+t)-1;

 // gamma in BMP is epsilon in BMJ

  return X;
} 
