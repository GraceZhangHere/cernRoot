#ifndef BMJ_TMC_H
#define BMJ_TMC_H

#include "TArrayD.h"
#include "TMath.h"

Double_t alpha=7.297352533e-3;
Double_t hbar_c=0.38937966e6;
Double_t e=1.6e-19;
Double_t M=0.938271998;
Double_t M_2=M*M;

Double_t prefactor_BH;
Double_t polv;
Double_t prefactor_DVCS;
Double_t prefactor_I;
Double_t Jacobian;
Double_t K;Double_t K_2;
Double_t Kp;Double_t Kp_2;
Double_t KT;
Double_t J;
Double_t yy;
Double_t P1; Double_t P2;
Double_t t_min,t_max;
Double_t F1_2;
Double_t F2_2;
Double_t F_sum_2;

Double_t x_2;
Double_t eps_2;
Double_t eps;

//Cross section result
Double_t result;
Double_t result_pol;

//If one want the details
Double_t DVCS, DVCS_pol, inter, inter_pol, BH_cont;
Double_t err_DVCS, err_Inter, err_DVCS_pol, err_Inter_pol;
// Definition of the CFF 
// 0=++, 1=0+, 2=-+
// ReH 0, ImH 1, ReE 2, ImE 3, ReHtilde 4,...
Double_t CFF[3][8];
Double_t CFF_BMP[3][8];

Bool_t Eval_error=false;
Double_t covariance[24][24];
Double_t JacobianTMC[24][24];
Double_t Jacobian_DVCS_BMJ[24];
Double_t Jacobian_Inter_BMJ[24];
Double_t Jacobian_DVCS_pol_BMJ[24];
Double_t Jacobian_Inter_pol_BMJ[24];
Double_t Jacobian_DVCS_BMP[24];
Double_t Jacobian_Inter_BMP[24];
Double_t Jacobian_DVCS_pol_BMP[24];
Double_t Jacobian_Inter_pol_BMP[24];

//Apply or not TMC at leading twist
Double_t tmin(Double_t q2,Double_t x);

Double_t GMp(Double_t q2);
Double_t GEp(Double_t q2);

Double_t F1p(Double_t q2);
Double_t F2p(Double_t q2);

Double_t KellyE(Double_t q2);
Double_t KellyM(Double_t q2);

Double_t BH(Double_t q2, Double_t x, Double_t t, Double_t kz, Double_t phig);
TArrayD BH_DVCS_interference(Double_t q2, Double_t x, Double_t t, Double_t kz, Double_t phig, Int_t helstate, Int_t part);
TArrayD Cprod(Int_t CFFa,Int_t CFFb,Int_t heltr,Int_t heltrf);
TArrayD Jacobian_Cprod(Int_t CFFa,Int_t CFFb,Int_t heltr,Int_t heltrf, Int_t JRow);
TArrayD CU_DVCS(Double_t q2, Double_t x, Double_t t,Int_t heltrans,Int_t heltransf);
TArrayD Jacobian_CU_DVCS(Double_t q2, Double_t x, Double_t t,Int_t heltrans,Int_t heltransf, Int_t JRow);
TArrayD CU_I(Double_t q2, Double_t x, Double_t t, Double_t kz, Double_t phi,Int_t heltrans,Int_t comp);
TArrayD Jacobian_CU_I(Double_t q2, Double_t x, Double_t t, Double_t kz, Double_t phi,Int_t heltrans,Int_t comp, Int_t JRow);
void XS(Double_t q2, Double_t x, Double_t t, Double_t kz, Double_t phi, Bool_t TMC);
Double_t chi0(Double_t q2, Double_t x, Double_t t);
Double_t chi(Double_t q2, Double_t x, Double_t t);
#endif
