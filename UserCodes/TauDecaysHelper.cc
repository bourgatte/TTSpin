#include "TauDecaysHelper.h"
#include <iostream>
TauDecaysHelper::TauDecaysHelper():
  mrho(0.773),
  mpi(0.13957018),
  mtau(1.776),
  ma1(1.251),
  debug(false)
{ }
TauDecaysHelper::TauDecaysHelper(vector<TLorentzVector> TauAndProd, string type): 
  mrho(0.773),
  mpi(0.13957018),
  mtau(1.776),
  ma1(1.251),
  debug(false)
 {
  if(type=="lepton" || type =="pion"){
    if(TauAndProd.size()!=2){ 
      std::cout<<" Warning!! Size of input vector  !=2  !! "<< " type:  "<<type<<std::endl;
    }
  }
  if(type=="rho"){
    if(TauAndProd.size()!=3){
      std::cout<<" Warning!! Size of input vector  !=3  !! "<< " type:  "<<type<<std::endl;
    }
  }

  Configure(TauAndProd,  type);
}


void 
TauDecaysHelper::Configure(vector<TLorentzVector> TauAndProd, string type){
  TauLV = TauAndProd.at(0);
  type_=type;
  if(type_ == "lepton" || type_=="pion"){
    ProductLV= TauAndProd.at(1);
  }
  if(type_ == "rho" ){
    TauRhoPi  = TauAndProd.at(1);
    TauRhoPi0= TauAndProd.at(2);
    ProductLV = TauRhoPi+TauRhoPi0;
    DPF_TauRhoPi = Boost(TauRhoPi,ProductLV);
    DPF_TauRhoPi0 =  Boost(TauRhoPi0,ProductLV);
  }
  InvisibleLV = TauLV - ProductLV;

  DPF_TauLV=  Boost(TauLV,ProductLV);
  DPF_InvisibleLV=  Boost(InvisibleLV,ProductLV);
  //  std::cout<<"TauDecaysHelper  Setup Pair  type 1 "<< type_ <<std::endl; ProductLV.Print();
}
bool  
TauDecaysHelper::isConfigured(){
  if(TauLV.E()!=0 && ProductLV.E()!=0) return true; return false;
}
TauDecaysHelper::~TauDecaysHelper(){
}
TLorentzVector 
TauDecaysHelper::Boost(TLorentzVector pB, TLorentzVector frame){
   TMatrixT<double> transform(4,4);
   TMatrixT<double> result(4,1);
   TVectorT<double> vec(4); 
   TVector3 b;
   if(frame.Vect().Mag()==0){ std::cout<<"TDH Boost is not set, perfrom calculation in the Lab Frame   "<<std::endl; return pB;}
    if(frame.E()==0){ std::cout<<" Caution: Please check that you perform boost correctly!  " <<std::endl; return pB;} 
   else   b=frame.Vect()*(1/frame.E());
   vec(0)  = pB.E();    vec(1)  = pB.Px();
   vec(2)  = pB.Py();  vec(3)  = pB.Pz();
   double gamma  = 1/sqrt( 1 - b.Mag2());
   transform(0,0)=gamma; transform(0,1) =- gamma*b.X() ;  transform(0,2) =  - gamma*b.Y();  transform(0,3) = - gamma*b.Z(); 
   transform(1,0)=-gamma*b.X(); transform(1,1) =(1+ (gamma-1)*b.X()*b.X()/b.Mag2()) ;  transform(1,2) = ((gamma-1)*b.X()*b.Y()/b.Mag2());  transform(1,3) = ((gamma-1)*b.X()*b.Z()/b.Mag2());
   transform(2,0)=-gamma*b.Y(); transform(2,1) = ((gamma-1)*b.Y()*b.X()/b.Mag2());  transform(2,2) = (1 + (gamma-1)*b.Y()*b.Y()/b.Mag2());  transform(2,3) =  ((gamma-1)*b.Y()*b.Z()/b.Mag2()); 
   transform(3,0)=-gamma*b.Z(); transform(3,1) =((gamma-1)*b.Z()*b.X()/b.Mag2()) ;  transform(3,2) = ((gamma-1)*b.Z()*b.Y()/b.Mag2());  transform(3,3) = (1 + (gamma-1)*b.Z()*b.Z()/b.Mag2()); 
   result=transform*convertToMatrix(vec);
   return TLorentzVector(result(1,0), result(2,0) ,result(3,0), result(0,0));
}

double
TauDecaysHelper::getOmega(){
  double omega=-999;
  if(type_=="pion" || type_=="lepton"){
    omega = 2*ProductLV.E()/TauLV.E() - 1;
  }
  if(type_=="rho"){
    double QQ =  ProductLV.M2();
    double Be = 0.5*(3*getCosbetaRho() -1);
    double Ps = 0.5*(3*getUltrarel_cospsiLF() -1);
    double RR = mtau*mtau/QQ;
    omega = ((-2 + RR + 2*(1+RR)*Ps*Be)*getCosthetaRho() + 3*sqrt(RR)*Be*getSinthetaRho()*2*getUltrarel_cospsiLF()*getSinpsiLF())  /  ( 2 +RR - 2*(1-RR)*Ps*Be);
    if(  isinf(fabs(omega)) ||  isnan(fabs(omega))) omega  = -999.;
  }
  return omega;
}

double
TauDecaysHelper::getOmegabar(){
  double omega=-999;
  if(type_=="pion" || type_=="lepton"){
    omega = 2*ProductLV.E()/TauLV.E() - 1;
  }
  if(type_=="rho"){
    double QQ =  ProductLV.M2();
    double Be = 0.5*(3*getCosbetaRho() -1);
    double Ps = 0.5*(3*getUltrarel_cospsiLF() -1);
    double RR = mtau*mtau/QQ;
    omega = (RR*getCosthetaRho() - sqrt(RR)*getSinthetaRho()*2* getSinbetaRho()*getCosbetaRho()*DPF_cosalpha() -   getCosthetaRho()* getSinbetaRho()*getSinbetaRho()*(1+RR) )  /  ( RR + (1-RR)*getSinbetaRho()*getSinbetaRho());
    if(  isinf(fabs(omega)) ||  isnan(fabs(omega))) omega  = -999.;
  }
  return omega;
}


double
TauDecaysHelper::getCosbetaRho(){
  double cb=-999;
  if(type_!="rho"){std::cout<<"Can not return cos beta, decay type is not rho!!  but: "<< type_ <<std::endl; return cb;}
  cb = (mrho/sqrt(mrho*mrho - 4*mpi*mpi))  *   (TauRhoPi.E() -   TauRhoPi0.E())/(TauRhoPi.E() +  TauRhoPi0.E()) ;
  if(fabs(cb) - 1 >0  && fabs(cb) < 1.01 &&  cb> 0  )  cb= 1;
  if(fabs(cb) - 1 >0  && fabs(cb) < 1.01 &&  cb< 0  )  cb=-1;

  if(fabs(cb) > 1 ){if(debug){std::cout<<"Warning! Cos beta > 1:  "<<cb <<std::endl;  }}
  return cb;
}



double
TauDecaysHelper::getSinbetaRho(){
  double sb=-999;
  if(type_!="rho"){std::cout<<"Can not return sin beta, decay type is not rho!!  but: "<< type_ <<std::endl; return sb;}
  if(fabs(getCosbetaRho()) > 1 ){if(debug){std::cout<<"Warning! Cos beta > 1:  "<< getCosbetaRho()<<std::endl;  }return 0;}
  sb = sqrt(1- getCosbetaRho()*getCosbetaRho());
  return sb;
}



double
TauDecaysHelper::getCosthetaRho(){
  double ct=-999;
  if(type_!="rho"){std::cout<<"Can not return cos theta, decay type is not rho!!  but: "<< type_ <<std::endl; return ct;}
  double QQ = ProductLV.M2();
  double x = ProductLV.E()/TauLV.E();
  double s = 4*TauLV.E()*TauLV.E();
  if(fabs(ct) - 1 >0  && fabs(ct) < 1.01 &&  ct> 0  ) ct = 1;
  if(fabs(ct) - 1 >0  && fabs(ct) < 1.01 &&  ct< 0  )  ct=-1;

  if( 1 - 4*mtau*mtau/s  <= 0 ){std::cout<<"Warning! In costheta root square <=0! return -999"<<std::endl; return ct;}
  ct= (2*x*mtau*mtau - mtau*mtau - QQ)/( (mtau*mtau - QQ)*sqrt(1 - 4*mtau*mtau/s) );
 
  if(fabs(ct) > 1 ){if(debug){std::cout<<"Warning! Cos theta > 1:  "<<ct<<std::endl; }}
  return ct;
}
double
TauDecaysHelper::getSinthetaRho(){
  double st=-999;
  if(type_!="rho"){std::cout<<"Can not return cos theta, decay type is not rho!!  but: "<< type_ <<std::endl; return st;}
   st = sqrt(1- getCosthetaRho()*getCosthetaRho());
   if(fabs(getCosthetaRho()) > 1 ){if(debug){std::cout<<"Warning! Cos theta > 1"<<std::endl; }}
  return st;
}


double 
TauDecaysHelper::getUltrarel_cospsiLF(){
  double cos=-999;
  if(type_!="rho"){std::cout<<"Can not return cos theta, decay type is not rho!!  but: "<< type_ <<std::endl; return cos;}
  double QQ = ProductLV.M2();
  cos = (getCosthetaRho()*(mtau*mtau  + QQ)   + (mtau*mtau  - QQ))/(getCosthetaRho()*(mtau*mtau  - QQ)   + (mtau*mtau  + QQ));
  if(fabs(cos) - 1 >0  && fabs(cos) < 1.01 && cos > 0  ) cos = 1;
  if(fabs(cos) - 1 >0  && fabs(cos) < 1.01 && cos < 0  ) cos =-1;
  if(fabs(cos) > 1 )if(debug){std::cout<<"Warning! Cos psi > 1:  "<<cos<<std::endl; }
  return cos;
}
double 
TauDecaysHelper::getSinpsiLF(){
  double sin = -999;
  if(type_!="rho"){std::cout<<"Can not return cos theta, decay type is not rho!!  but: "<< type_ <<std::endl; return sin;}
  sin = sqrt(1 - getUltrarel_cospsiLF()*getUltrarel_cospsiLF());
  if(getUltrarel_cospsiLF()*getUltrarel_cospsiLF() > 1  ){if(debug){std::cout<<"Warning! getultrarel_cospsiLF > 1"<<std::endl;}}
  return    sin;
}


TMatrixT<double> TauDecaysHelper::convertToMatrix(TVectorT<double> V){
  TMatrixT<double> M(V.GetNrows(),1);
  for(int i=0; i < M.GetNrows(); i++){
    M(i,0)=V(i);
  } return M;
}
 
 TLorentzVector
 TauDecaysHelper::sLV(){
   if(type_!="rho"){std::cout<<"Warninng in sLV(). This is not rho decay!"<< type_ <<std::endl; }
   double QQ = ProductLV.M2();
   double l0 = 0.5*(mtau*mtau + QQ)/sqrt(QQ);
   double l   = 0.5*(mtau*mtau  - QQ)/sqrt(QQ);
   return TLorentzVector(getSinpsiLF(),0,-l0*getCosthetaRho()/mtau,-l*getCosthetaRho()/mtau);
 }


 TVector3
 TauDecaysHelper::nPerp(){
   if(type_!="rho"){std::cout<<"Warninng in nPerp(). This is not rho decay!  "<< type_ <<std::endl; return TVector3(0,0,0);}
   return   DPF_TauRhoPi.Vect()*(1/DPF_TauRhoPi.Vect().Mag());


 }

 TVector3
 TauDecaysHelper::ns(){
  if(type_!="rho"){std::cout<<"Warninng in ns(). This is not rho decay!  "<< type_ <<std::endl; return TVector3(0,0,0);}
   return   sLV().Vect()*(1/sLV().Vect().Mag());
 }
 TVector3
 TauDecaysHelper::nL(){
  if(type_!="rho"){std::cout<<"Warninng in nL(). This is not rho decay!  "<< type_ <<std::endl; return TVector3(0,0,0);}
   return   -ProductLV.Vect()*(1/ProductLV.Vect().Mag());
 }
 TVector3
 TauDecaysHelper::nT(){
  if(type_!="rho"){std::cout<<"Warninng in nT(). This is not rho decay!  "<< type_ <<std::endl; return TVector3(0,0,0);}
   return   DPF_TauLV.Vect()*(1/DPF_TauLV.Vect().Mag());
 }


 double  
 TauDecaysHelper::DPF_cosalpha(){
   if(type_!="rho"){std::cout<<"Warninng in DPF_cosalpha(). This is not rho decay!  "<< type_ <<std::endl; return 0;}
    TVector3 nTCrossns  = nT().Cross(ns());
    TVector3 nTCrossnPerp  = nT().Cross(nPerp());

    if(nTCrossns.Mag() ==0 || nTCrossnPerp.Mag() ==0){if(debug){std::cout<<" Can not compute cos alpha, one denominator is 0, return DPF cos alpha =0  "<< std::endl; }return 0;}
   return nTCrossns.Dot(nTCrossnPerp)/nTCrossns.Mag()/nTCrossnPerp.Mag();
 }
 double  
 TauDecaysHelper::DPF_sinalpha(){
  if(type_!="rho"){std::cout<<"Warninng in DPF_sinalpha(). This is not rho decay!  "<< type_ <<std::endl; return 0;}
    TVector3 nTCrossns  = nT().Cross(ns());
    TVector3 nTCrossnPerp  = nT().Cross(nPerp());

    if(nTCrossns.Mag() ==0 || nTCrossnPerp.Mag() ==0){if(debug){std::cout<<" Can not compute sin alpha, one denominator is 0, return DPF sin alpha =0  "<< std::endl; }return 0;}
   return -ns().Dot(nTCrossnPerp)/nTCrossns.Mag()/nTCrossnPerp.Mag();

 }
