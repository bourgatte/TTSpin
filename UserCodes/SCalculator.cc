#include "SCalculator.h"
#include <iostream>

SCalculator::SCalculator(){
}

SCalculator::SCalculator(vector<TLorentzVector> TauAndProd)
{
  if(TauAndProd.size()!=2){
    std::cout<<" Warning!! Size of rho input vector != 2 !! "<<std::endl;
  }
  TLorentzVector LabFrame = TauAndProd.at(1);
  Setup(TauAndProd, LabFrame);
}

void 
SCalculator::Setup(vector<TLorentzVector> TauAndProd, TLorentzVector ReferenceFrame){
   debug=false;
   for(unsigned int i=0; i<TauAndProd.size(); i++){
     TauAndProd_LF.push_back(TauAndProd.at(i));
   }



 

}

void 
SCalculator::Configure(vector<TLorentzVector> TauAndProd){
  if(TauAndProd.size()!=2){
    std::cout<<" Warning!! Size of input vector != 4 !! "<<std::endl;
  }
  TLorentzVector LabFrame = TauAndProd.at(1);
  Setup(TauAndProd,LabFrame);
}

bool
SCalculator::isConfigured(){
  if(TauAndProd_LF.size()!=2){ std::cout<<"Error:   SCalculator is not Configured! Check  the size of input vector!  Size =  "<< TauAndProd_LF.size() <<std::endl; return false;} return true;
}

SCalculator::~SCalculator(){
}

TLorentzVector 
SCalculator::Boost(TLorentzVector pB, TLorentzVector frame){
   TMatrixT<double> transform(4,4);
   TMatrixT<double> result(4,1);
   TVectorT<double> vec(4); 
   TVector3 b;
   if(frame.Vect().Mag()==0){ std::cout<<"RH Boost is not set, perfrom calculation in the Lab Frame   "<<std::endl; return pB;}
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
TMatrixT<double> SCalculator::convertToMatrix(TVectorT<double> V){
  TMatrixT<double> M(V.GetNrows(),1);
  for(int i=0; i < M.GetNrows(); i++){
    M(i,0)=V(i);
  } return M;
}

TVector3
SCalculator::Rotate(TVector3 LVec, TVector3 Rot){
  TVector3 vec = LVec;
  vec.RotateZ(0.5*TMath::Pi() - Rot.Phi());  // not 0.5, to avoid warnings about 0 pT
  vec.RotateX(Rot.Theta());
  return vec;
}

TVector3
SCalculator::pvec(){
   TLorentzVector piLab = TauAndProd_LF.at(1);
   TLorentzVector tauLab = TauAndProd_LF.at(0);
   TVector3 Rot1 = tauLab.Vect();

   TLorentzVector tauLabR1    = tauLab;
   TLorentzVector piLabR1     = piLab ;
   tauLabR1.SetVect(Rotate(tauLabR1.Vect(),Rot1));
   piLabR1.SetVect(Rotate(piLabR1.Vect(),Rot1));

   TLorentzVector tauTau= Boost(tauLabR1,tauLabR1);
   TLorentzVector piTau= Boost(piLabR1,tauLabR1);

   return  piTau.Vect()*(1/piTau.Vect().Mag());
}
