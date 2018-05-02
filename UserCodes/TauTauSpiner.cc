#include "TauTauSpiner.h"
#include <iostream>
TauTauSpiner::TauTauSpiner():
  mrho(0.773),
  mpi(0.13957018),
  mtau(1.776),
  ma1(1.251),
  debug(false)
{ }  
TauTauSpiner::TauTauSpiner(vector<TLorentzVector> TauAndProd, string type, int TauCharge): 
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
  if(type=="a1"){
    if(TauAndProd.size()!=4){
      std::cout<<" Warning!! Size of 2nd  input vector  !=4 !! "<< " type:  "<<type<<std::endl;
    }
  }
  
  Configure(TauAndProd,  type,TauCharge);
 }
TauTauSpiner::TauTauSpiner(vector<TLorentzVector> TauAndProd1, string type1,vector<TLorentzVector> TauAndProd2, string type2,int TauCharge1,int TauCharge2): 
  mrho(0.773),
  mpi(0.13957018),
  mtau(1.776),
  ma1(1.251),
  debug(false)
 {
  if(type1=="lepton" || type1 =="pion"){
    if(TauAndProd1.size()!=2){ 
      std::cout<<" Warning!! Size of 1st  input vector  !=2  !! "<< " type:  "<<type1<<std::endl;
    }
  }
  if(type1=="rho"){
    if(TauAndProd1.size()!=3){
      std::cout<<" Warning!! Size of 1st  input vector  !=3  !! "<< " type:  "<<type1<<std::endl;
    }
  }
  if(type1=="a1"){
    if(TauAndProd1.size()!=4){
      std::cout<<" Warning!! Size of 2nd  input vector  !=4 !! "<< " type:  "<<type1<<std::endl;
    }
  }
  if(type2=="lepton" || type2 =="pion"){
    if(TauAndProd2.size()!=2){ 
      std::cout<<" Warning!! Size of 2nd  input vector  !=2  !! "<< " type:  "<<type2<<std::endl;
    }
  }
  if(type2=="rho"){
    if(TauAndProd2.size()!=3){
      std::cout<<" Warning!! Size of 2nd  input vector  !=3  !! "<< " type:  "<<type2<<std::endl;
    }
  }
  if(type2=="a1"){
    if(TauAndProd2.size()!=4){
      std::cout<<" Warning!! Size of 2nd  input vector  !=4 !! "<< " type:  "<<type2<<std::endl;
    }
  }
  
  ConfigurePair(TauAndProd1,type1,TauAndProd2,type2,TauCharge1,TauCharge2);
 }

void 
TauTauSpiner::Configure(vector<TLorentzVector> TauAndProd, string type,int TauCharge){
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
  if(type_=="a1"){
    TauA1OSPi  = TauAndProd.at(1);
    TauA1SSPi1 = TauAndProd.at(2);
    TauA1SSPi2 = TauAndProd.at(3);
    ProductLV = TauA1OSPi + TauA1SSPi1 + TauA1SSPi2;  
    taucharge_ = TauCharge;
  }
  InvisibleLV = TauLV - ProductLV;
  DPF_TauLV=  Boost(TauLV,ProductLV);
  DPF_InvisibleLV=  Boost(InvisibleLV,ProductLV);
  //  std::cout<<"Configure type 1 "<< type_; ProductLV.Print();

}

void 
TauTauSpiner::ConfigurePair(vector<TLorentzVector> TauAndProd1, string type1,  vector<TLorentzVector> TauAndProd2, string type2,int TauCharge1,int TauCharge2){
  TauLV1 = TauAndProd1.at(0);
  type1_=type1;
  if(type1_ == "lepton" || type1_=="pion"){
    ProductLV1= TauAndProd1.at(1);
  }
  if(type1_ == "rho" ){
    TauRhoPi1      = TauAndProd1.at(1);
    TauRhoPi01     = TauAndProd1.at(2);
    ProductLV1     = TauRhoPi1+TauRhoPi01;
    DPF_TauRhoPi1  = Boost(TauRhoPi1,ProductLV1);
    DPF_TauRhoPi01 = Boost(TauRhoPi01,ProductLV1);
  }

  if(type1_=="a1"){
    TauA1OSPi1  = TauAndProd1.at(1);
    TauA1SSPi11  = TauAndProd1.at(2);
    TauA1SSPi21  = TauAndProd1.at(3);
    ProductLV1 = TauA1OSPi1 + TauA1SSPi11 + TauA1SSPi21;  
    taucharge_=TauCharge1;
  }
  InvisibleLV1 = TauLV1 - ProductLV1;
  DPF_TauLV1=  Boost(TauLV1,ProductLV1);
  DPF_InvisibleLV1=  Boost(InvisibleLV1,ProductLV1);


  TauLV2 = TauAndProd2.at(0);
  type2_=type2;
  if(type2_ == "lepton" || type2_=="pion"){
    ProductLV2= TauAndProd2.at(1);
  }
  if(type2_ == "rho" ){
    TauRhoPi2      = TauAndProd2.at(1);
    TauRhoPi02     = TauAndProd2.at(2);
    ProductLV2     = TauRhoPi2+TauRhoPi02;
    DPF_TauRhoPi2  = Boost(TauRhoPi2,ProductLV2);
    DPF_TauRhoPi02 = Boost(TauRhoPi02,ProductLV2);
  }
  if(type2_=="a1"){
    TauA1OSPi2  = TauAndProd2.at(1);
    TauA1SSPi12  = TauAndProd2.at(2);
    TauA1SSPi22  = TauAndProd2.at(3);
    ProductLV2 = TauA1OSPi2 + TauA1SSPi12 + TauA1SSPi22;  
    taucharge_=TauCharge2;
  }
  InvisibleLV2 = TauLV2 - ProductLV2;
  DPF_TauLV2=  Boost(TauLV2,ProductLV2);
  DPF_InvisibleLV2=  Boost(InvisibleLV2,ProductLV2);

  TauLV2 = TauAndProd2.at(0);
  type2_=type2;

}


void 
TauTauSpiner::SetupLeg(string which){
    if(which=="first")
    {
      TauLV = TauLV1;
      type_=type1_;
      if(type_ == "lepton" || type_=="pion" ){ProductLV= ProductLV1;  }


      if(type_ == "rho" )
	{
	  TauRhoPi      = TauRhoPi1;
	  TauRhoPi0     = TauRhoPi01;
	  ProductLV     = TauRhoPi+TauRhoPi0;
	  DPF_TauRhoPi  = Boost(TauRhoPi,ProductLV);
	  DPF_TauRhoPi0 = Boost(TauRhoPi0,ProductLV);
	}
      if(type_=="a1")
	{
	  TauA1OSPi  =   TauA1OSPi1;
	  TauA1SSPi1  =  TauA1SSPi11;
	  TauA1SSPi2  =  TauA1SSPi21;
	  ProductLV = TauA1OSPi + TauA1SSPi1 + TauA1SSPi2;  
	}
      InvisibleLV = TauLV - ProductLV;
      DPF_TauLV=  Boost(TauLV,ProductLV);
      DPF_InvisibleLV=  Boost(InvisibleLV,ProductLV);
      //      std::cout<<"type 1 "<< type_ <<std::endl; ProductLV.Print();
    }
  if(which=="second")
    {
      TauLV = TauLV2;
      type_=type2_;
      if(type_ == "lepton" || type_=="pion" ){ProductLV= ProductLV2;  }
      if(type_ == "rho" )
	{
	  TauRhoPi      = TauRhoPi2;
	  TauRhoPi0     = TauRhoPi02;
	  ProductLV     = TauRhoPi+TauRhoPi0;
	  DPF_TauRhoPi  = Boost(TauRhoPi,ProductLV);
	  DPF_TauRhoPi0 = Boost(TauRhoPi0,ProductLV);
	}

      if(type_=="a1")
	{
	  TauA1OSPi  =   TauA1OSPi2;
	  TauA1SSPi1  =  TauA1SSPi12;
	  TauA1SSPi2  =  TauA1SSPi22;
	  ProductLV = TauA1OSPi + TauA1SSPi1 + TauA1SSPi2;  
	}

      InvisibleLV = TauLV - ProductLV;
      DPF_TauLV=  Boost(TauLV,ProductLV);
      DPF_InvisibleLV=  Boost(InvisibleLV,ProductLV);
      //      std::cout<<"type 2 "<< type_ <<std::endl; ProductLV.Print();

    }
}

bool  
TauTauSpiner::isConfigured(){if(TauLV.E()!=0 && ProductLV.E()!=0) return true; return false;}

bool  
TauTauSpiner::isPairConfigured(){if(TauLV1.E()!=0 && ProductLV1.E()!=0 && TauLV2.E()!=0 && ProductLV2.E()!=0) return true; return false;}


TauTauSpiner::~TauTauSpiner(){}



TLorentzVector 
TauTauSpiner::Boost(TLorentzVector pB, TLorentzVector frame)
   {
     TMatrixT<double> transform(4,4);
     TMatrixT<double> result(4,1);
     TVectorT<double> vec(4); 
     TVector3 b;
     if(frame.Vect().Mag()==0){ std::cout<<"TPI Boost is not set, perfrom calculation in the Lab Frame   "<<std::endl; return pB;}
     if(frame.E()==0){ std::cout<<" Caution: Please check that you perform boost correctly!  " <<std::endl; return pB;} 
     else   b=frame.Vect()*(1/frame.E());
     vec(0)  = pB.E();    vec(1)  = pB.Px();
     vec(2)  = pB.Py();   vec(3)  = pB.Pz();
     double gamma  = 1/sqrt( 1 - b.Mag2());
     transform(0,0)=gamma; transform(0,1) =- gamma*b.X() ;  transform(0,2) =  - gamma*b.Y();  transform(0,3) = - gamma*b.Z(); 
     transform(1,0)=-gamma*b.X(); transform(1,1) =(1+ (gamma-1)*b.X()*b.X()/b.Mag2()) ;  transform(1,2) = ((gamma-1)*b.X()*b.Y()/b.Mag2());  transform(1,3) = ((gamma-1)*b.X()*b.Z()/b.Mag2());
     transform(2,0)=-gamma*b.Y(); transform(2,1) = ((gamma-1)*b.Y()*b.X()/b.Mag2());  transform(2,2) = (1 + (gamma-1)*b.Y()*b.Y()/b.Mag2());  transform(2,3) =  ((gamma-1)*b.Y()*b.Z()/b.Mag2()); 
     transform(3,0)=-gamma*b.Z(); transform(3,1) =((gamma-1)*b.Z()*b.X()/b.Mag2()) ;  transform(3,2) = ((gamma-1)*b.Z()*b.Y()/b.Mag2());  transform(3,3) = (1 + (gamma-1)*b.Z()*b.Z()/b.Mag2()); 
     result=transform*convertToMatrix(vec);
     return TLorentzVector(result(1,0), result(2,0) ,result(3,0), result(0,0));
   }

double
TauTauSpiner::getOmega(string which)
   { 
     double omega=-999;
     SetupLeg(which);
     if(type_=="pion" || type_=="lepton")
       {
	 omega = 2*ProductLV.E()/TauLV.E() - 1;
       }
     if(type_=="rho")
       {
	 std::vector<TLorentzVector> particles; // tau, pi, pi0
	 particles.push_back(TauLV);
	 particles.push_back(TauRhoPi);
	 particles.push_back(TauRhoPi0);
	 rhoHelper rho;
	 rho.Configure(particles);
	 omega=rho.getOmegaRho();
       }
     if(type_=="a1")
       {
	 std::vector<TLorentzVector> particles;// tau, os,ss1,ss2
	 particles.push_back(TauLV);
	 particles.push_back(TauA1OSPi);
	 particles.push_back(TauA1SSPi1);
	 particles.push_back(TauA1SSPi2);
	 PolarimetricA1 a1;
	 a1.Configure(particles,taucharge_);
 	 omega = a1.getOmegaA1();
       }
     return omega;
   }

double
TauTauSpiner::getVisibleOmega(string which)
   {  
     double omega(-999.);
     SetupLeg(which);
     if(type_!="rho") {
       std::cout<<"This observable  is available for rho decay only (beta angle). a1 will be implemented later" << std::endl;
     }
     if(type_=="rho")
       {
	 std::vector<TLorentzVector> particles; // tau, pi, pi0
	 particles.push_back(TauLV);
	 particles.push_back(TauRhoPi);
	 particles.push_back(TauRhoPi0);
	 rhoHelper rho;
	 rho.Configure(particles);
	 omega=rho.getCosbetaRho();
       }
     return omega;
   }


double
TauTauSpiner::getOmegabar(string which){
  double omega=-999;
  SetupLeg(which);

  if(type_=="pion" || type_=="lepton"){
        omega = 2*ProductLV.E()/TauLV.E() - 1;
  }
  if(type_=="rho")
      {
	std::vector<TLorentzVector> particles; // tau, pi, pi0
	particles.push_back(TauLV);
	particles.push_back(TauRhoPi);
	particles.push_back(TauRhoPi0);
	rhoHelper rho;
	rho.Configure(particles);
	omega=rho.getOmegaRhoBar();
      }

     if(type_=="a1")
       {
	 std::vector<TLorentzVector> particles;// tau, os,ss1,ss2
	 particles.push_back(TauLV);
	 particles.push_back(TauA1OSPi);
	 particles.push_back(TauA1SSPi1);
	 particles.push_back(TauA1SSPi2);
	 PolarimetricA1 a1;
	 a1.Configure(particles,taucharge_);
	 omega = a1.getOmegaA1Bar();
       }
  return omega;
}


double
TauTauSpiner::getCombOmega(){
  double omega1 = getOmega("first");
  double omega2 = getOmega("second");
  double Omega=999.;
  Omega = (omega1 + omega2)/(1 + omega1*omega2);
  if(  std::isinf(std::fabs(Omega)) ||  std::isnan(std::fabs(Omega))) Omega  = -999.;
  return Omega;
}



double
TauTauSpiner::getCombOmegaBar(){
  double omega1 = getOmegabar("first");
  double omega2 = getOmegabar("second");
  double Omega=999.;
  Omega = (omega1 + omega2)/(1 + omega1*omega2);
  if(  std::isinf(std::fabs(Omega)) ||  std::isnan(std::fabs(Omega))) Omega  = -999.;
  return Omega;
}

double
TauTauSpiner::getCombVisibleOmega(){
  double omega1 = getVisibleOmega("first");
  double omega2 = getVisibleOmega("second");
  double Omega=999.;
  Omega = (omega1 + omega2)/(1 + omega1*omega2);
  if(  std::isinf(std::fabs(Omega)) ||  std::isnan(std::fabs(Omega))) Omega  = -999.;
  return Omega;
}


TLorentzVector
TauTauSpiner::getVisiblePairLV()
   { 
     if(!isPairConfigured()) return TLorentzVector(0,0,0,0);
     return ProductLV1+ProductLV2;
   }

TMatrixT<double> TauTauSpiner::convertToMatrix(TVectorT<double> V){
  TMatrixT<double> M(V.GetNrows(),1);
  for(int i=0; i < M.GetNrows(); i++){
    M(i,0)=V(i);
  } return M;
}
