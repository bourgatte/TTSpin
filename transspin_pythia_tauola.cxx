/**
 * Example of use of tauola C++ interfate. Pythia events are
 * generated with a stable tau. Taus are subseuently decay via
 * tauola, plots of polatization observables for tau-> mununu and tau-> pinu
 * are produced.
 *

 */ 
     
#include "Tauola/Log.h"
#include "Tauola/Tauola.h"
#include "Tauola/TauolaHepMCEvent.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h" 
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "UserCodes/a1Helper.h"
#include "UserCodes/TauDecaysHelper.h"
#include "UserCodes/TauPolInterface.h"
#include "UserCodes/PolarimetricA1.h"
#include "TLorentzVector.h"
#include "TComplex.h"
#include "TMatrixT.h"
#include "TVectorT.h"
#include "TMatrixTSym.h"
#include "TMath.h"

//pythia header files
#ifdef PYTHIA8180_OR_LATER
#include "Pythia8/Pythia.h" 
#include "Pythia8/Pythia8ToHepMC.h"
#else 
#include "Pythia.h"
#include "HepMCInterface.h"
#endif

//MC-TESTER header files
#include "Generate.h"
#include "HepMCEvent.H"
#include "Setup.H"
 
#include "tauola_print_parameters.h"
using namespace std;
using namespace Pythia8; 
using namespace Tauolapp;

int NumberOfEvents =20000; 
bool ApplyCut(false);
double pt_cut = 30;  // GeV - approximately correspond to CMS trigger
double eta_cut = 10; // - very large value, all events pass - at the moment switched off at all
int EventsToCheck=5;

// elementary test of HepMC typically executed before
// detector simulation based on http://home.fnal.gov/~mrenna/HCPSS/HCPSShepmc.html
// similar test was performed in Fortran
// we perform it before and after Tauola (for the first several events)



void checkMomentumConservationInEvent(HepMC::GenEvent *evt)
{
 	cout<<"List of stable particles: "<<endl;

	double px=0.0,py=0.0,pz=0.0,e=0.0;

	int barcodetau1vertex(0),barcodetau2vertex(0);
	for ( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin();  p != evt->particles_end(); ++p )
	  {
	    if( (*p)->status() == 1 )
	      {
		HepMC::FourVector m = (*p)->momentum();
		HepMC::GenVertex *productProductionvertex = (*p)->production_vertex();
		px+=m.px();
		py+=m.py();
		pz+=m.pz();
		e +=m.e();
	      }
	  }
	cout.precision(6);
	cout.setf(ios_base::floatfield);
}
TMatrixT<double> convertToMatrix(TVectorT<double> V){
  TMatrixT<double> M(V.GetNrows(),1);
  for(int i=0; i < M.GetNrows(); i++){
    M(i,0)=V(i);
  } return M;
}


TLorentzVector 
BoostR(TLorentzVector pB, TLorentzVector frame){
   TMatrixT<double> transform(4,4);
   TMatrixT<double> result(4,1);
   TVectorT<double> vec(4); 
   TVector3 b;
   if(frame.Vect().Mag()==0){ std::cout<<"RH Boost is not set, perfrom calculation in the Lab Frame   "<<std::endl; return pB;}
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

TVector3
Rotate(TVector3 LVec, TVector3 Rot){
  TVector3 vec = LVec;
  vec.RotateZ(0.5*TMath::Pi() - Rot.Phi());  // not 0.5, to avoid warnings about 0 pT
  vec.RotateX(Rot.Theta());
  return vec;
}


void redMinus(TauolaParticle *minus) // this is JAK1
{
  //   
  // this method can be used to redefine branching ratios in decay of tau-
  // either generally, or specific to  tau- with pointer *minus.
  //
  // Pointer *minus can be used to define properties of decays for taus
  // at specific point(s) in the event tree. Example: 
  // vector<TauolaParticle*> x=minus->getMothers();
  // and define special versions depending on x. 
  //
  // Any combination of methods
  // Tauola::setTauBr(int mode, double br);
  //Tauola::setTaukle(double bra1, double brk0, double brk0b,double brks);
   Tauola::setTaukle(1, 0,0,0);
  // can be called here 


   for(unsigned int dec=1; dec <23; dec++){
      double br =0.0; 
      if(dec ==2|| dec ==3|| dec ==4|| dec ==5) br=0.25;
      Tauola::setTauBr(dec, br);
   }

}

void redPlus(TauolaParticle *plus) // this is JAK2
{
  //   
  // this method can be used to redefine branching ratios in decay of tau+
  // either generally, or specific to  tau+ with pointer *plus.
  //
  // Pointer *plus can be used to define properties of decays for tau
  // at specific point(s) in the event tree. Example: 
  // vector<TauolaParticle*> x=plus->getMothers();
  // and define special versions depending on x. 
  //
  // Any combination of methods
  // Tauola::setTauBr(int mode, double br);
  //Tauola::setTaukle(double bra1, double brk0, double brk0b,double brks);
   Tauola::setTaukle(1, 0,0,0);
  // can be called here 
  for(unsigned int dec=1; dec <23; dec++){
     double br =0.0;
      if(dec ==2|| dec ==3|| dec ==4|| dec ==5) br=0.25;
     Tauola::setTauBr(dec, br);
   }
}

void SortPions(std::vector<HepMC::GenParticle > pionsvec)
{

  int npim(0),npip(0);
  int    OSMCPionIndex(0);
  int    SSMCPion1Index(0);
  int    SSMCPion2Index(0);

    HepMC::GenParticle os;
    HepMC::GenParticle ss1;
    HepMC::GenParticle ss2;
    for(int i=0; i<pionsvec.size(); i++){
      if( pionsvec.at(i).pdg_id()== 211) npip++;
      if( pionsvec.at(i).pdg_id()==-211) npim++;
    }
    if(npip == 1 && npim==2){
      int nss=0;
      for(int i=0; i<pionsvec.size(); i++){
	if(pionsvec.at(i).pdg_id()== 211){
	  OSMCPionIndex=i;
	}
	if(pionsvec.at(i).pdg_id()==-211 && nss ==0){
	  nss++;
	  SSMCPion1Index=i;
	}
	if(pionsvec.at(i).pdg_id()==-211 &&  nss == 1){
	  SSMCPion2Index=i;
	}
      }
    }
    if( npip== 2 && npim==1){
      int nss=0;
      for(int i=0; i<pionsvec.size(); i++){
	if(pionsvec.at(i).pdg_id()== -211){
	  
	  OSMCPionIndex=i;
	}
	if(pionsvec.at(i).pdg_id()==211 && nss ==0){
	  nss++;
	  SSMCPion1Index=i;
	}
	if(pionsvec.at(i).pdg_id()==211 &&  nss == 1){
	  SSMCPion2Index=i;
	}
      }
    }
    os=pionsvec.at(OSMCPionIndex);
    ss1=pionsvec.at(SSMCPion1Index);
    ss2=pionsvec.at(SSMCPion2Index);
    
  
    pionsvec.clear();
    pionsvec.push_back(os);
    pionsvec.push_back(ss1);
    pionsvec.push_back(ss2);
}

int main(int argc,char **argv){

  Log::SummaryAtExit();

  // Initialization of pythia
  Pythia pythia;
  Event& event = pythia.event;

  TString path= (TString)std::getenv("PWD") +"/output/";
  TString FileName  = path+"TauolaHelicity_" + TString(argv[1]) + ".root";
  TFile *file = new TFile(FileName,"RECREATE");


  //TFile *file = new TFile("HelicityVals_up.root","RECREATE");


  TH1F *rhobeta_plus= new TH1F("rhobeta_plus","#rho^{+}",50,-1.1,1.1);
  TH1F *rhobeta_minus= new TH1F("rhobeta_minus","#rho^{-}",50,-1.1,1.1);

  TH1F *pi_plus= new TH1F("pi_plus","#pi^{+}",50,-1.1,1.1);
  TH1F *pi_minus= new TH1F("pi_minus","#pi^{-} ",50,-1.1,1.1);

  TH1F *pi_pt= new TH1F("pi_pt","#pi^{-} ",100,0,100);
  TH1F *pi_ptcut= new TH1F("pi_ptcut","#pi^{-} ",100,0,100);
 
  //BT2  (-0.11994077    , 3.76959657E-03)   0.0000000 

  TH1F *pip_plus= new TH1F("pip_plus","#pi^{+}",50,-1.1,1.1);
  TH1F *pip_minus= new TH1F("pip_minus","#pi^{-} ",50,-1.1,1.1);


  TH1F *mu_plus= new TH1F("mu_plus","#mu^{+}",50,-1.1,1.1);
  TH1F *mu_minus= new TH1F("mu_minus","#mu^{-}",50,-1.1,1.1);

  TH1F *mu_polar_plus= new TH1F("mu_polar_plus","#mu^{+}",50,-1.1,1.1);
  TH1F *mu_polar_minus= new TH1F("mu_polar_minus","#mu^{-}",50,-1.1,1.1);

  
  TH1F *OmegaMuPi_plus= new TH1F("OmegaMuPi_plus","#omega_{#pi#mu}^{+}",50,-1.1,1.1);
  TH1F *OmegaMuPi_minus= new TH1F("OmegaMuPi_minus","#omega_{#pi#mu}^{-}",50,-1.1,1.1);
 
  TH1F *Omegapipi_plus= new TH1F("Omegapipi_plus","#omega_{#pi#pi}^{+}",50,-1.1,1.1);
  TH1F *Omegapipi_minus= new TH1F("Omegapipi_minus","#omega_{#pi#pi}^{-}",50,-1.1,1.1);
  
  TH1F *pipi_mass_plus= new TH1F("pipi_mass_plus","M_{#pi#pi}^{+}",60,0,120);
  TH1F *pipi_mass_minus= new TH1F("pipi_mass_minus","M_{#pi#pi}^{-}",60,0,120);
  
  TH1F *mupi_mass_plus= new TH1F("mupi_mass_plus","M_{#mu#pi}^{+}",60,0,120);
  TH1F *mupi_mass_minus= new TH1F("mupi_mass_minus","M_{#mu#pi}^{-}",60,0,120);
 
  TH1F *omega_murho_plus= new TH1F("omega_murho_plus","#omega_{#mu#rho}^{+}",50,-1.1,1.1);
  TH1F *omega_murho_minus= new TH1F("omega_murho_minus","#omega_{#mu#rho}^{-}",50,-1.1,1.1);


 TH1F *omega_rho_plus= new TH1F("omega_rho_plus","#omega_{#rho}^{+}",50,-1.1,1.1);
 TH1F *omega_rho_minus= new TH1F("omega_rho_minus","#omega_{#rho}^{-}",50,-1.1,1.1);

 TH1F *omegabar_rho_plus= new TH1F("omegabar_rho_plus","#bar{#omega}_{#rho}^{+}",50,-1.1,1.1);
 TH1F *omegabar_rho_minus= new TH1F("omegabar_rho_minus","#bar{#omega}_{#rho}^{-}",50,-1.1,1.1);


 TH1F *omega_a1_plus= new TH1F("omega_a1_plus","#omega_{a1}^{+}",50,-1.1,1.1);
 TH1F *omega_a1_minus= new TH1F("omega_a1_minus","#omega_{a1}^{-}",50,-1.1,1.1);

 // TH1F *omega_a1_plustest= new TH1F("omega_a1_plustest","#omega_{a1}^{+}",40,-1.1,1.1);
 // TH1F *omega_a1_minustest= new TH1F("omega_a1_minustest","#omega_{a1}^{+}",40,-1.1,1.1);

 TH1F *omega_a1p_plus= new TH1F("omega_a1p_plus","#omega_{a1}^{+}",50,-1.1,1.1);
 TH1F *omega_a1p_minus= new TH1F("omega_a1p_minus","#omega_{a1}^{-}",50,-1.1,1.1);



 TH1F *omegabar_a1_plus= new TH1F("omegabar_a1_plus","#bar{#omega}_{a1}^{+}",50,-1.1,1.1);
 TH1F *omegabar_a1_minus= new TH1F("omegabar_a1_minus","#bar{#omega}_{a1}^{-}",50,-1.1,1.1);
 
 TH1F *TRFomegabar_a1_plus= new TH1F("TRFomegabar_a1_plus","TRF #omega_{a1}^{+}",50,-1.1,1.1);
 TH1F *TRFomegabar_a1_minus= new TH1F("TRFomegabar_a1_minus","TRF #omega_{a1}^{-}",50,-1.1,1.1);

  TH1F *TRFomegabar_a1scalar_plus= new TH1F("TRFomegabar_a1scalar_plus","TRF scalar #omega a1",50,-1.1,1.1);
  TH1F *TRFomegabar_a1scalar_minus= new TH1F("TRFomegabar_a1scalar_minus","TRF scalar  #omega a1",50,-1.1,1.1);
 

 TH2F *cosbetacostheta_plus= new TH2F("cosbetacostheta_plus","cos#beta  cos#theta   {a1}^{+}",50,-1.1,1.1,50,-1.1,1.1);
 TH2F *cosbetacostheta_minus= new TH2F("cosbetacostheta_minus","cos#beta  cos#theta  {a1}^{-}",50,-1.1,1.1,50,-1.1,1.1);
 cosbetacostheta_plus->GetXaxis()->SetTitle("cos#beta"); cosbetacostheta_plus->GetYaxis()->SetTitle("cos#theta");
 cosbetacostheta_minus->GetXaxis()->SetTitle("cos#beta"); cosbetacostheta_minus->GetYaxis()->SetTitle("cos#theta");

 TH2F *cosbetacosthetarho_plus= new TH2F("cosbetacosthetarho_plus","cos#beta  cos#theta  #rho^{+}",50,-1.1,1.1,50,-1.1,1.1);
 TH2F *cosbetacosthetarho_minus= new TH2F("cosbetacosthetarho_minus","cos#beta  cos#theta  #rho^{-}" ,50,-1.1,1.1,50,-1.1,1.1);
 cosbetacosthetarho_plus->GetXaxis()->SetTitle("cos#beta"); cosbetacosthetarho_plus->GetYaxis()->SetTitle("cos#theta");
 cosbetacosthetarho_minus->GetXaxis()->SetTitle("cos#beta"); cosbetacosthetarho_minus->GetYaxis()->SetTitle("cos#theta");

 TH2F *TRFcosbetacostheta_plus= new TH2F("TRFcosbetacostheta_plus","cos#beta  cos#theta  {a1}^{+}",50,-1.1,1.1,50,-1.1,1.1);
 TH2F *TRFcosbetacostheta_minus= new TH2F("TRFcosbetacostheta_minus","cos#beta  cos#theta  {a1}^{-}",50,-1.1,1.1,50,-1.1,1.1);
 TRFcosbetacostheta_plus->GetXaxis()->SetTitle("cos#beta"); TRFcosbetacostheta_plus->GetYaxis()->SetTitle("cos#theta");
 TRFcosbetacostheta_minus->GetXaxis()->SetTitle("cos#beta"); TRFcosbetacostheta_minus->GetYaxis()->SetTitle("cos#theta");


 TH1F *omega_a1pi_plus= new TH1F("omega_a1pi_plus","#omega_{a1#pi}^{+}",50,-1.1,1.1);
 TH1F *omega_a1pi_minus= new TH1F("omega_a1pi_minus","#omega_{a1#pi}^{-}",50,-1.1,1.1);

 TH1F *omega_a1mu_plus= new TH1F("omega_a1mu_plus","#omega_{a1#mu}^{+}",50,-1.1,1.1);
 TH1F *omega_a1mu_minus= new TH1F("omega_a1mu_minus","#omega_{a1#mu}^{-}",50,-1.1,1.1);

 TH1F *omega_pirho_plus= new TH1F("omega_pirho_plus","#omega_{#pi#rho}^{+}",50,-1.1,1.1);
 TH1F *omega_pirho_minus= new TH1F("omega_pirho_minus","#omega_{#pi#rho}^{-}",50,-1.1,1.1);
  
 TH1F *omega_rhorho_plus= new TH1F("omega_rhorho_plus","#omega_{#rho#rho}^{+}",50,-1.1,1.1);
 TH1F *omega_rhorho_minus= new TH1F("omega_rhorho_minus","#omega_{#rho#rho}^{-}",50,-1.1,1.1);
  

 TH1F *mass_pirho_plus= new TH1F("mass_pirho_plus","M_{#pi#rho}^{+}",60,0,120);
 TH1F *mass_pirho_minus= new TH1F("mass_pirho_minus","M_{#pi#rho}^{-}",60,0,120);

 TH1F *mass_murho_plus= new TH1F("mass_murho_plus","M_{#mu#rho}^{+}",60,0,120);
 TH1F *mass_murho_minus= new TH1F("mass_murho_minus","M_{#mu#rho}^{-}",60,0,120);

 TH1F *mass_rhorho_plus= new TH1F("mass_rhorho_plus","M_{#rho#rho}^{+}",60,0,120);
 TH1F *mass_rhorho_minus= new TH1F("mass_rhorho_minus","M_{#rho#rho}^{-}",60,0,120);



 TH1F *a1_mcos2gamma_plus= new TH1F("a1_mcos2gamma_plus"," <cos2#gamma>^{+}",50,-0.5,0.5);
 TH1F *a1_mcos2gamma_minus= new TH1F("a1_mcos2gamma_minus"," <cos2#gamma>^{-}",50,-0.5,0.5);



 // TH1F *omega_a1rho_plus= new TH1F("omega_a1rho_plus","#omega_{a1#rho}^{+}",50,-1.1,1.1);
 // TH1F *omega_a1rho_minus= new TH1F("omega_a1rho_minus","#omega_{a1#rho}^{-}",50,-1.1,1.1);


 // TH1F *omega_pirho_plus= new TH1F("omega_pirho_plus","#omega_{#pi#rho}^{+}",50,-1.1,1.1);
 // TH1F *omega_pirho_minus= new TH1F("omega_pirho_minus","#omega_{#pi#rho}^{-}",50,-1.1,1.1);



TH1F *s1= new TH1F("s1","s2",40,0,2);
TH1F *s2= new TH1F("s2","s2",40,0,2);
TH1F *qq= new TH1F("qq","qq",40,0,3);
TH1F *hmag= new TH1F("hmag","hmag",40,0.5,1.5);


 TH2F *s1s2= new TH2F("s1s2","s1s2",20,0,2,20,0,2);




  // Pythia8 HepMC interface depends on Pythia8 version
#ifdef PYTHIA8180_OR_LATER
  HepMC::Pythia8ToHepMC ToHepMC;
#else
  HepMC::I_Pythia8 ToHepMC;
#endif

  //pythia.readString("HadronLevel:all = off");
  pythia.readString("HadronLevel:Hadronize = off");
  pythia.readString("SpaceShower:QEDshowerByL = off");
  pythia.readString("SpaceShower:QEDshowerByQ = off");
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = off");

  // Tauola is currently set to undecay taus. Otherwise, uncomment this line.
  // Uncommenting it will speed up the test significantly
  //  pythia.particleData.readString("15:mayDecay = off");

  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  pythia.readString("23:onMode = off"); 
  pythia.readString("23:onIfAny = 15");
  //  pythia.readString("23:onIfMatch = 15 -15");

  pythia.init( 11, -11, 92.);          //electron positron collisions

  // Set up Tauola

  // Set Tauola decay mode (if needed)
  //Tauola::setSameParticleDecayMode(2);     //19 and 22 contains K0 
    //   Tauola::setOppositeParticleDecayMode(3); // 20 contains eta

  // Set Higgs scalar-pseudoscalar mixing angle
  //  Tauola::setHiggsScalarPseudoscalarMixingAngle(0.7853);
  //  Tauola::setHiggsScalarPseudoscalarPDG(25);

  Tauola::initialize();
  //  const char* str="hhu1";
  std::cout<<"  "<< TMath::Hash(argv[1])<<std::endl;
  std::cout<<" ---------------  "<< time(NULL)<<std::endl;
 

  //  Tauola::setSeed(time(NULL), 0, 0);
  Tauola::setSeed(TMath::Hash(argv[1]), 0, 0);
  tauola_print_parameters(); // Prints TAUOLA  parameters (residing inside its library): e.g. to test user interface

  // Our default units are GEV and MM, that will be outcome  units after TAUOLA
  // if HepMC unit variables  are correctly set. 
  // with the following coice you can fix the units for final outcome:
  //  Tauola::setUnits(Tauola::GEV,Tauola::MM); 
  //  Tauola::setUnits(Tauola::MEV,Tauola::CM); 

  // Other usefull settings:
  //  Tauola::setEtaK0sPi(0,0,0);  // switches to decay eta K0_S and pi0 1/0 on/off. 
  //  Tauola::setTauLifetime(0.0); //new tau lifetime in mm
    Tauola::spin_correlation.setAll(true);

    Log::LogDebug(true);

      Tauola::setRedefineTauMinus(redMinus);  // activates execution of routine redMinus in TAUOLA interface
      Tauola::setRedefineTauPlus(redPlus);    // activates execution of routine redPlus  in TAUOLA interface
      // Tauola::setRedefineTauMinus(redPlus);  // activates execution of routine redMinus in TAUOLA interface
      // Tauola::setRedefineTauPlus(redMinus);    // activates execution of routine redPlus  in TAUOLA interface



  MC_Initialize();

  // Begin event loop. Generate event.
  for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent){

    if(iEvent%1000==0) Log::Info()<<"Event: "<<iEvent<<endl;
    if (!pythia.next()) continue;
    //    std::cout<<"-------------- "<<std::endl;
    // Convert event record to HepMC
    HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();

    // Conversion needed if HepMC use different momentum units
    // than Pythia. However, requires HepMC 2.04 or higher.
    HepMCEvt->use_units(HepMC::Units::GEV,HepMC::Units::MM);

    ToHepMC.fill_next_event(event, HepMCEvt);

    if(iEvent<EventsToCheck)
    {
      cout<<"                                          "<<endl;
      cout<<"Momentum conservation chceck BEFORE/AFTER Tauola"<<endl;
      checkMomentumConservationInEvent(HepMCEvt);
    }

    // Run TAUOLA on the event
    TauolaHepMCEvent * t_event = new TauolaHepMCEvent(HepMCEvt);

    // Since we let Pythia decay taus, we have to undecay them first.
    t_event->undecayTaus();
    t_event->decayTaus();
    delete t_event; 

    if(iEvent<EventsToCheck)
    {
      checkMomentumConservationInEvent(HepMCEvt);
    }

    int JAK1(0); int SubJAK1(0);
    HepMC::GenParticle *FirstTau;
    std::vector<HepMC::GenParticle > FirstTauProducts;
    int JAK2(0); int SubJAK2(0);
    HepMC::GenParticle *SecondTau;
    std::vector<HepMC::GenParticle > SecondTauProducts;
    std::vector<HepMC::GenParticle > A1Pions;
    std::vector<HepMC::GenParticle > A1Pions1;
    std::vector<HepMC::GenParticle > A1Pions2;
    std::vector<HepMC::GenParticle > SortA1Pions;  //os, ss1, ss2
    for ( HepMC::GenEvent::particle_const_iterator p =HepMCEvt->particles_begin();  p != HepMCEvt->particles_end(); ++p ){  
      if((*p)->pdg_id()==15){
	FirstTau = *p;
	for ( HepMC::GenEvent::particle_const_iterator d =HepMCEvt->particles_begin();  d != HepMCEvt->particles_end(); ++d ){  
	  if((*d)->pdg_id()!=15){
	    if((*p)->end_vertex() == (*d)->production_vertex()){
	      FirstTauProducts.push_back(**d);
	      if(abs((*d)->pdg_id()) ==  12) {JAK1 =1; 
	      }else if(abs((*d)->pdg_id()) ==  14){ JAK1=2;
	      }else if(abs((*d)->pdg_id())==  211){ JAK1 = 3;/*std::cout<<"pion of negative tau "<< (*d)->pdg_id() <<std::endl;*/}
	      

	      if( abs((*d)->pdg_id())==213 ){
		JAK1 = 4;
		for ( HepMC::GenEvent::particle_const_iterator dd =HepMCEvt->particles_begin();  dd != HepMCEvt->particles_end(); ++dd ){  
		  if( abs((*dd)->pdg_id())!=213  ){
		    if((*d)->end_vertex() == (*dd)->production_vertex()){

		      FirstTauProducts.push_back(**dd);
		    }
		  }
		}
	      }
		if( abs((*d)->pdg_id())==20213 ){
		  JAK1 = 5; int npi(0);
		  for ( HepMC::GenEvent::particle_const_iterator dd =HepMCEvt->particles_begin();  dd != HepMCEvt->particles_end(); ++dd ){  
		    if( abs((*dd)->pdg_id())!=20213  ){
		      if((*d)->end_vertex() == (*dd)->production_vertex()){
		       
			FirstTauProducts.push_back(**dd);
			if(abs((*dd)->pdg_id())==  211)   {
			  A1Pions1.push_back(**dd); npi++;

			}
		      }
		    }
		  }
		   // std::cout<<"a1 pions  of negative tau " << A1Pions.size()<<std::endl;
		   // std::cout<<      A1Pions1.at(2).pdg_id()<< "  " <<A1Pions1.at(2).momentum().px()<<std::endl;
		   // std::cout<<      A1Pions1.at(0).pdg_id()<< "  " <<A1Pions1.at(0).momentum().px()<<std::endl;
		   // std::cout<<      A1Pions1.at(1).pdg_id()<< "  " <<A1Pions1.at(1).momentum().px()<<std::endl;
		   

		  if(npi==3) SubJAK1=51; else SubJAK1=52;
		}
	    }
	  } 
	}
      }
    
      
      if((*p)->pdg_id()==-15){
	SecondTau = *p;
	for ( HepMC::GenEvent::particle_const_iterator d =HepMCEvt->particles_begin();  d != HepMCEvt->particles_end(); ++d )
	  {  
	    if((*d)->pdg_id()!=-15){
	     
	      if((*p)->end_vertex() == (*d)->production_vertex()){
		SecondTauProducts.push_back(**d);
		if(abs((*d)->pdg_id()) ==  12) {JAK2 =1; 
		}else if(abs((*d)->pdg_id()) ==  14){ JAK2=2;
		}else if(abs((*d)->pdg_id())==  211){ JAK2=3;/*	std::cout<<"pion of positive tau "<< (*d)->pdg_id() <<std::endl;*/}
	
		
 
		if( abs((*d)->pdg_id())==213 ){
		  JAK2 = 4;
		  for ( HepMC::GenEvent::particle_const_iterator dd =HepMCEvt->particles_begin();  dd != HepMCEvt->particles_end(); ++dd ){  
		    if( abs((*dd)->pdg_id())!=213  ){
		      if((*d)->end_vertex() == (*dd)->production_vertex()){
			SecondTauProducts.push_back(**dd);
		      }
		    }
		  }
		}

		if( abs((*d)->pdg_id())==20213 ){
		  JAK2 = 5; int npi(0);
		  for ( HepMC::GenEvent::particle_const_iterator dd =HepMCEvt->particles_begin();  dd != HepMCEvt->particles_end(); ++dd ){  
		    if( abs((*dd)->pdg_id())!=20213  ){
		      if((*d)->end_vertex() == (*dd)->production_vertex()){
			SecondTauProducts.push_back(**dd);
			if(abs((*dd)->pdg_id())==  211){
			  A1Pions2.push_back(**dd); npi++;
			}
		      }
		    }
		  }

		  if(npi==3) SubJAK2=51; else SubJAK2=52;
		}
		
	      }
	    } 
	  }
      }
    }

    bool HelPlus=false;
    bool HelMinus=false;
    if(Tauola::getHelPlus() == 1 )HelMinus=true;
    if(Tauola::getHelPlus() ==-1)HelPlus=true;

    bool HelPlus2=false;
    bool HelMinus2=false;
    if(Tauola::getHelMinus() == 1 )HelMinus2=true;
    if(Tauola::getHelMinus() ==-1)HelPlus2=true;


    int HelWeightPlus = HelPlus;
    int HelWeightMinus = HelMinus;
    // int HelWeightPlus = 0;
    // int HelWeightMinus = 0;


    // std::cout<<"  two helicities  "<< HelMinus <<"  "<<HelPlus <<std::endl;
    // std::cout<<"  s  "<< HelMinus2 <<"  "<<HelPlus2 <<std::endl;
    int tauHelicity  = Tauola::getHelPlus();

   
    TLorentzVector tau1(0,0,0,0);
    TLorentzVector mu1(0,0,0,0);
    TLorentzVector numu1(0,0,0,0);
    TLorentzVector nutau1(0,0,0,0);
    TLorentzVector nutau2(0,0,0,0);
    TLorentzVector pi2(0,0,0,0);
    TLorentzVector pi1(0,0,0,0);
    TLorentzVector tau2(0,0,0,0);

    TLorentzVector rhopi2(0,0,0,0);
    TLorentzVector rhopi02(0,0,0,0);

    TLorentzVector a1ospi(0,0,0,0);
    TLorentzVector a1ss1pi(0,0,0,0);
    TLorentzVector a1ss2pi(0,0,0,0);
    TLorentzVector a1(0,0,0,0);




    TauPolInterface Mu1;
    TauPolInterface Mu2;
    TauPolInterface Pi1;
    TauPolInterface Pi2;
    TauPolInterface Rho2;
    TauPolInterface Rho1;
    TauPolInterface A2;
    TauPolInterface A1;

    int taucharge1;
    int taucharge2;

    TauPolInterface TauPolPiPi;
    TauPolInterface TauPolMuPi;
    TauPolInterface TauPolMuRho;
    TauPolInterface TauPolPiRho;
    TauPolInterface TauPolPiA1;
    TauPolInterface TauPolMuA1;
    TauPolInterface TauPolRhoRho;



    a1Helper a1h;
    a1Helper a1hh;
    rhoHelper RhoHelp;
    PolarimetricA1 Polarimetr;
    PolarimetricA1 Polarimetr1;
    TauPolInterface TauPolPi1;
    TauPolInterface TauPolPi2;
   
    TauPolInterface TauPolRho1;

    TauPolInterface TauPolMu1;
    TauPolInterface TauPolRho2;
    TauPolInterface TauPolA1;
 
 
    

    vector<TLorentzVector> tauandprod1,tauandprod2, tauandprodMuon2;
    vector<TLorentzVector> tauandprodMuon1,tauandprodRho,tauandprodRho2,tauandprodA1;
    tau1.SetPxPyPzE(FirstTau->momentum().px(), FirstTau->momentum().py(), FirstTau->momentum().pz(), FirstTau->momentum().e());
    tau2.SetPxPyPzE(SecondTau->momentum().px(), SecondTau->momentum().py(), SecondTau->momentum().pz(), SecondTau->momentum().e());

    bool passed1(true);
    bool passed2(true);
    if(ApplyCut){
      passed1 = false;
      passed2 = false;
      if(JAK1==2){

	for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
	  if(abs(a->pdg_id())==13){
	    TLorentzVector  prod(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  );
	    if(prod.Pt() > pt_cut ) passed1=true;
	  }
	}
      }
      if(JAK1==3){
	for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
	  if(abs(a->pdg_id())==211){
	    TLorentzVector  prod(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  );
	    if(prod.Pt() > pt_cut ) passed1=true;
	  }
	}
      }

      if(JAK1==4){
	TLorentzVector prod(0,0,0,0);
	for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
	  if(abs(a->pdg_id())==211){prod += TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  );}
	  if(abs(a->pdg_id())==111){prod += TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  );}
	}
	if(prod.Pt() > pt_cut ) passed1=true;
      }
      
      if(JAK1==5){
	TLorentzVector prod(0,0,0,0);

	prod+=TLorentzVector(A1Pions1.at(0).momentum().px(), A1Pions1.at(0).momentum().py(), A1Pions1.at(0).momentum().pz(), A1Pions1.at(0).momentum().e());
	prod+=TLorentzVector(A1Pions1.at(1).momentum().px(), A1Pions1.at(1).momentum().py(), A1Pions1.at(1).momentum().pz(), A1Pions1.at(1).momentum().e());
	prod+=TLorentzVector(A1Pions1.at(2).momentum().px(), A1Pions1.at(2).momentum().py(), A1Pions1.at(2).momentum().pz(), A1Pions1.at(2).momentum().e());
	
	if(prod.Pt() > pt_cut ) passed1=true;
      }

      if(JAK2==2){
	for(std::vector<HepMC::GenParticle>::const_iterator a = SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
	  if(abs(a->pdg_id())==13){
	    TLorentzVector  prod(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  );
	    if(prod.Pt() > pt_cut) passed2=true;
	  }
	}
      }

      if(JAK2==3){
	for(std::vector<HepMC::GenParticle>::const_iterator a = SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
	  if(abs(a->pdg_id())==211){
	    TLorentzVector  prod(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  );
	    if(prod.Pt() > pt_cut) passed2=true;
	  }
	}
      }

      if(JAK2==4){
	TLorentzVector prod(0,0,0,0);
	for(std::vector<HepMC::GenParticle>::const_iterator a = SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
	  if(abs(a->pdg_id())==211){prod += TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  );}
	  if(abs(a->pdg_id())==111){prod += TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  );}
	}
	if(prod.Pt() > pt_cut ) passed2=true;
      }

      if(JAK2==5){
	TLorentzVector prod(0,0,0,0);
	prod+=TLorentzVector(A1Pions2.at(0).momentum().px(), A1Pions2.at(0).momentum().py(), A1Pions2.at(0).momentum().pz(), A1Pions2.at(0).momentum().e());
	prod+=TLorentzVector(A1Pions2.at(1).momentum().px(), A1Pions2.at(1).momentum().py(), A1Pions2.at(1).momentum().pz(), A1Pions2.at(1).momentum().e());
	prod+=TLorentzVector(A1Pions2.at(2).momentum().px(), A1Pions2.at(2).momentum().py(), A1Pions2.at(2).momentum().pz(), A1Pions2.at(2).momentum().e());
	if(prod.Pt() > pt_cut ) passed2=true;
      }
    }

    //---------------------------------------------------------------------------
    //--------------------  tau-
    if(JAK1==2){
      vector<TLorentzVector> tauandprod;
      vector<TLorentzVector> tauandprod_neutrinos; // first goes tau neutrino
      tauandprod.push_back(TLorentzVector(FirstTau->momentum().px(), FirstTau->momentum().py(), FirstTau->momentum().pz(), FirstTau->momentum().e()));
      for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
	
	if(abs(a->pdg_id())==16){tauandprod_neutrinos.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );} 
	if(abs(a->pdg_id())==14){tauandprod_neutrinos.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );} 
	if(abs(a->pdg_id())==13){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );} }
      
      tauandprod1=tauandprod;
      Mu1.Configure(tauandprod1,"lepton");
      
      mu_plus->Fill(Mu1.getOmega(),HelWeightPlus);
      mu_minus->Fill(Mu1.getOmega(),HelWeightMinus);


      //-------------- muon polarimetric vector 


      TLorentzVector PLab = tauandprod.at(0);
      TLorentzVector q1Lab= tauandprod.at(1);
      TLorentzVector q2Lab = tauandprod_neutrinos.at(1);
      TLorentzVector NLab = tauandprod_neutrinos.at(0);

     
      TVector3 Rot1 = PLab.Vect();

      TLorentzVector PLabR1    = PLab;
      TLorentzVector q1LabR1     = q1Lab;
      TLorentzVector q2LabR1     = q2Lab;
      TLorentzVector NLabR1     = NLab;


      PLabR1.SetVect(Rotate(PLabR1.Vect(),Rot1));
      q1LabR1.SetVect(Rotate(q1LabR1.Vect(),Rot1));
      q2LabR1.SetVect(Rotate(q2LabR1.Vect(),Rot1));
      NLabR1.SetVect(Rotate(NLabR1.Vect(),Rot1));
      //
      // PLabR1.Print();
      TLorentzVector PTau = BoostR(PLabR1,PLabR1);
      TLorentzVector q1Tau= BoostR(q1LabR1,PLabR1);
      TLorentzVector q2Tau =BoostR(q2LabR1,PLabR1);
      TLorentzVector NTau = BoostR(NLabR1,PLabR1);



      // PTau.Print();
      // q2Tau.Print();
      //std::cout<<" Pq2  "<< q2Tau*PTau << std::endl;
      TVector3 TauLabDir = PLabR1.Vect()*(1/PLabR1.Vect().Mag());
      //denom = 
      TVector3 MuPolar   =-PTau.M()* q2Tau.Vect()* (1/(q2Tau*PTau));
      // MuPolar.Print();
      // TauLabDir.Print();
      // double denom = (tauTau*q1Tau)*(q2Tau*muTau);
      // TVector3 muTauDir =-tauTau.M()* q1Tau.Vect()*(muTau*q2Tau) * (1/denom);




      mu_polar_plus->Fill(MuPolar*TauLabDir,HelWeightPlus);
      mu_polar_minus->Fill(MuPolar*TauLabDir,HelWeightMinus);


    
      //-------------- muon polarimetric vector 


    }
    
    if(JAK1==3){
      vector<TLorentzVector> tauandprod;
      tauandprod.push_back(TLorentzVector(FirstTau->momentum().px(), FirstTau->momentum().py(), FirstTau->momentum().pz(), FirstTau->momentum().e()));
       for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
	 if(abs(a->pdg_id())==211){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );}}
       tauandprod1 = tauandprod;
       Pi1.Configure(tauandprod1,"pion");
       if(ApplyCut){
	 if(passed1){
	   pip_plus->Fill(Pi1.getOmega(),HelWeightPlus);
	   pip_minus->Fill(Pi1.getOmega(),HelWeightMinus);
	   pi_ptcut->Fill((tauandprod.at(1) ).Pt());
	 }
       }else{
	 pip_plus->Fill(Pi1.getOmega(),HelWeightPlus);
	 pip_minus->Fill(Pi1.getOmega(),HelWeightMinus);
	 pi_pt->Fill((tauandprod.at(1) ) .Pt());
       }


    }
    
  
    

      if(JAK1==4){
       vector<TLorentzVector> tauandprod;
       tauandprod.push_back(TLorentzVector(FirstTau->momentum().px(), FirstTau->momentum().py(), FirstTau->momentum().pz(), FirstTau->momentum().e()));
       for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
     	if(abs(a->pdg_id())==211){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );}
     	if(abs(a->pdg_id())==111){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );}
       }
       tauandprod1=tauandprod;  
       Rho1.Configure(tauandprod1,"rho");
       TauPolRho1.Configure(tauandprod1,"rho");
       RhoHelp.Configure(tauandprod1);
    	      


       if(ApplyCut){
	 if(passed1){
	   rhobeta_plus->Fill(Rho1.getVisibleOmega(),HelWeightPlus);
	   rhobeta_minus->Fill(Rho1.getVisibleOmega(),HelWeightMinus);
	   
	   omega_rho_plus->Fill(Rho1.getOmega(),HelWeightPlus);
	   omega_rho_minus->Fill(Rho1.getOmega(),HelWeightMinus);
	   
	   omegabar_rho_plus->Fill(Rho1.getOmegabar(),HelWeightPlus);
	   omegabar_rho_minus->Fill(Rho1.getOmegabar(),HelWeightMinus);
	   
	   cosbetacosthetarho_plus->Fill(RhoHelp.getCosbetaRho(),RhoHelp.getCosthetaRho(),HelWeightPlus);
	   cosbetacosthetarho_minus->Fill(RhoHelp.getCosbetaRho(),RhoHelp.getCosthetaRho(),HelWeightMinus);
	 }
       }else{
	 rhobeta_plus->Fill(Rho1.getVisibleOmega(),HelWeightPlus);
	 rhobeta_minus->Fill(Rho1.getVisibleOmega(),HelWeightMinus);
	 
	 omega_rho_plus->Fill(Rho1.getOmega(),HelWeightPlus);
	 omega_rho_minus->Fill(Rho1.getOmega(),HelWeightMinus);
	 
	 omegabar_rho_plus->Fill(Rho1.getOmegabar(),HelWeightPlus);
	 omegabar_rho_minus->Fill(Rho1.getOmegabar(),HelWeightMinus);
	 
	 cosbetacosthetarho_plus->Fill(RhoHelp.getCosbetaRho(),RhoHelp.getCosthetaRho(),HelWeightPlus);
	 cosbetacosthetarho_minus->Fill(RhoHelp.getCosbetaRho(),RhoHelp.getCosthetaRho(),HelWeightMinus);
       }

     }

     if(JAK1==5 && SubJAK1==51){
       vector<TLorentzVector> particles;
       particles.clear();
       SortPions(A1Pions1);
       int taucharge =  (A1Pions1.at(0).pdg_id()+A1Pions1.at(1).pdg_id()+A1Pions1.at(2).pdg_id() > 0) ? 1 : -1;
       taucharge1=taucharge;

       a1ss1pi.SetPxPyPzE(A1Pions1.at(0).momentum().px(), A1Pions1.at(0).momentum().py(), A1Pions1.at(0).momentum().pz(), A1Pions1.at(0).momentum().e());
       a1ss2pi.SetPxPyPzE(A1Pions1.at(1).momentum().px(), A1Pions1.at(1).momentum().py(), A1Pions1.at(1).momentum().pz(), A1Pions1.at(1).momentum().e());
       a1ospi.SetPxPyPzE(A1Pions1.at(2).momentum().px(), A1Pions1.at(2).momentum().py(), A1Pions1.at(2).momentum().pz(), A1Pions1.at(2).momentum().e());
       particles.push_back(tau1);
       particles.push_back(a1ospi);
       particles.push_back(a1ss1pi);
       particles.push_back(a1ss2pi);
 


       Polarimetr1.Configure(particles, taucharge);
       A1.Configure(particles,"a1",taucharge);

       omega_a1p_minus->Fill(Polarimetr1.getOmegaA1Bar(),HelWeightMinus);  omega_a1p_plus->Fill(Polarimetr1.getOmegaA1Bar(),HelWeightPlus);
       //  std::cout<<"omegabar   "<<Polarimetr1.getOmegaA1Bar() <<std::endl;
       a1hh.Configure(particles, a1ospi+a1ss1pi+a1ss2pi);
       TauPolA1.Configure(particles,"a1");
       tauandprod1=particles;
       //      omega_a1_minus->Fill(a1h.getA1omega(),HelWeightMinus);                                                          omegabar_a1_plus->Fill(a1h.getA1omega(),HelWeightPlus);

       omegabar_a1_minus->Fill(a1hh.getA1omega(),HelWeightMinus);                                                          omegabar_a1_plus->Fill(a1hh.getA1omega(),HelWeightPlus);

     }





     //--------------------  tau+
     if(JAK2==2){
       vector<TLorentzVector> tauandprodmu2;
       tauandprodmu2.push_back(TLorentzVector(SecondTau->momentum().px(), SecondTau->momentum().py(), SecondTau->momentum().pz(), SecondTau->momentum().e()));
       for(std::vector<HepMC::GenParticle>::const_iterator a = SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
     	if(abs(a->pdg_id())==13){tauandprodmu2.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );} }
       tauandprod2=tauandprodmu2;
       Mu2.Configure(tauandprod2,"lepton");
     }
     
     if(JAK2==3){
       vector<TLorentzVector> tauandprod;
       tauandprod.push_back(TLorentzVector(SecondTau->momentum().px(), SecondTau->momentum().py(), SecondTau->momentum().pz(), SecondTau->momentum().e()));
       for(std::vector<HepMC::GenParticle>::const_iterator a = SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
	 if(abs(a->pdg_id())==211){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );} }
       tauandprod2 = tauandprod;
       Pi2.Configure(tauandprod2,"pion");
       if(ApplyCut){
	 if(passed1){
	   pi_plus->Fill(Pi2.getOmega(),HelWeightPlus);
	   pi_minus->Fill(Pi2.getOmega(),HelWeightMinus);
	 }
       }else{
	 pi_plus->Fill(Pi2.getOmega(),HelWeightPlus);
	 pi_minus->Fill(Pi2.getOmega(),HelWeightMinus);
       }
     }
     
     if(JAK2==4){
       vector<TLorentzVector> tauandprod;
       tauandprod.push_back(TLorentzVector(SecondTau->momentum().px(), SecondTau->momentum().py(), SecondTau->momentum().pz(), SecondTau->momentum().e()));
       //      std::cout<<" ---- "<<std::endl;
       for(std::vector<HepMC::GenParticle>::const_iterator a =SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
	 //	std::cout<<" rho decays:   "<< a->pdg_id() << std::endl;
	 if(abs(a->pdg_id())==211){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );}
	 if(abs(a->pdg_id())==111){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );}
       }
       tauandprod2=tauandprod;  
       Rho2.Configure(tauandprod2,"rho");


     }

    
    if(JAK2==5 && SubJAK2==51){
      vector<TLorentzVector> particles;
      particles.clear();
      SortPions(A1Pions2);
      
      int taucharge =  (A1Pions2.at(0).pdg_id()+A1Pions2.at(1).pdg_id()+A1Pions2.at(2).pdg_id() > 0) ? 1 : -1;
      taucharge2=taucharge;
      a1ss1pi.SetPxPyPzE(A1Pions2.at(0).momentum().px(), A1Pions2.at(0).momentum().py(), A1Pions2.at(0).momentum().pz(), A1Pions2.at(0).momentum().e());
      a1ss2pi.SetPxPyPzE(A1Pions2.at(1).momentum().px(), A1Pions2.at(1).momentum().py(), A1Pions2.at(1).momentum().pz(), A1Pions2.at(1).momentum().e());
      a1ospi.SetPxPyPzE(A1Pions2.at(2).momentum().px(), A1Pions2.at(2).momentum().py(), A1Pions2.at(2).momentum().pz(), A1Pions2.at(2).momentum().e());
      
      particles.push_back(tau2);
      particles.push_back(a1ospi);
      particles.push_back(a1ss1pi);
      particles.push_back(a1ss2pi);
      a1h.Configure(particles, a1ospi+a1ss1pi+a1ss2pi);
      TauPolA1.Configure(particles,"a1");
      tauandprod2=particles;
      Polarimetr.Configure(particles, taucharge);
      A2.Configure(particles,"a1", taucharge);
      
      if(ApplyCut){
	if(passed2){
	  
	  TRFomegabar_a1_minus->Fill(a1h.TRF_vgetA1omega(),HelWeightMinus);                                      TRFomegabar_a1_plus->Fill(a1h.TRF_vgetA1omega(),HelWeightPlus);
	  TRFomegabar_a1scalar_minus->Fill(a1h.TRF_vgetA1omega("scalar"),HelWeightMinus);                  TRFomegabar_a1scalar_plus->Fill(a1h.TRF_vgetA1omega("scalar"),HelWeightPlus);
	  cosbetacostheta_minus->Fill(a1h.cosbeta(),a1h.costhetaLF(),HelWeightMinus);                              cosbetacostheta_plus->Fill(a1h.cosbeta(),a1h.costhetaLF(),HelWeightPlus);
	  TRFcosbetacostheta_minus->Fill(a1h.TRF_cosbeta(),a1h.costhetaLF(),HelWeightMinus);                  TRFcosbetacostheta_plus->Fill(a1h.TRF_cosbeta(),a1h.costhetaLF(),HelWeightPlus);
	  
	  s1->Fill((a1ospi +a1ss1pi).M() );
	  s2->Fill((a1ospi +a1ss2pi).M() );
	  qq->Fill((a1ospi+a1ss1pi+a1ss2pi).M());

	  omega_a1_minus->Fill(Polarimetr.getOmegaA1Bar(),HelWeightMinus);  omega_a1_plus->Fill(Polarimetr.getOmegaA1Bar(),HelWeightPlus);
	  hmag->Fill(Polarimetr.PVC().Vect().Mag());
	  s1s2->Fill((a1ospi+a1ss2pi).M2(),(a1ospi+a1ss1pi).M2());
	  int hel(0);
	  if(HelWeightMinus==1)     hel = -1;
	  if(HelWeightPlus==1)      hel =  1;
	  if(HelWeightPlus==1)      a1_mcos2gamma_plus->Fill(a1h.getMoment(a1h.costhetaLF(),"c2g", 1));                 
	  if(HelWeightMinus==1)     a1_mcos2gamma_minus->Fill(a1h.getMoment(a1h.costhetaLF(),"c2g", -1));                 
	}
      }else{
	TRFomegabar_a1_minus->Fill(a1h.TRF_vgetA1omega(),HelWeightMinus);                                      TRFomegabar_a1_plus->Fill(a1h.TRF_vgetA1omega(),HelWeightPlus);
	TRFomegabar_a1scalar_minus->Fill(a1h.TRF_vgetA1omega("scalar"),HelWeightMinus);                  TRFomegabar_a1scalar_plus->Fill(a1h.TRF_vgetA1omega("scalar"),HelWeightPlus);
	cosbetacostheta_minus->Fill(a1h.cosbeta(),a1h.costhetaLF(),HelWeightMinus);                              cosbetacostheta_plus->Fill(a1h.cosbeta(),a1h.costhetaLF(),HelWeightPlus);
	TRFcosbetacostheta_minus->Fill(a1h.TRF_cosbeta(),a1h.costhetaLF(),HelWeightMinus);                  TRFcosbetacostheta_plus->Fill(a1h.TRF_cosbeta(),a1h.costhetaLF(),HelWeightPlus);
	
	s1->Fill((a1ospi +a1ss1pi).M() );
	s2->Fill((a1ospi +a1ss2pi).M() );
	qq->Fill((a1ospi+a1ss1pi+a1ss2pi).M());
	
	omega_a1_minus->Fill(Polarimetr.getOmegaA1Bar(),HelWeightMinus);  omega_a1_plus->Fill(Polarimetr.getOmegaA1Bar(),HelWeightPlus);
	hmag->Fill(Polarimetr.PVC().Vect().Mag());
	s1s2->Fill((a1ospi+a1ss2pi).M2(),(a1ospi+a1ss1pi).M2());
	int hel(0);
	if(HelWeightMinus==1)     hel = -1;
	if(HelWeightPlus==1)      hel =  1;
	if(HelWeightPlus==1)      a1_mcos2gamma_plus->Fill(a1h.getMoment(a1h.costhetaLF(),"c2g", 1));                 
	if(HelWeightMinus==1)     a1_mcos2gamma_minus->Fill(a1h.getMoment(a1h.costhetaLF(),"c2g", -1));                 
      }
    }
    
    //----------------------------- pairs -----------------

    
    if(JAK1==4 &&JAK2==3 ){
      if(ApplyCut){
	if(passed1 && passed2){
	  mass_pirho_plus->Fill( (tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)).M(), HelWeightPlus);
	  mass_pirho_minus->Fill((tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)).M(), HelWeightMinus);
	}
      }else{
	mass_pirho_plus->Fill( (tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)).M(), HelWeightPlus);
	mass_pirho_minus->Fill((tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)).M(), HelWeightMinus);
      }
    }
    
    if(JAK1==4 &&JAK2==2 ){
      if(ApplyCut){
	if(passed1 && passed2){
	  mass_murho_plus->Fill( (tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)).M(), HelWeightPlus);
	  mass_murho_minus->Fill((tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)).M(),HelWeightMinus);
	}
      }else{
	  mass_murho_plus->Fill( (tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)).M(), HelWeightPlus);
	  mass_murho_minus->Fill((tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)).M(),HelWeightMinus);
      }
    }
    
    if(JAK1==3 &&JAK2==2 ){
      if(ApplyCut){
	if(passed1 && passed2){
	  mupi_mass_plus->Fill( (tauandprod1.at(1)  + tauandprod2.at(1)).M(), HelWeightPlus);
	  mupi_mass_minus->Fill( (tauandprod1.at(1)  + tauandprod2.at(1)).M(),HelWeightMinus );
	}
      }else{
	  mupi_mass_plus->Fill( (tauandprod1.at(1)  + tauandprod2.at(1)).M(), HelWeightPlus);
	  mupi_mass_minus->Fill( (tauandprod1.at(1)  + tauandprod2.at(1)).M(),HelWeightMinus );
      }
    }
    
    
    if(JAK1 ==3 && JAK2 == 3){
      if(ApplyCut){
	if(passed1 && passed2){
	  if(Pi1.isConfigured() && Pi2.isConfigured()){
	    double OmPiPi=  (Pi1.getOmega() +Pi2.getOmega() )/(1 + Pi1.getOmega()*Pi2.getOmega());
	    TauPolPiPi.ConfigurePair(tauandprod1,"pion",tauandprod2,"pion");
	    Omegapipi_plus->Fill(TauPolPiPi.getCombOmegaBar(),HelWeightPlus);
	    Omegapipi_minus->Fill(TauPolPiPi.getCombOmegaBar(),HelWeightMinus);       
	    pipi_mass_plus->Fill((tauandprod1.at(1)+tauandprod2.at(1)).M(),HelWeightPlus);
	    pipi_mass_minus->Fill((tauandprod1.at(1)+tauandprod2.at(1)).M(),HelWeightMinus);
	  }
	}
      }else{
	double OmPiPi=  (Pi1.getOmega() +Pi2.getOmega() )/(1 + Pi1.getOmega()*Pi2.getOmega());
	TauPolPiPi.ConfigurePair(tauandprod1,"pion",tauandprod2,"pion");
	Omegapipi_plus->Fill(TauPolPiPi.getCombOmegaBar(),HelWeightPlus);
	Omegapipi_minus->Fill(TauPolPiPi.getCombOmegaBar(),HelWeightMinus);       
	pipi_mass_plus->Fill((tauandprod1.at(1)+tauandprod2.at(1)).M(),HelWeightPlus);
	pipi_mass_minus->Fill((tauandprod1.at(1)+tauandprod2.at(1)).M(),HelWeightMinus);
      } 
    }
    
     
     if(JAK1 ==2 && JAK2 == 3){
       if(Mu1.isConfigured() && Pi2.isConfigured()){
        double OmMuPi=  (Mu1.getOmega() +Pi2.getOmega() )/(1 + Mu1.getOmega()*Pi2.getOmega());
        TauPolMuPi.ConfigurePair(tauandprod1,"lepton",tauandprod2,"pion");
	OmegaMuPi_plus->Fill(TauPolMuPi.getCombOmega(),HelWeightPlus);
        OmegaMuPi_minus->Fill(TauPolMuPi.getCombOmega(),HelWeightMinus);       
       }
     }
     
      if(JAK1 ==2 && JAK2 == 4){
        if(Mu1.isConfigured() && Rho2.isConfigured()){
	  double OmMuRho=  (Mu1.getOmega() +Rho2.getOmega() )/(1 + Mu1.getOmega()*Rho2.getOmega());
	  TauPolMuRho.ConfigurePair(tauandprod1,"lepton",tauandprod2,"rho");
 
	  omega_murho_plus->Fill(TauPolMuRho.getCombOmegaBar(),HelWeightPlus);
	  omega_murho_minus->Fill(TauPolMuRho.getCombOmegaBar(),HelWeightMinus);    
        }
      }

      if(JAK1 ==3 && JAK2 == 4){
	if(Pi1.isConfigured() && Rho2.isConfigured()){
	  double OmPiRho=  (Pi1.getOmega() +Rho2.getOmega() )/(1 + Pi1.getOmega()*Rho2.getOmega());
	   TauPolPiRho.ConfigurePair(tauandprod1,"pion",tauandprod2,"rho");
	   
	   omega_pirho_plus->Fill(TauPolPiRho.getCombOmegaBar(),HelWeightPlus);
	   omega_pirho_minus->Fill(TauPolPiRho.getCombOmegaBar(),HelWeightMinus);       
         }
      }
      
      if(JAK1 ==3 && JAK2 == 5 &&   SubJAK2==51){
	if(Pi1.isConfigured() && A2.isConfigured()){
	  TauPolPiA1.ConfigurePair(tauandprod1,"pion",tauandprod2,"a1",taucharge2);
	   omega_a1pi_plus->Fill(TauPolPiA1.getCombOmegaBar(),HelWeightPlus);
	   omega_a1pi_minus->Fill(TauPolPiA1.getCombOmegaBar(),HelWeightMinus);       
	}
       }
      
      
      if(JAK1 ==2 && JAK2 == 5 &&   SubJAK2==51){
	if(Mu1.isConfigured() && A2.isConfigured()){
	  TauPolMuA1.ConfigurePair(tauandprod1,"lepton",tauandprod2,"a1",taucharge2);
	  omega_a1mu_plus->Fill(TauPolMuA1.getCombOmegaBar(),HelWeightPlus);
	   omega_a1mu_minus->Fill(TauPolMuA1.getCombOmegaBar(),HelWeightMinus);       
	}
      }
      

      if(JAK1 ==4 && JAK2 == 4){
	if(Rho1.isConfigured() && Rho2.isConfigured()){
	  TauPolRhoRho.ConfigurePair(tauandprod1,"rho",tauandprod2,"rho");
	  omega_rhorho_plus->Fill(TauPolRhoRho.getCombVisibleOmega(),HelWeightPlus);
	  omega_rhorho_minus->Fill(TauPolRhoRho.getCombVisibleOmega(),HelWeightMinus);       
	  mass_rhorho_plus->Fill((tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)+ tauandprod2.at(2)).M(),HelWeightPlus  );
	  mass_rhorho_minus->Fill((tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)+ tauandprod2.at(2)).M(),HelWeightMinus );
	}
      }
      
      
      // Run MC-TESTER on the event
      HepMCEvent temp_event(*HepMCEvt,false);
      MC_Analyze(&temp_event);
      
      // Print some events at the end of the run
      if(iEvent>=NumberOfEvents-5){  
	// pythia.event.list();
	HepMCEvt->print();
      }

      // Clean up HepMC event
      delete HepMCEvt;  
  }
  
  
  

  
  file->Write();
  file->Close();
  
  pythia.statistics();
  MC_Finalize();
  
  // This is an access to old FORTRAN info on generated tau sample. 
  // That is why it refers to old version number (eg. 2.7) for TAUOLA.
  //Tauola::summary();
}


