/**
 * Example of use of tauola C++ interface. Pythia events are
 * generated with a stable tau. Taus are subsequently decay via
 * tauola.
 *
 * @author Nadia Davidson
 * @date 17 June 2008
 */

#include "Tauola/Log.h"
#include "Tauola/Plots.h"
#include "Tauola/Tauola.h"
#include "Tauola/TauolaHepMCEvent.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h" 
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
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

using namespace std;
using namespace Pythia8; 
using namespace Tauolapp;

unsigned int NumberOfEvents = 1000;
unsigned int EventsToCheck  = 20;

// elementary test of HepMC typically executed before
// detector simulation based on http://home.fnal.gov/~mrenna/HCPSS/HCPSShepmc.html
// similar test was performed in Fortran
// we perform it before and after Tauola (for the first several events)
void checkMomentumConservationInEvent(HepMC::GenEvent *evt)
{
	//cout<<"List of stable particles: "<<endl;

	double px=0.0,py=0.0,pz=0.0,e=0.0;
	
	for ( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin();
	      p != evt->particles_end(); ++p )
	{
		if( (*p)->status() == 1 )
		{
			HepMC::FourVector m = (*p)->momentum();
			px+=m.px();
			py+=m.py();
			pz+=m.pz();
			e +=m.e();
			//(*p)->print();
		}
	}
  cout.precision(6);
  cout.setf(ios_base::floatfield);
	cout<<endl<<"Vector Sum: "<<px<<" "<<py<<" "<<pz<<" "<<e<<endl;
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

int main(int argc,char **argv){

 
  TString FileName  = "TTSpin_plots.root";
  TFile *file = new TFile(FileName,"RECREATE");
  TH1F *phiang= new TH1F("phiang","#phi ",50,0,3.14);
  // Program needs at least 4 parameters
  if(argc<5)
  {
    cout<<endl<<"Usage: "<<argv[0]<<" <pythia_conf> <pythia_mode> <no_events> <tauola_mode> <mixing_angle>"<<endl;
    cout<<endl<<"   eg. "<<argv[0]<<" pythia_H.conf 0 10000 4 0.7853"<<endl;
    cout<<endl;
    return -1;
  }

  Log::SummaryAtExit();

  // Initialisation of pythia
  Pythia pythia;
  Event& event = pythia.event;
  
  // Pythia8 HepMC interface depends on Pythia8 version
#ifdef PYTHIA8180_OR_LATER
  HepMC::Pythia8ToHepMC ToHepMC;
#else
  HepMC::I_Pythia8 ToHepMC;
  // we keep it for backward compatibility
  //pythia.readString("HadronLevel:all = off");
  pythia.readString("HadronLevel:Hadronize = off");
  pythia.readString("SpaceShower:QEDshowerByL = off");
  pythia.readString("SpaceShower:QEDshowerByQ = off");
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = off");

#endif

  // Initial pythia configuration
  pythia.particleData.readString("15:mayDecay = off");

  /*
    Read input parameters from console. List of parameters:
    1. Pythia configuration filename
    2. Are we using pp collisions? (If not - e+ e- collisions)
    3. Number of events
    4. Tauola decay mode (refer to documentation)
    5. Higgs scalar-pseudoscalar mixing angle

    Example where all input parameters are used:

    ./taumain_pythia_example.exe pythia_H.conf 0 100000 4 0.7853
      - use pythia_H.conf
      - initialize using e+ e- collisions
      - generate 100 000 events
      - fix TAUOLA decay to channel 4 (RHO_MODE)
      - Higgs scalar-pseudoscalar mixing angle set to 0.7853
  */

  // 1. Load pythia configuration file (argv[1], from console)
  if(argc>1) pythia.readFile(argv[1]);

  // // 2. Initialize pythia to pp or e+e- collisions (argv[2], from console)
   if(atoi(argv[2])==1) pythia.init( -2212, -2212, 14000.0); // p_bar  p_bar  collisions
   else                 pythia.init( 11, -11, 500);          // e+ e- collisions
  //pythia.init( -2212, -2212, 14000.0);
  // 3. Get number of events (argv[3], from console)
  if(argc>3) NumberOfEvents=atoi(argv[3]);

  // 4. Set Tauola decay mode (argv[4], from console)
  if(argc>4)
  {
    Tauola::setSameParticleDecayMode(atoi(argv[4]));
    Tauola::setOppositeParticleDecayMode(atoi(argv[4]));
  }

  // 5. Set Higgs scalar-pseudoscalar mixing angle (argv[5], from console)
  if(argc>5)
  {
    Tauola::setHiggsScalarPseudoscalarMixingAngle(atof(argv[5]));
    Tauola::setHiggsScalarPseudoscalarPDG(25);
  }

  Tauola::setRadiation(false); // turn off radiation in leptionic decays
  Tauola::initialize();

  MC_Initialize();

  // Begin event loop
  for(unsigned int iEvent = 0; iEvent < NumberOfEvents; ++iEvent)
  {
    if (iEvent%1000==0) Log::Info()<<"Event: "<<iEvent<<endl;
    if (!pythia.next()) continue;

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

    // We do not let Pythia decay taus, so we don't have to undecay them
    //t_event->undecayTaus();
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
    TLorentzVector tau1(0,0,0,0);
    TLorentzVector tau2(0,0,0,0);
    vector<TLorentzVector> tauandprod1,tauandprod2;
    tau1.SetPxPyPzE(FirstTau->momentum().px(), FirstTau->momentum().py(), FirstTau->momentum().pz(), FirstTau->momentum().e());
    tau2.SetPxPyPzE(SecondTau->momentum().px(), SecondTau->momentum().py(), SecondTau->momentum().pz(), SecondTau->momentum().e());



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
    }
    
    if(JAK1==3){
      vector<TLorentzVector> tauandprod;
      tauandprod.push_back(TLorentzVector(FirstTau->momentum().px(), FirstTau->momentum().py(), FirstTau->momentum().pz(), FirstTau->momentum().e()));
      for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
	if(abs(a->pdg_id())==211){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );}}
      tauandprod1 = tauandprod;
    }
    
    if(JAK1==4){
      vector<TLorentzVector> tauandprod;
      tauandprod.push_back(TLorentzVector(FirstTau->momentum().px(), FirstTau->momentum().py(), FirstTau->momentum().pz(), FirstTau->momentum().e()));
      for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
	if(abs(a->pdg_id())==211){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );}
	if(abs(a->pdg_id())==111){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );}
      }
      tauandprod1=tauandprod;  
    }
    
    if(JAK1==5 && SubJAK1==51){
      vector<TLorentzVector> particles;
      particles.clear();
      SortPions(A1Pions1);
      int taucharge =  (A1Pions1.at(0).pdg_id()+A1Pions1.at(1).pdg_id()+A1Pions1.at(2).pdg_id() > 0) ? 1 : -1;
      //      taucharge1=taucharge;
      tauandprod1=particles;
    }
    
     //--------------------  tau+
     if(JAK2==2){
       vector<TLorentzVector> tauandprodmu2;
       tauandprodmu2.push_back(TLorentzVector(SecondTau->momentum().px(), SecondTau->momentum().py(), SecondTau->momentum().pz(), SecondTau->momentum().e()));
       for(std::vector<HepMC::GenParticle>::const_iterator a = SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
     	if(abs(a->pdg_id())==13){tauandprodmu2.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );} }
       tauandprod2=tauandprodmu2;
     }
     
     if(JAK2==3){
       vector<TLorentzVector> tauandprod;
       tauandprod.push_back(TLorentzVector(SecondTau->momentum().px(), SecondTau->momentum().py(), SecondTau->momentum().pz(), SecondTau->momentum().e()));
       for(std::vector<HepMC::GenParticle>::const_iterator a = SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
	 if(abs(a->pdg_id())==211){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );} }
       tauandprod2 = tauandprod;
     }
     
     if(JAK2==4){
       vector<TLorentzVector> tauandprod;
       tauandprod.push_back(TLorentzVector(SecondTau->momentum().px(), SecondTau->momentum().py(), SecondTau->momentum().pz(), SecondTau->momentum().e()));
       for(std::vector<HepMC::GenParticle>::const_iterator a =SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
	 if(abs(a->pdg_id())==211){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );}
	 if(abs(a->pdg_id())==111){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );}
       }
       tauandprod2=tauandprod;  
     }

    
    if(JAK2==5 && SubJAK2==51){
      vector<TLorentzVector> particles;
      particles.clear();
      SortPions(A1Pions2);
      
      int taucharge =  (A1Pions2.at(0).pdg_id()+A1Pions2.at(1).pdg_id()+A1Pions2.at(2).pdg_id() > 0) ? 1 : -1;
      //      taucharge2=taucharge;
      tauandprod2=particles;
    }




    if(JAK1 ==3 && JAK2 == 3){

      TLorentzVector tau1Lab = tauandprod1.at(0);
      TLorentzVector tau2Lab = tauandprod2.at(0);
      TLorentzVector ttRF = tau1Lab + tau2Lab;
      TVector3 Rot1 = tau1Lab.Vect();
      TLorentzVector HiigsLV = tau1Lab + tau2Lab;
      //      std::cout<<"  mnass    "<< ( tauandprod1.at(0) +tauandprod2.at(0) ).M() <<std::endl;

       TLorentzVector pi1 = tauandprod1.at(1);
       TLorentzVector pi2 = tauandprod2.at(1);

       TLorentzVector pi1R1 = pi1;
       TLorentzVector pi2R1 = pi2;
       TLorentzVector tau1LabR1 = tau1Lab;
       TLorentzVector tau2LabR1 = tau2Lab;

       tau1LabR1.SetVect(Rotate(tau1LabR1.Vect(),Rot1));
       tau2LabR1.SetVect(Rotate(tau2LabR1.Vect(),Rot1));

       pi2R1.SetVect(Rotate(pi2R1.Vect(),Rot1));
       pi1R1.SetVect(Rotate(pi1R1.Vect(),Rot1));


       TLorentzVector Tau1HRF = BoostR(tau1Lab,HiigsLV);
       TLorentzVector Tau2HRF = BoostR(tau2Lab,HiigsLV);

       TLorentzVector Pi1HRF = BoostR(pi1,HiigsLV);
       TLorentzVector Pi2HRF = BoostR(pi2,HiigsLV);


       TVector3 n1 = Pi1HRF.Vect().Cross(Tau1HRF.Vect());
       TVector3 n2 = Pi2HRF.Vect().Cross(Tau2HRF.Vect());
 
          if(n1.Mag()!=0 && n2.Mag()!=0){
            double mag = n1.Mag()*n2.Mag();
            double thetaangle = acos(n1*n2/mag);
	    //	    std::cout<<"  theta angle "<< thetaangle <<std::endl;

	    phiang->Fill(thetaangle);

	  }


       // std::cout<<"  TT "<< polpi1.X()*polpi2.X() + polpi1.Y()*polpi2.Y()  <<std::endl;
       // std::cout<<"  TN "<< polpi1.X()*polpi2.Y() + polpi1.Y()*polpi2.X()  <<std::endl;
       // phiang->Fill(acos(v1*v2));
}




    // Run MC-TESTER on the event
    HepMCEvent temp_event(*HepMCEvt,false);
    MC_Analyze(&temp_event);

    // Clean up HepMC event
    delete HepMCEvt;
  }
  pythia.statistics();
  Tauola::summary();
  MC_Finalize();

  
  file->Write();
  file->Close();


}

