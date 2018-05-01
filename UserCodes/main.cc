#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "TLorentzVector.h"
#include "a1Helper.h"
#include <vector>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TComplex.h"

void PrintLV(TLorentzVector a){
  std::cout<<" px:    "<< a.Px() << "  py:   "<< a.Py() <<"  pz:   "<< a.Pz() <<"  E:   "<< a.E() << "  mass:   "<< a.M()<< std::endl;
}

double momentOneP(Double_t *x, Double_t *par)
{
TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

std::vector<TLorentzVector> particles;
particles.push_back(ta1);
particles.push_back(os);
particles.push_back(ss1);
particles.push_back(ss2);

a1Helper Helper(particles,a1);


 Float_t xx =x[0];
 // Double_t f = par[0]*Helper.BreitWigner(xx).Rho2();
 Double_t f = par[0]*Helper.getMoment(xx,"one",1);


 return f;
}

double momentOneM(Double_t *x, Double_t *par)
{
TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

std::vector<TLorentzVector> particles;
particles.push_back(ta1);
particles.push_back(os);
particles.push_back(ss1);
particles.push_back(ss2);

a1Helper Helper(particles,a1);


 Float_t xx =x[0];
 // Double_t f = par[0]*Helper.BreitWigner(xx).Rho2();
 Double_t f = par[0]*Helper.getMoment(xx,"one",-1);


 return f;
}

double momentc2gM(Double_t *x, Double_t *par)
{
TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

std::vector<TLorentzVector> particles;
particles.push_back(ta1);
particles.push_back(os);
particles.push_back(ss1);
particles.push_back(ss2);

a1Helper Helper(particles,a1);


 Float_t xx =x[0];
 // Double_t f = par[0]*Helper.BreitWigner(xx).Rho2();
 Double_t f = par[0]*Helper.getMoment(xx,"c2g",-1);


 return f;
}
double momentc2gP(Double_t *x, Double_t *par)
{
TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

std::vector<TLorentzVector> particles;
particles.push_back(ta1);
particles.push_back(os);
particles.push_back(ss1);
particles.push_back(ss2);

a1Helper Helper(particles,a1);


 Float_t xx =x[0];
 // Double_t f = par[0]*Helper.BreitWigner(xx).Rho2();
 Double_t f = par[0]*Helper.getMoment(xx,"c2g",1);


 return f;
}





double momentbetaM(Double_t *x, Double_t *par)
{
TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

std::vector<TLorentzVector> particles;
particles.push_back(ta1);
particles.push_back(os);
particles.push_back(ss1);
particles.push_back(ss2);

a1Helper Helper(particles,a1);


 Float_t xx =x[0];
 // Double_t f = par[0]*Helper.BreitWigner(xx).Rho2();
 Double_t f = par[0]*Helper.getMoment(xx,"beta",-1);


 return f;
}



double momentbetaP(Double_t *x, Double_t *par)
{
TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

std::vector<TLorentzVector> particles;
particles.push_back(ta1);
particles.push_back(os);
particles.push_back(ss1);
particles.push_back(ss2);

a1Helper Helper(particles,a1);


 Float_t xx =x[0];
 // Double_t f = par[0]*Helper.BreitWigner(xx).Rho2();
 Double_t f = par[0]*Helper.getMoment(xx,"beta",1);


 return f;
}



double momentcbM(Double_t *x, Double_t *par)
{
TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

std::vector<TLorentzVector> particles;
particles.push_back(ta1);
particles.push_back(os);
particles.push_back(ss1);
particles.push_back(ss2);

a1Helper Helper(particles,a1);


 Float_t xx =x[0];
 // Double_t f = par[0]*Helper.BreitWigner(xx).Rho2();
 Double_t f = par[0]*Helper.getMoment(xx,"cb",-1);


 return f;
}

double momentcbP(Double_t *x, Double_t *par)
{
TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

std::vector<TLorentzVector> particles;
particles.push_back(ta1);
particles.push_back(os);
particles.push_back(ss1);
particles.push_back(ss2);

a1Helper Helper(particles,a1);


 Float_t xx =x[0];
 // Double_t f = par[0]*Helper.BreitWigner(xx).Rho2();
 Double_t f = par[0]*Helper.getMoment(xx,"cb",1);


 return f;
}





double myfunctionBRWGRho(Double_t *x, Double_t *par)
{
TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

std::vector<TLorentzVector> particles;
particles.push_back(ta1);
particles.push_back(os);
particles.push_back(ss1);
particles.push_back(ss2);

a1Helper Helper(particles,a1);


 Float_t xx =x[0];
 // Double_t f = par[0]*Helper.BreitWigner(xx).Rho2();
 Double_t f = par[0]*Helper.BRho(xx).Rho2();


 return f;
}

double myfunctionBRWGA1(Double_t *x, Double_t *par)
{
TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

std::vector<TLorentzVector> particles;
particles.push_back(ta1);
particles.push_back(os);
particles.push_back(ss1);
particles.push_back(ss2);

a1Helper Helper(particles,a1);


 Float_t xx =x[0];
  Double_t f = par[0]*Helper.BreitWigner(xx,"a1").Rho2();

 return f;
}

double myfunctionWA(Double_t *x, Double_t *par)
{
TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

std::vector<TLorentzVector> particles;
particles.push_back(ta1);
particles.push_back(os);
particles.push_back(ss1);
particles.push_back(ss2);

a1Helper Helper(particles,a1);


 Float_t xx =x[0];
 Double_t f = par[0]*Helper.MomentSFunction(xx,"WA");
 return f;
}


double myfunctionWCWA(Double_t *x, Double_t *par)
{
TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

std::vector<TLorentzVector> particles;
particles.push_back(ta1);
particles.push_back(os);
particles.push_back(ss1);
particles.push_back(ss2);

a1Helper Helper(particles,a1);


 Float_t xx =x[0];
 Double_t f = par[0]*Helper.MomentSFunction(xx,"WC")/Helper.MomentSFunction(xx,"WA");
 return f;
}

double myfunctionWSDWA(Double_t *x, Double_t *par)
{
TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

std::vector<TLorentzVector> particles;
particles.push_back(ta1);
particles.push_back(os);
particles.push_back(ss1);
particles.push_back(ss2);

a1Helper Helper(particles,a1);


 Float_t xx =x[0];
 Double_t f = par[0]*Helper.MomentSFunction(xx,"WSD")/Helper.MomentSFunction(xx,"WA");
 return f;
}

double myfunctionWSBWA(Double_t *x, Double_t *par)
{
TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

std::vector<TLorentzVector> particles;
particles.push_back(ta1);
particles.push_back(os);
particles.push_back(ss1);
particles.push_back(ss2);

a1Helper Helper(particles,a1);


 Float_t xx =x[0];
 Double_t f = par[0]*Helper.MomentSFunction(xx,"WSB")/Helper.MomentSFunction(xx,"WA");
 return f;
}
double myfunctionWSAWA(Double_t *x, Double_t *par)
{
TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

std::vector<TLorentzVector> particles;
particles.push_back(ta1);
particles.push_back(os);
particles.push_back(ss1);
particles.push_back(ss2);

a1Helper Helper(particles,a1);


 Float_t xx =x[0];
 Double_t f = par[0]*Helper.MomentSFunction(xx,"WSA")/Helper.MomentSFunction(xx,"WA");
 return f;
}



double myfunctionWEWA(Double_t *x, Double_t *par)
{
TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

std::vector<TLorentzVector> particles;
particles.push_back(ta1);
particles.push_back(os);
particles.push_back(ss1);
particles.push_back(ss2);

a1Helper Helper(particles,a1);


 Float_t xx =x[0];
 Double_t f = par[0]*Helper.MomentSFunction(xx,"WE")/Helper.MomentSFunction(xx,"WA");
 return f;
}


double myfunctionWDWA(Double_t *x, Double_t *par)
{
TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

std::vector<TLorentzVector> particles;
particles.push_back(ta1);
particles.push_back(os);
particles.push_back(ss1);
particles.push_back(ss2);

a1Helper Helper(particles,a1);


 Float_t xx =x[0];
 Double_t f = par[0]*Helper.MomentSFunction(xx,"WD")/Helper.MomentSFunction(xx,"WA");
 return f;
}
int main(int argc, const char* argv[]) {
  TFile *file = new TFile("Figures.root","RECREATE");
  


  // TLorentzVector rho(-13.491274,-4.967516,1.272212,14.454947);
  // TLorentzVector pi0(-8.682081,-2.818832,0.750293,9.160063) ;
  // TLorentzVector pi(-4.809193,-2.148683,0.521918,5.294883);

  // TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
  // TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
  // TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
  // TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
  // TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);

  // std::vector<TLorentzVector> particles;
  // particles.push_back(ta1);
  // particles.push_back(os);
  // particles.push_back(ss1);
  // particles.push_back(ss2);

  // a1Helper Helper(particles,a1);

  // std::cout<<" Helper.MomentSFunction(xx; " << Helper.MomentSFunction(0.45,"WA") << std::endl;

    TF1 *f1brwga1 = new TF1("myfuncBRWGA1",myfunctionBRWGA1,0.45,6,1);
    f1brwga1->SetParameter(0,1);
    TH1 * f1hbrwga1 = f1brwga1->GetHistogram();
    f1hbrwga1->SetName("brwga1");
    f1hbrwga1->Write();





    TF1 *f1brwgrho = new TF1("myfuncBRWGRho",myfunctionBRWGRho,0.45,6,1);
    f1brwgrho->SetParameter(0,1);
    TH1 * f1hbrwgrho = f1brwgrho->GetHistogram();
    f1hbrwgrho->SetName("brwgrho");
    f1hbrwgrho->Write();



    TF1 *f1wa = new TF1("myfuncWA",myfunctionWA,0.45,6,1);
     f1wa->SetParameter(0,1);
    TH1 * f1hwa = f1wa->GetHistogram();
    f1hwa->SetName("wa");
    f1hwa->Write();



    TF1 *f1wcwa = new TF1("myfuncWCWA",myfunctionWCWA,0.45,6,1);
    f1wcwa->SetParameter(0,1);
    TH1 * f1hwcwa = f1wcwa->GetHistogram();
    f1hwcwa->SetName("wcwa");
    f1hwcwa->Write();

    TF1 *f1wdwa = new TF1("myfuncWDWA",myfunctionWDWA,0.45,6,1);
    f1wdwa->SetParameter(0,1);
    TH1 * f1hwdwa = f1wdwa->GetHistogram();
    f1hwdwa->SetName("wdwa");
    f1hwdwa->Write();

    TF1 *f1wewa = new TF1("myfuncWEWA",myfunctionWEWA,0.45,6,1);
    f1wewa->SetParameter(0,1);
    TH1 * f1hwewa = f1wewa->GetHistogram();
    f1hwewa->SetName("wewa");
    f1hwewa->Write();

    TF1 *f1wsdwa = new TF1("myfuncWSDWA",myfunctionWSDWA,0.45,6,1);
    f1wsdwa->SetParameter(0,1);
    TH1 * f1hwsdwa = f1wsdwa->GetHistogram();
    f1hwsdwa->SetName("wsdwa");
    f1hwsdwa->Write();
    TF1 *f1wsbwa = new TF1("myfuncWSBWA",myfunctionWSBWA,0.45,6,1);
    f1wsbwa->SetParameter(0,1);
    TH1 * f1hwsbwa = f1wsbwa->GetHistogram();
    f1hwsbwa->SetName("wsbwa");
    f1hwsbwa->Write();

    TF1 *f1wsawa = new TF1("myfuncWSAWA",myfunctionWSAWA,0.45,6,1);
    f1wsawa->SetParameter(0,1);
    TH1 * f1hwsawa = f1wsawa->GetHistogram();
    f1hwsawa->SetName("wsawa");
    f1hwsawa->Write();

    //-----------------------------------------------------------------------------------------------------------------
    TF1 *fmoneP = new TF1("momentOneP",momentOneP,-1,1,1);
    fmoneP->SetParameter(0,1);
    TH1 * hfmoneP = fmoneP->GetHistogram();
    hfmoneP->SetName("MomentOneP");
    hfmoneP ->Write();

     //-----------------------------------------------------------------------------------------------------------------
    TF1 *fmoneM = new TF1("momentOneM",momentOneM,-1,1,1);
    fmoneM->SetParameter(0,1);
    TH1 * hfmoneM = fmoneM->GetHistogram();
    hfmoneM->SetName("MomentOneM");
    hfmoneM ->Write();


     //-----------------------------------------------------------------------------------------------------------------
    TF1 *fmc2gM = new TF1("momentc2gM",momentc2gM,-1,1,1);
    fmc2gM->SetParameter(0,1);
    TH1 * hfmc2gM = fmc2gM->GetHistogram();
    hfmc2gM->SetName("Momentc2gM");
    hfmc2gM ->Write();

    //-----------------------------------------------------------------------------------------------------------------
    TF1 *fmc2gP = new TF1("momentc2gP",momentc2gP,-1,1,1);
    fmc2gP->SetParameter(0,1);
    TH1 * hfmc2gP = fmc2gP->GetHistogram();
    hfmc2gP->SetName("Momentc2gP");
    hfmc2gP ->Write();

     //-----------------------------------------------------------------------------------------------------------------
    TF1 *fmbetaM = new TF1("momentbetaM",momentbetaM,-1,1,1);
    fmbetaM->SetParameter(0,1);
    TH1 * hfmbetaM = fmbetaM->GetHistogram();
    hfmbetaM->SetName("MomentbetaM");
    hfmbetaM ->Write();

    //-----------------------------------------------------------------------------------------------------------------
    TF1 *fmbetaP = new TF1("momentbetaP",momentbetaP,-1,1,1);
    fmbetaP->SetParameter(0,1);
    TH1 * hfmbetaP = fmbetaP->GetHistogram();
    hfmbetaP->SetName("MomentbetaP");
    hfmbetaP ->Write();

    //-----------------------------------------------------------------------------------------------------------------
    TF1 *fmcbM = new TF1("momentcbM",momentcbM,-1,1,1);
    fmcbM->SetParameter(0,1);
    TH1 * hfmcbM = fmcbM->GetHistogram();
    hfmcbM->SetName("MomentcbM");
    hfmcbM ->Write();

    //-----------------------------------------------------------------------------------------------------------------
    TF1 *fmcbP = new TF1("momentcbP",momentcbP,-1,1,1);
    fmcbP->SetParameter(0,1);
    TH1 * hfmcbP = fmcbP->GetHistogram();
    hfmcbP->SetName("MomentcbP");
    hfmcbP ->Write();




  file->Write();
  file->Close();

}














// ---------
// (x,y,z,t)=(-26.788939,2.687771,-37.247365,45.993420) (P,eta,phi,E)=(45.959086,-1.128328,3.041596,45.993420)
// (x,y,z,t)=(-15.012722,1.198184,-21.891933,26.584986) (P,eta,phi,E)=(26.572057,-1.168748,3.061950,26.584986)
// (x,y,z,t)=(-3.664335,0.233661,-5.183266,6.353556) (P,eta,phi,E)=(6.352023,-1.144735,3.077913,6.353556)
// (x,y,z,t)=(-7.876100,0.720873,-12.047975,14.412696) (P,eta,phi,E)=(14.412020,-1.207630,3.050320,14.412696)
// (x,y,z,t)=(-3.472285,0.243649,-4.660688,5.818730) (P,eta,phi,E)=(5.817056,-1.101985,3.071538,5.818730)
// ---------
// (x,y,z,t)=(-32.822804,29.511723,-12.752927,45.979044) (P,eta,phi,E)=(45.944700,-0.285049,2.409263,45.979044)
// (x,y,z,t)=(-16.949156,15.002583,-7.017585,23.721868) (P,eta,phi,E)=(23.698056,-0.305267,2.417042,23.721868)
// (x,y,z,t)=(-5.037089,4.346876,-2.412471,7.078636) (P,eta,phi,E)=(7.077260,-0.355084,2.429615,7.078636)
// (x,y,z,t)=(-0.923205,0.855447,-0.441767,1.341169) (P,eta,phi,E)=(1.333887,-0.344161,2.394271,1.341169)
// (x,y,z,t)=(-10.988859,9.800258,-4.163346,15.302059) (P,eta,phi,E)=(15.301423,-0.279118,2.413307,15.302059)
// ---------
// (x,y,z,t)=(-38.702364,-21.588178,10.946964,45.682771) (P,eta,phi,E)=(45.648203,0.244574,-2.632781,45.682771)
// (x,y,z,t)=(-23.317578,-13.107163,7.037263,27.685411) (P,eta,phi,E)=(27.659180,0.260141,-2.629496,27.685411)
// (x,y,z,t)=(-11.382100,-6.650994,3.389942,13.612462) (P,eta,phi,E)=(13.611746,0.254395,-2.612769,13.612462)
// (x,y,z,t)=(-3.030154,-1.805861,1.235359,3.740128) (P,eta,phi,E)=(3.737523,0.343422,-2.604147,3.740128)
// (x,y,z,t)=(-8.905323,-4.650308,2.411961,10.332820) (P,eta,phi,E)=(10.331878,0.237834,-2.660348,10.332820)
// ---------
// (x,y,z,t)=(16.862003,42.748406,0.932397,45.997607) (P,eta,phi,E)=(45.963276,0.020288,1.195086,45.997607)
// (x,y,z,t)=(13.909979,36.136450,1.181131,38.752321) (P,eta,phi,E)=(38.739200,0.030499,1.203349,38.752321)
// (x,y,z,t)=(3.471275,9.245016,0.396963,9.884186) (P,eta,phi,E)=(9.883201,0.040187,1.211609,9.884186)
// (x,y,z,t)=(5.788296,15.810430,0.630423,16.849064) (P,eta,phi,E)=(16.848486,0.037435,1.219846,16.849064)
// (x,y,z,t)=(4.650407,11.081000,0.153744,12.019067) (P,eta,phi,E)=(12.018256,0.012793,1.173446,12.019067)
// ---------
// (x,y,z,t)=(-18.658958,41.971850,1.745631,45.999970) (P,eta,phi,E)=(45.965641,0.037995,1.989116,45.999970)
// (x,y,z,t)=(-18.431158,41.073397,1.637837,45.066871) (P,eta,phi,E)=(45.049018,0.036373,1.992600,45.066871)
// (x,y,z,t)=(-8.572447,19.482827,0.453507,21.290668) (P,eta,phi,E)=(21.290210,0.021304,1.985303,21.290668)
// (x,y,z,t)=(-2.981254,7.248164,0.512469,7.855308) (P,eta,phi,E)=(7.854068,0.065342,1.961016,7.855308)
// (x,y,z,t)=(-6.877456,14.342405,0.671861,15.920894) (P,eta,phi,E)=(15.920282,0.042227,2.017925,15.920894)
// ---------
// (x,y,z,t)=(-21.363342,-31.640470,-25.580263,45.989332) (P,eta,phi,E)=(45.954995,-0.627947,-2.164677,45.989332)
// (x,y,z,t)=(-8.665741,-13.056724,-10.431898,18.858951) (P,eta,phi,E)=(18.825451,-0.624333,-2.156742,18.858951)
// (x,y,z,t)=(-0.701739,-1.383549,-0.893003,1.795432) (P,eta,phi,E)=(1.789999,-0.547820,-2.040189,1.795432)
// (x,y,z,t)=(-2.848837,-4.717945,-3.658006,6.616295) (P,eta,phi,E)=(6.614823,-0.622695,-2.114027,6.616295)
// (x,y,z,t)=(-5.115166,-6.955231,-5.880890,10.447225) (P,eta,phi,E)=(10.446292,-0.637162,-2.204915,10.447225)
// ---------

