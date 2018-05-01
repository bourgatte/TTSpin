#ifndef TauDecaysHelper_h
#define TauDecaysHelper_h
#include <vector>
#include "TLorentzVector.h"
#include "TComplex.h"
#include "TMatrixT.h"
#include "TVectorT.h"
#include "TMatrixTSym.h"
#include <string.h>
#include <vector> 
#include "TLorentzVector.h"
using namespace std;
class TauDecaysHelper {
 public:
  TauDecaysHelper();
  TauDecaysHelper(vector<TLorentzVector> TauAndProd, string type);
  ~TauDecaysHelper();

  void Configure(vector<TLorentzVector> TauAndProd, string type );
  bool  isConfigured();
  TLorentzVector Boost(TLorentzVector pB, TLorentzVector frame);
  double getOmega();
  double getOmegabar();
  double getCosbetaRho();
  double getSinbetaRho();
  double getCosthetaRho();
  double getSinthetaRho();

  double getUltrarel_cospsiLF();
  double getSinpsiLF();

  



  TLorentzVector sLV(); 
  TVector3 nPerp();
  TVector3 ns(); 
  TVector3 nL(); 
  TVector3 nT();
  double  DPF_cosalpha(); 
  double  DPF_sinalpha(); 



 private:
  double mrho;
  double mpi;
  double mtau;
  double ma1;
  bool debug;
  TMatrixT<double> convertToMatrix(TVectorT<double> V);
  TLorentzVector TauLV;
  TLorentzVector ProductLV;
  TLorentzVector TauRhoPi;
  TLorentzVector TauRhoPi0;
  TLorentzVector InvisibleLV;

  TLorentzVector DPF_TauLV;
  TLorentzVector DPF_TauRhoPi;
  TLorentzVector DPF_TauRhoPi0;
  TLorentzVector DPF_InvisibleLV;

  string type_;
};
#endif
