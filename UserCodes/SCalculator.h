// -*- C++ -*-
//
// 
/**\class SCalculator.h SCalculator.cc
 Description: 
*/
//
// Original Author:  Vladimir Cherepanov 
//         Created:  Mon Sep 4 13:49:02 CET 2017
//
//

#ifndef SCalculator_h
#define SCalculator_h


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


class SCalculator {
 
 public:
  SCalculator();
  SCalculator(vector<TLorentzVector> TauAndProd);
  ~SCalculator();
  void Configure(vector<TLorentzVector> TauAndProd);
  bool isConfigured();
  void Setup(vector<TLorentzVector> TauAndProd, TLorentzVector ReferenceFrame );
  void Initialize(TLorentzVector t, TLorentzVector mu);
  std::vector<TLorentzVector> getBoosted(){return TauAndProd_LF;}
  TLorentzVector Boost(TLorentzVector pB, TLorentzVector frame);
  TVector3 Rotate(TVector3 LVec, TVector3 Rot);
  TVector3 pvec();

  //====================


 private:
 
  TVector3 polvec();
  vector<TLorentzVector> TauAndProd_LF;
  TLorentzVector TauLV;
  bool debug;
  TMatrixT<double> convertToMatrix(TVectorT<double> V);
};
#endif
