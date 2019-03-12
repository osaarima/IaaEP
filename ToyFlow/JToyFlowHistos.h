//===========================================================
// JToyFlowHistos.h
//===========================================================

#ifndef JTOYFLOWHISTOS_H
#define JToyFLOWHISTOS_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TFile.h>
#include "JConst.h"


class JToyFlowHistos {

public:
  JToyFlowHistos(); 
  virtual ~JToyFlowHistos(){;}	  //destructor
  
  // ALICE methods =====================================================
  void CreateHistos();

private:
  char  hname[40], htit[40];
  float b1, b2, pb1, pb2;

public:
  //===================================================
  // Toy Flow histograms
  //===================================================
  TH1D *pah[R_COUNT][NC], *pbh[R_COUNT][NC], *pch[R_COUNT][NC];
  TH1D *evph[R_COUNT][NC];
  TH2D *contami2d[R_COUNT][NC];
  TH2D *highcontami2d[R_COUNT][NC];
  TH2D *evpcorr2d[R_COUNT][NC];
  TH2D *evpcorrvsdet2d[NC];
  TH1D *evpdifference[R_COUNT][NC];

};

#endif


