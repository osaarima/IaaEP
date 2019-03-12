#include <stdio.h>
#include <vector>
#include <TRandom.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TComplex.h>
#include "JConst.h"
#include "JToyFlowHistos.h"

typedef unsigned int uint;


int main(int argc, char **pargv){
	uint seed = argc > 1?atol(pargv[1]):1000;
	uint evtc = argc > 2?atol(pargv[2]):1000;
	printf("seed:\t%u\nevents:\t%u\n",seed,evtc);
	TFile *pfo = new TFile(argc > 3?pargv[3]:"results.root","recreate");

	// Define PDF based on F.A
	double pi = TMath::Pi();
	TF1 *pdf = new TF1("df","[0]*(1+2*[1]*cos(x-[5])+2*[2]*cos(2*(x-[6]))+2*[3]*cos(3*(x-[7]))+2*[4]*cos(4*(x-[8])))",-pi,pi);
	// initialize normalization and v1
	pdf->SetParameter(0,100.0); // v0
	pdf->SetParameter(1,0.0); // v1

	TRandom *prng = new TRandom(seed);

	JToyFlowHistos *jhistos = new JToyFlowHistos();
	jhistos->CreateHistos();

	JToyFlowInputs *jflowinputs = new JToyFlowInputs();
	jflowinputs->Load(); // dN/deta and vn as a function of centrality
	// jflowinputs->GetMultiplicity(cent);
	// jflowinputs->GetVn(cert, n);

	// Need to get Multiplicity and vn based on the measured data.
	// ALICE dN/deta vn
	double vn[2];
    // Event Loop
	for(uint evt = 0; evt < evtc; ++evt){
		//Event generation ----------------------------
		double cent = prng->Uniform(0,50.0);
		double Multiplicity = jflowinputs->GetMultiplicity(cent);
		v[0]=  jflowinputs->GetVn(cent, 2); 
		v[1]=  jflowinputs->GetVn(cent, 3);
		pdf->SetParameter(2,v[0]); //v2 = pgr_v[0]->Eval(cent);
		pdf->SetParameter(3,v[1]); //v3 [1]
		pdf->SetParameter(4,0.01); //v4

		pdf->SetParameter(5,prng->Uniform(-pi,pi));  //EP for v1
		pdf->SetParameter(6,prng->Uniform(-pi/2.0,pi/2.0)); //EP for v2
		pdf->SetParameter(7,prng->Uniform(-pi/3.0,pi/3.0)); //EP for v3
		pdf->SetParameter(8,prng->Uniform(-pi/4.0,pi/4.0)); //EP for v4

		uint cid = 0;
		for(; cid < NC; ++cid)
			if(cent < CentBins[cid+1])
				break;

		//Q-vectors -----------------------------------
    	double trueevp = pdf->GetParameter(6);

    	std::vector <double> trackphi[D_COUNT];
		TComplex Qsd[D_COUNT];
		for(uint s = 0; s < D_COUNT; ++s){
			uint ntracks = (uint)pgr_nch[s]->Eval(cent) *0.9;

			TComplex Qa2 = TComplex(0,0);
			for(uint i = 0; i < ntracks; ++i){
				double tphi = pdf->GetRandom();
        		double phi = CheckDetectorPhi(tphi);
        		trackphi[s].push_back(phi);
				Qa2 += TComplex(TMath::Cos(2.0*phi),TMath::Sin(2.0*phi));
			}
			Qa2 /= (double) (ntracks);
			Qsd[s] = Qa2/TComplex::Abs(Qa2);
		}

    	//Calculate Evp
		for (uint s = 0; s < R_COUNT; ++s) {
      		double recoevp = TMath::ATan2(Qsd[s].Im(), Qsd[s].Re())/2;
      		double evpdiff = trueevp - recoevp;
      
      		evph[s][cid]->Fill(recoevp);
      		evpcorr2d[s][cid]->Fill(trueevp, recoevp);
      		evpdifference[s][cid]->Fill(evpdiff);
    	}

		//Calculate the resolution components
		TComplex ab, ac, bc;
		ab = Qsd[D_V0A]*TComplex::Conjugate(Qsd[D_TPC_ETAA]);
		ac = Qsd[D_V0A]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
		bc = Qsd[D_TPC_ETAA]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
		pah[R_V0A][cid]->Fill(ab.Re());
		pbh[R_V0A][cid]->Fill(ac.Re());
		pch[R_V0A][cid]->Fill(bc.Re());

		ab = Qsd[D_V0C]*TComplex::Conjugate(Qsd[D_TPC_ETAA]);
		ac = Qsd[D_V0C]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
		//bc = Qsd[D_TPC_ETAA]*TComplex::Conjugate(Qsd[D_TPC]);
		pah[R_V0C][cid]->Fill(ab.Re());
		pbh[R_V0C][cid]->Fill(ac.Re());
		pch[R_V0C][cid]->Fill(bc.Re());

		ab = Qsd[D_V0P]*TComplex::Conjugate(Qsd[D_TPC_ETAA]);
		ac = Qsd[D_V0P]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
		bc = Qsd[D_TPC_ETAA]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
		pah[R_V0P][cid]->Fill(ab.Re());
		pbh[R_V0P][cid]->Fill(ac.Re());
		pch[R_V0P][cid]->Fill(bc.Re());

	} // end of event loop
	return 0;
}

