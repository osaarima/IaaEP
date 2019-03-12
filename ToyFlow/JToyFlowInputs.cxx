#include "JToyFlowInputs.h"

for(uint i = 0; i < N_VN; ++i){
	TH1D *pgr[N_VN];//[3];
	for(uint j = 0; j < 2; ++j)
		pgr[j] = (TH1D*)pff->Get(Form("Table %u/%s",i+1,pgn[j]));

	uint n = pgr[0]->GetNbinsX();
	pgr_v[i] = new TGraphErrors(n);
	for(uint ic = 0; ic < n; ++ic){
		pgr_v[i]->SetPoint(ic,
				pgr[0]->GetBinCenter(ic+1),pgr[0]->GetBinContent(ic+1));
		pgr_v[i]->SetPointError(ic,0,pgr[1]->GetBinContent(ic+1));
	}
}

for(uint i = 0; i < D_COUNT; ++i)
pgr_nch[i] = new TGraph(NC);

for(uint i = 0; i < NC; ++i){
	TGraph gr_eta(ETADST_N,etadst,etanch[i]);
	TF1 f("etach",[&](double *px, double *pp)->double{
			return gr_eta.Eval(px[0]);
			},-3.5,5.1,0);

	for(uint j = 0; j < D_COUNT; ++j){
		double nch = f.Integral(cov[j][0],cov[j][1])/(cov[j][1]-cov[j][0]);
		pgr_nch[j]->SetPoint(i,0.5*(CentBins[i]+CentBins[i+1]),nch);
	}
}

