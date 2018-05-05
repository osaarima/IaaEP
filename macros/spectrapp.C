#include "include/Filipad.h"

double lowx=-1.;
double highx=100.;
double ly = 1e-8;
double hy = 1.2e1;
double lowRAA = 0.0;
double highRAA = 1.3;

TLatex latexRun;
TString strRunPP = "pp #sqrt{#it{s}} = 2.76 TeV";

int gMarkers[] = {20,24,21,25,23,27,29,30};
int gColors[]={kBlack, kRed, kBlue, kDeepSea, kPink, kGray, kRed, kBlack};
void LoadALICEPub();
void LoadPP();
void DrawPPOnly();
double sigmaPYTHIA = 65; // ??????



TGraphAsymmErrors	*gr_alice_pp_5000GeV;


const int NPP = 3;
TH1D *hChargedPtPPOnly[NPP];
TH1D *hChargedPtPPOnly_Ratio[NPP];

TFile *fPP[NPP];
int iRef = 0;
TString filePP[NPP] = {
	"legotrain_JCIaa/data/JCIaaJt_legotrain_CF_pp_MC-531_20180505-0110-LHC17l3b_cent_woSDD.root", // 5TeV pp kineonly
	"legotrain_JCIaa/data/JCIaaJt_legotrain_CF_pp_MC-532_20180505-0110-LHC17l3b_cent_woSDD.root", // 5TeV pp reco
	"legotrain_JCIaa/data/JCIaaJt_legotrain_CF_pp-1773_20180423-1806-LHC17p_pass1_CENT_woSDD.root"
};

TString dirPP[NPP] = {
	"JCIAA_KineOnly",
	"JCIAA_Reco",
	"JCIAA_GlobalSDD"
};
TString commentPP[NPP] = {
	"LHC17l3b_cent_woSDD_KineOnly",
	"LHC17l3b_cent_woSDD_Reco",
	"LHC17p_pass1_CENT_woSDD"
};


void RunPP(){
	LoadALICEPub();
	LoadPP();
	DrawPPOnly();
}

void LoadPP() {
	
	double noEventspp;
	double sc = 2*TMath::Pi()*1.6;
	TH1D *hiCentrpp;
	for(int i=0;i<NPP;i++) {
		fPP[i] = TFile::Open(filePP[i]);
		hChargedPtPPOnly[i] = (TH1D *)fPP[i]->Get(Form("%s/AliJHistos/hChargedPt/hChargedPtC%02d",dirPP[i].Data(),0)); 
		hiCentrpp = (TH1D *)fPP[i]->Get(Form("%s/AliJHistos/hiCentr",dirPP[i].Data()));
		//if(i==2 || i==3 ) hiCentrpp = (TH1D *)fPP[1]->Get(Form("%s/AliJHistos/hiCentr",dirPP[1].Data()));
		noEventspp = hiCentrpp->GetBinContent(1); // with vertex cut
		cout << i<<"\t"<< dirPP[i] << "\t"<< noEventspp << endl;
		hChargedPtPPOnly[i]->Scale(1./noEventspp/sc);
	}
}

void LoadALICEPub() {
	TFile *faliceNew = TFile::Open("spectra/alice_pp.root");
	gr_alice_pp_5000GeV = (TGraphAsymmErrors*)faliceNew->Get("gr_alice_pp_5000GeV");
}


void DrawPPOnly() {
	// pp comparison
	Filipad *fpad;
	fpad = new Filipad(1, 1.5, 0.3, 100, 100, 0.7,2);
	fpad->Draw();
	//==== Upper pad
	TPad *p = fpad->GetPad(1); //upper pad
	p->SetTickx(); p->SetLogx(1); p->SetLogy(1); p->cd();
	//hy = hChargedPtPPOnly[iRef]->GetMaximum()*1.2;
	TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
	hset( *hfr, "p_{T} [GeV/c]", "1/N_{eve}1/(2#pip_{T})dN/dp_{T}[(GeV/c)^{-2}]",1.1,1.15, 0.06,0.06, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
	//Ed^{3}#sigma/dp^{3} [mb GeV^{-2}c^{3}] ??
	hfr->Draw();
	//Legend definition
	TLegend *leg = new TLegend(0.35,0.4,0.85,0.78,"","brNDC");
	leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

	for(int i=0;i<NPP;i++) {
		hChargedPtPPOnly[i]->SetMarkerStyle(gMarkers[i]);
		hChargedPtPPOnly[i]->Draw("p,same");
		leg->AddEntry(hChargedPtPPOnly[i],commentPP[i],"pl");
		hChargedPtPPOnly_Ratio[i] = (TH1D*)hChargedPtPPOnly[i]->Clone();
		hChargedPtPPOnly_Ratio[i]->Divide(hChargedPtPPOnly[iRef]);
	}
	gr_alice_pp_5000GeV->SetMarkerStyle(31);
	gr_alice_pp_5000GeV->Draw("p,same");
	leg->AddEntry(gr_alice_pp_5000GeV,gr_alice_pp_5000GeV->GetTitle(),"pl");
	
	leg->Draw();
	// Calculation Various Ratios


	//==== Lower pad
	p = fpad->GetPad(2);
	p->SetTickx(); p->SetGridy(1); p->SetLogx(1), p->SetLogy(0); p->cd();
	TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, 0.0, 5);
	hset( *hfr1, "p_{T} [GeV/c]", "Ratio",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
	hfr1->Draw();
	for(int i=0;i<NPP;i++) {
		hChargedPtPPOnly_Ratio[i]->Draw("p,same");
	}
	gPad->GetCanvas()->SaveAs("figs_iaa/ppOnlyspectra_comp.pdf");

}


