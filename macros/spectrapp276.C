#include "include/Filipad.h"

double lowx=-0.1;
double highx=100.;
double ly = 1e-9;
double hy = 2.2e3;
double lowRAA = -0.2;
double highRAA = 5;

TLatex latexRun;
TString strRunPP = "pp #sqrt{#it{s}} = 2.76 TeV";

int gMarkers[] = {20,24,21,25,23,27,29,30};
int gColors[]={kBlack, kRed, kBlue, kDeepSea, kPink, kGray, kRed, kBlack};
void LoadALICEPub();
void LoadPP();
void DrawPPOnly();
enum energy {k900,k2760,k7000, kPb};
double triggEff[4] = {0.91, 0.881,0.852,1.0}; // MB_AND http://arxiv.org/pdf/1208.4968v1.pdf table 3.
double vtxEff[4] = {0.763, 0.9,0.9,1.0}; // MB_AND http://arxiv.org/pdf/1208.4968v1.pdf

double sigmaData[2] = { 62.2, 73.2 };
double sigmaPYTHIA[2] = { 62.04, 71.39 };

TH1D *hChargedPtpp_ALICEpub;
TGraphAsymmErrors	*gr_alice_pp_2760GeV;



const int NPP = 3;
TH1D *hChargedPtPPOnly[NPP];
TH1D *hChargedPtPPOnly_Ratio[NPP];

TFile *fPP[NPP];
int iRef = 1;
TString filePP[NPP] = {
	"legotrain_JCIaa/data/JCIaa_legotrain_CF_pp-1708_20180405-0222-2760GeV_LHC11a_p4_AOD113_noSDD.root",
	"legotrain_JCIaa/data/JCIaaJt_legotrain_CF_pp_MC-512_20180427-1104-LHC12f1a_Pythia_2760GeV.root", // 2.76TeV pp
	"legotrain_JCIaa/data/JCIaa_legotrain_CF_pp_MC-521_20180503-1035-LHC12f1a_Pythia_2760GeV.root"
};

TString dirPP[NPP] = {
  "JCIAA_GlobalSDD_H0_T0",
  "JCIAA_KineOnly",
  "JCIAA_Reco"
};
TString commentPP[NPP] = {
	"LHC11a_p4_AOD113_noSDD",
	"LHC12f1a_Pythia_2760GeV_KineOnly",
	"LHC12f1a_Pythia_2760GeV_Reco"
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
		hChargedPtPPOnly[i] = (TH1D *)fPP[i]->Get(Form("%s/AliJHistos/hChargedPtPublished/hChargedPtPublishedC%02d",dirPP[i].Data(),0)); 
		hiCentrpp = (TH1D *)fPP[0]->Get(Form("%s/AliJHistos/hiCentr",dirPP[0].Data()));
		//if(i==2 || i==3 ) hiCentrpp = (TH1D *)fPP[1]->Get(Form("%s/AliJHistos/hiCentr",dirPP[1].Data()));
		noEventspp = hiCentrpp->GetBinContent(1); // with vertex cut
		cout << i<<"\t"<< dirPP[i] << "\t"<< noEventspp << endl;
		//hChargedPtPPOnly[i]->Scale(sigmaData[0],"width");
		hChargedPtPPOnly[i]->Scale(1.,"width");
		double norm =   triggEff[k2760] * vtxEff[k2760] / sc; 
		hChargedPtPPOnly[i]->Scale(norm/noEventspp);
	}
}

void LoadALICEPub() {
	TString infileALICE = "jacek_Feb_14_2013/2012-11-16_dNdPt_all/dNdPt_pp_2.76TeV.root";
	TFile *finA = TFile::Open(infileALICE);
	hChargedPtpp_ALICEpub = (TH1D *)finA->Get("dNdPt");
	//hChargedPtpp_ALICEpub->Scale(sigmaData[0]);
	TFile *faliceNew = TFile::Open("spectra/alice_pp.root");
	gr_alice_pp_2760GeV = (TGraphAsymmErrors*)faliceNew->Get("gr_alice_pp_2760GeV");
}


void DrawPPOnly() {
	// pp comparison
	Filipad *fpad;
	fpad = new Filipad(1, 1.1, 0.3, 100, 100, 0.7,2);
	fpad->Draw();
	//==== Upper pad
	TPad *p = fpad->GetPad(1); //upper pad
	p->SetTickx(); p->SetLogx(1); p->SetLogy(1); p->cd();
	hy = hChargedPtPPOnly[iRef]->GetMaximum()*1.2;
	TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
	hset( *hfr, "p_{T} [GeV/c]", "1/N_{eve}1/(2#pip_{T})dN/dp_{T}[(GeV/c)^{-2}]",1.1,1.2, 0.06,0.06, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
	//hset( *hfr, "p_{T} [GeV/c]", "Ed^{3}#sigma/dp^{3} [mb GeV^{-2}c^{3}]",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);
	//Ed^{3}#sigma/dp^{3} [mb GeV^{-2}c^{3}] ??
	hfr->Draw();
	//Legend definition
	TLegend *leg = new TLegend(0.45,0.6,0.8,0.85,"","brNDC");
	leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

	for(int i=0;i<NPP;i++) {
		hChargedPtPPOnly[i]->SetMarkerStyle(gMarkers[i]);
		hChargedPtPPOnly[i]->SetMarkerColor(gColors[i]);
		hChargedPtPPOnly[i]->Draw("p,same");
		leg->AddEntry(hChargedPtPPOnly[i],commentPP[i],"pl");
		hChargedPtPPOnly_Ratio[i] = (TH1D*)hChargedPtPPOnly[i]->Clone();
		hChargedPtPPOnly_Ratio[i]->Divide(hChargedPtPPOnly[iRef]);
	}
	hChargedPtpp_ALICEpub->SetMarkerStyle(30);
	hChargedPtpp_ALICEpub->Draw("p,same");
	leg->AddEntry(hChargedPtpp_ALICEpub,"ALICE","pl");
	gr_alice_pp_2760GeV->SetMarkerStyle(31);
	gr_alice_pp_2760GeV->Draw("p,same");
	leg->AddEntry(gr_alice_pp_2760GeV,gr_alice_pp_2760GeV->GetTitle(),"pl");

	leg->Draw();
	// Calculation Various Ratios


	//==== Lower pad
	p = fpad->GetPad(2);
	p->SetTickx(); p->SetGridy(1); p->SetLogx(1), p->SetLogy(0); p->cd();
	TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, lowRAA, highRAA);
	hset( *hfr1, "p_{T} [GeV/c]", "Ratio",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
	hfr1->Draw();
	for(int i=0;i<NPP;i++) {
		hChargedPtPPOnly_Ratio[i]->Draw("p,same");
	}
	//gPad->GetCanvas()->SaveAs("figs_iaa/ppOnlyspectra_comp.pdf");

}


