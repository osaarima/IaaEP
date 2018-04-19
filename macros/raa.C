#include "include/Filipad.h"

double lowx=-1.;
double highx=50.;
double ly = 1e-7;
double hy = 2.2e5;
double lowRAA = 0.0;
double highRAA = 1.3;

TLatex latexRun;
TString strRunPP = "pp #sqrt{#it{s}} = 2.76 TeV";
TString strRunAA = "PbPb #sqrt{#it{s_{NN}}} = 2.76 TeV";

int gMarkers[] = {20,24,21,25,23,27,29,30};
int gColors[]={kBlack, kRed, kBlue, kDeepSea, kPink, kGray, kRed, kBlack};
void LoadALICEPub();
void CalculateRaa();
void DrawRaa();
void DrawPP();

const int NC_alice = 10; 
const int cent_alice[NC_alice+1] = {0,5,10,20,30,40,50,60,70,80,90};
const double Ncoll_alice[NC_alice] = {1686.87,1319.89,923.26,558.68,321.20,171.67,85.13,38.51,15.78,6.32}; //CentBinBorders          0 5 10 20 30 40 50 60 70 80 90 #
const double eNcoll_alice[NC_alice] = {197.7,153.7,99.6,56.4,31.0,15.2,8.0,3.8,1.3,0.5}; //CentBinBorders          0 5 10 20 30 40 50 60 70 80 90
const double  Npart_alice[NC_alice] = {382.8,329.7,260.5,186.4,128.9,85.0,52.8,30.0,15.8,7.52};
const double eNpart_alice[NC_alice] = { 3.1,4.6,4.4,3.9,3.3,2.6,2.0,1.3,0.6,0.4};

const int NC = 6;
const int cent[NC+1] =    {0,      5,    10,     20,    40,    60,       80};
const double  Ncoll[NC] = { 1686.87,1319.89, 923.26, 438.8, 128.2,  19.995 }; //60-90%
//const double  Ncoll[NC] = { 1502.7, 923.26, 438.8, 128.2,  26.82 }; // 60-80%

enum dataType { AA, pp };
TString inFile = "legotrain_JCIaa/data/JCIaa_legotrain_CF_PbPb_MC-964_20180407-1152-AMPT_LHC13f3c.root";
TString ppInFile = "legotrain_JCIaa/mc/JCIaaGF_pythia8230_pp2.76TeV_GF0-configSoftQCD.root";
TFile *fin[2];
TVector *CentBinBorders;
TH1D *hChargedPtpp_ALICEpub;
TH1D *hChargedPtAA[NC];
TH1D *hChargedPtpp;
TH1D *hChargedRAA[NC];

void Run(){
	CalculateRaa();
	LoadALICEPub();
	DrawRaa();
	//DrawPP();
}
void CalculateRaa() {
	cout <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	cout <<"PbPb : "<< inFile << endl;
	cout <<"pp :"<< ppInFile << endl;
	TString dirAA="JCIAA_EPInclusive";
	TString dirPP="JCIaa";

	fin[AA] = TFile::Open(inFile);
	fin[pp] = TFile::Open(ppInFile);

	TString TopDir[2] = {dirAA,dirPP};
	cout << TopDir[0] << endl;
	cout << TopDir[1] << endl;

	TH1D *hiCentr;
	TH1D *hiCentrpp;
	double noEvents[NC];
	double noEventspp;

	
	CentBinBorders              = (TVector*) fin[AA]->Get(dirAA+"/JCard/CentBinBorders");
	CentBinBorders->Print();

	double sc = 2*TMath::Pi()*1.6;

	int NumCent = CentBinBorders->GetNoElements()-1;
	hiCentr = (TH1D *)fin[AA]->Get(Form("%s/AliJHistos/hiCentr",dirAA.Data()));
	for(int ic=0;ic<NC;ic++) {
		hChargedPtAA[ic] = (TH1D *)fin[AA]->Get(Form("%s/AliJHistos/hChargedPt/hChargedPtC%02d",dirAA.Data(),ic));	
		cout << ic << endl;
		noEvents[ic] = hiCentr->GetBinContent(ic+1); // with vertex cut
		hChargedPtAA[ic]->Scale(1./noEvents[ic]/sc,"width");
	}
	hChargedPtpp = (TH1D *)fin[pp]->Get(Form("%s/AliJHistos/hChargedPt/hChargedPtC%02d",dirPP.Data(),0)); 
	hiCentrpp = (TH1D *)fin[pp]->Get(Form("%s/AliJHistos/hiCentr",dirPP.Data()));
	noEventspp = hiCentrpp->GetBinContent(1); // with vertex cut
	hChargedPtpp->Scale(1./noEventspp/sc,"width");

	// Calculate Raa
	for(int ic=0;ic<NC;ic++) {
		hChargedRAA[ic] = (TH1D*)hChargedPtAA[ic]->Clone();
		hChargedRAA[ic]->Divide(hChargedPtpp);
		hChargedRAA[ic]->Scale(1./Ncoll[ic]);
	}
}
void LoadALICEPub() {
	TString infileALICE = "jacek_Feb_14_2013/2012-11-16_dNdPt_all/dNdPt_pp_2.76TeV.root";
	double sigmaData = 62.2;
	TFile *finA = TFile::Open(infileALICE);
	hChargedPtpp_ALICEpub = (TH1D *)finA->Get("dNdPt_stat");
	hChargedPtpp_ALICEpub->Scale(sigmaData);
}


void DrawRaa() {
	// 1/N_{eve} 1/(2#pip_{T}) dN/dp_{T} [ (GeV/c)^{-2} ]

	Filipad *fpad;
	fpad = new Filipad(1, 1.1, 0.5, 100, 100, 0.7,2);
	fpad->Draw();
	//==== Upper pad
	TPad *p = fpad->GetPad(1); //upper pad
	p->SetTickx(); p->SetLogx(0); p->SetLogy(1); p->cd();
	hy = hChargedPtAA[0]->GetMaximum()*1.2;
	TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
	hset( *hfr, "p_{T} [GeV/c]", "1/N_{eve}1/(2#pip_{T})dN/dp_{T}[(GeV/c)^{-2}]",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
	hfr->Draw();
	//Legend definition
	TLegend *leg = new TLegend(0.45,0.4,0.85,0.78,"AMPT String Melting","brNDC");
	leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

	for(int i=0;i<NC;i++){
		hChargedPtAA[i]->SetMarkerStyle(gMarkers[i]);
		hChargedPtAA[i]->SetMarkerColor(gColors[i]);
		hChargedPtAA[i]->SetLineColor(gColors[i]);
		hChargedPtAA[i]->Draw("p,same");
		leg->AddEntry(hChargedPtAA[i],Form("%2.0f-%2.0f%%",(*CentBinBorders)[i+1], (*CentBinBorders)[i+2] ),"pl");
	}

	leg->Draw();
	
	//==== Lower pad
	p = fpad->GetPad(2);
	p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
	TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, lowRAA, highRAA);
	hset( *hfr1, "p_{T} [GeV/c]", "R_{AA}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
	hfr1->Draw();
	for(int i=0;i<NC;i++){
		hChargedRAA[i]->SetMarkerStyle(gMarkers[i]);
		hChargedRAA[i]->SetMarkerColor(gColors[i]);
		hChargedRAA[i]->SetLineColor(gColors[i]);			
		hChargedRAA[i]->Draw("p,same");
	}
	gPad->GetCanvas()->SaveAs("figs_iaa/Raa_AMPT13c.pdf");
}

void DrawPP() {
	// pp comparison
	Filipad *fpad;
	fpad = new Filipad(1, 1.1, 0.5, 100, 100, 0.7,2);
	fpad->Draw();
	//==== Upper pad
	TPad *p = fpad->GetPad(1); //upper pad
	p->SetTickx(); p->SetLogx(0); p->SetLogy(1); p->cd();
	hy = hChargedPtpp->GetMaximum()*1.2;
	TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, 1e-7, 3e3); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
	hset( *hfr, "p_{T} [GeV/c]", "1/N_{eve}1/(2#pip_{T})dN/dp_{T}[(GeV/c)^{-2}]",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
	//Ed^{3}#sigma/dp^{3} [mb GeV^{-2}c^{3}] ??
	hfr->Draw();
	//Legend definition
	TLegend *leg = new TLegend(0.45,0.4,0.85,0.78,"","brNDC");
	leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

	hChargedPtpp->SetMarkerStyle(20);
	hChargedPtpp->Draw("p,same");
	hChargedPtpp_ALICEpub->SetMarkerStyle(24);
	hChargedPtpp_ALICEpub->Draw("p,same");

	leg->AddEntry(hChargedPtpp,"Pythia SoftQCD","pl");
	leg->AddEntry(hChargedPtpp_ALICEpub,"ALICE","pl");

	leg->Draw();
	// Calculation Various Ratios


	//==== Lower pad
	p = fpad->GetPad(2);
	p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
	TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, lowRAA, highRAA);
	hset( *hfr1, "p_{T} [GeV/c]", "Ratio",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
	hfr1->Draw();
	gPad->GetCanvas()->SaveAs("figs_iaa/ppspectra_comp.pdf");

}


