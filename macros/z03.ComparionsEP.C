
#include "include/Filipad.h"

void LoadData();
void compare();
void DrawSignal(int padid, int iPTT, int iPTA);
void DrawIAA(int padID, int iPTT, int iPTA);

double lowx=-0.1;
double highx=0.6;
double ly = -0.05;
double hy = 5.6;
double lowIAA = -0.2;
double highIAA = 2.2;

TLatex latexRun;
TString strRun = "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV";

const int Nsets = 4;
TString infiles[Nsets] = {
	"sysErrors/Signal_LHC10h_AOD86_MgFpMgFm_5217_JCIAA_TPCOnly_H0_T0_LHC11a_p4_AOD113_noSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	//"sysErrors/Signal_LHC10h_AOD86_MgFpMgFm_5317_JCIAA_V0C_E90_LHC11a_p4_AOD113_noSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root"};
	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_LHC12f1a_Pythia_2760GeV_KineOnly_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_V0C_E00_LHC12f1a_Pythia_2760GeV_KineOnly_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_V0C_E90_LHC12f1a_Pythia_2760GeV_KineOnly_Iaa_R0.2_1.0_1.60_Near_Wing0.root"};
	/*
	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_pythia8230_pp2.76TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_TPC_E00_pythia8230_pp2.76TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_TPC_E90_pythia8230_pp2.76TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_V0A_E00_pythia8230_pp2.76TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_V0A_E90_pythia8230_pp2.76TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_V0P_E00_pythia8230_pp2.76TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_V0P_E90_pythia8230_pp2.76TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
};
*/
TFile *fin[Nsets];

TString sLeg[Nsets] = {
	"lhc10h Inclusive",
	//"lhc10h Out"};
	"AMPT V0C Inclusive",
	"AMPT V0C In",
	"AMPT V0C Out"};
/*
	"AMPT String melting",
	"AMPT TPC In",
	"AMPT TPC Out",
	"AMPT V0A In",
	"AMPT V0A Out",
	"AMPT V0P In",
	"AMPT V0P Out"
};
*/
int gMarkers[] = {20,24,21,25,23,27,29,30};
int gColors[]={kBlack, kRed, kBlue, kDeepSea, kPink, kGray, kRed, kBlack};

const int kMAXD       = 20; //maximal number of pT trigger bins
const int kCENT       = 10; //maximal number of pT trigger bins
enum dataType { AA, pp };


TH1D *hDeltaEtaSig[Nsets][2][kCENT][kMAXD][kMAXD]; // background substracted signal based on fit
TH1D *hIAADeltaEtaSig[Nsets][kCENT][kMAXD][kMAXD]; // background substracted signal IAA

TH1D *hRatioEP_DeltaEtaSig[Nsets][kCENT][kMAXD][kMAXD]; //Data TPC V0A V0P to AMPT inclusive
TH1D *hRatioEP_IAA[Nsets][kCENT][kMAXD][kMAXD]; //Data TPC V0A V0P to AMPT inclusive

TVector *TriggPtBorders;
TVector *AssocPtBorders;
TVector *CentBinBorders;
int NumCent[2];
int NPTT;
int NPTA;
int iRef=0;

//------------------------------------------------------------------------------------------------
void LoadData() {
	
	for(int i=0;i<Nsets;i++){
		fin[i] = TFile::Open(infiles[i]);
	}

	int irefD = 0;
	TriggPtBorders             = (TVector*) fin[irefD]->Get("TriggPtBorders");
	AssocPtBorders             = (TVector*) fin[irefD]->Get("AssocPtBorders");
	CentBinBorders             = (TVector*) fin[irefD]->Get("CentBinBorders");
	NumCent[AA]    = CentBinBorders->GetNoElements()-2; // for 5TeV
	NumCent[pp]    = 1; 
	NPTT     = TriggPtBorders->GetNoElements()-1;
	NPTA     = AssocPtBorders->GetNoElements()-1;
	cout <<"PbPb"<<endl;
	cout <<"bins:  c="<<  NumCent[AA] <<" ptt="<< NPTT <<" pta="<< NPTA  << endl; 
	cout <<"+++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	//------------ R e a d    D a t a ------------    
	cout <<"Reading data...."<<endl;
	for(int i=0;i<Nsets;i++){
		for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp  //(*fCentralityBinBorders)[i+1]
			for(int ic=0; ic<NumCent[idtyp]; ic++){
				for(int iptt=0; iptt<NPTT; iptt++){
					for(int ipta=0;ipta<NPTA;ipta++) {
						hDeltaEtaSig[i][idtyp][ic][iptt][ipta] = (TH1D *)fin[i]->Get(Form("hDeltaEtaSig%02dC%02dT%02dA%02d",idtyp,ic,iptt,ipta));
						if(idtyp==AA) hIAADeltaEtaSig[i][ic][iptt][ipta] = (TH1D *)fin[i]->Get(Form("hIAADeltaEtaSigC%02dT%02dA%02d",ic,iptt,ipta));
						printf("Set%d C%d T%d A %d\n",i, ic, iptt, ipta);
						hDeltaEtaSig[i][idtyp][ic][iptt][ipta]->Print();
					} // ipta
				} // iptt 
			} // ic
		} // pp or AA
	} // iset
	// Calculation Various Ratios
	// 1. Out/In ratio for each data setfor(int i=0;i<Nsets;i++){
	cout <<"Calculating ratios...."<<endl;
	for(int ic=0; ic<NumCent[AA]; ic++){
		for(int iptt=0; iptt<NPTT; iptt++){
			for(int ipta=0;ipta<NPTA;ipta++) {
				for(int i=0;i<Nsets;i++){
					hRatioEP_DeltaEtaSig[i][ic][iptt][ipta] = (TH1D*)hDeltaEtaSig[i][AA][ic][iptt][ipta]->Clone();
					hRatioEP_DeltaEtaSig[i][ic][iptt][ipta]->Divide(hDeltaEtaSig[iRef][AA][ic][iptt][ipta]);

					hRatioEP_IAA[i][ic][iptt][ipta] = (TH1D*)hIAADeltaEtaSig[i][ic][iptt][ipta]->Clone();
					hRatioEP_IAA[i][ic][iptt][ipta]->Divide(hIAADeltaEtaSig[iRef][ic][iptt][ipta]);
				}
			} // ipta
		} // iptt 
	} // ic
}

//------------------------------------------------------------------------------------------------
void Compare(){
	LoadData();

	int ic=0;
	for(int iptt=3; iptt<NPTT; iptt++){
		for(int ipta=2;ipta<NPTA;ipta++) {
			DrawSignal(ic++,iptt, ipta);
			//DrawIAA(ic++,iptt, ipta);
	 	}
	}

}




//------------------------------------------------------------------------------------------------
void DrawSignal(int padID, int iPTT, int iPTA) {
	Filipad *fpad[NumCent[AA]];
	lowx = -0.01;
	for(int ic=0;ic<NumCent[AA];ic++){
		fpad[ic] = new Filipad(padID+ic+1, 1.1, 0.5, 100, 100, 0.7,NumCent[AA]);
		fpad[ic]->Draw();
		//==== Upper pad
		TPad *p = fpad[ic]->GetPad(1); //upper pad
		p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
		hy = hDeltaEtaSig[iRef][AA][ic][iPTT][iPTA]->GetMaximum()*1.2;
		TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
		hset( *hfr, "|#Delta#eta|", "1/N_{trigg} dN/d|#Delta#eta|",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
		hfr->Draw();
		//Legend definition
		TLegend *leg = new TLegend(0.45,0.4,0.85,0.78,"","brNDC");
		leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

		latexRun.DrawLatexNDC( 0.25, 0.85 ,strRun);

		leg->AddEntry((TObject*)NULL,hDeltaEtaSig[0][AA][ic][iPTT][iPTA]->GetTitle(),"");

		for(int iS=0;iS<Nsets;iS++) {
			hDeltaEtaSig[iS][AA][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[iS]);
			hDeltaEtaSig[iS][AA][ic][iPTT][iPTA]->SetMarkerColor(gColors[iS]);
			hDeltaEtaSig[iS][AA][ic][iPTT][iPTA]->SetLineColor(gColors[iS]);
			hDeltaEtaSig[iS][AA][ic][iPTT][iPTA]->Draw("p,same");
			leg->AddEntry(hDeltaEtaSig[iS][AA][ic][iPTT][iPTA],sLeg[iS],"pl");
		}

		
		leg->Draw();

		//==== Lower pad
		p = fpad[ic]->GetPad(2);
		p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
		TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, lowIAA, highIAA);
		hset( *hfr1, "|#Delta#eta|", Form("Ratio to %s",sLeg[iRef].Data()),1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
		hfr1->Draw();
		for(int i=0;i<Nsets;i++) {
			hRatioEP_DeltaEtaSig[i][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[i]);
			hRatioEP_DeltaEtaSig[i][ic][iPTT][iPTA]->SetMarkerColor(gColors[i]);
			hRatioEP_DeltaEtaSig[i][ic][iPTT][iPTA]->SetLineColor(gColors[i]);
			hRatioEP_DeltaEtaSig[i][ic][iPTT][iPTA]->Draw("p,same");
		}
		gPad->GetCanvas()->SaveAs(Form("figs_iaa/PbPbCompEP_C%02dT%02dA%02d.pdf",ic,iPTT,iPTA));
	}
}

//------------------------------------------------------------------------------------------------
void DrawIAA(int padID, int iPTT, int iPTA) {
	Filipad *fpad[NumCent[AA]];
	lowx = -0.01;
	for(int ic=0;ic<NumCent[AA];ic++){
		fpad[ic] = new Filipad(padID+ic+1, 1.1, 0.5, 100, 100, 0.7,NumCent[AA]);
		fpad[ic]->Draw();
		//==== Upper pad
		TPad *p = fpad[ic]->GetPad(1); //upper pad
		p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
		TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, lowIAA, highIAA); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
		hset( *hfr, "|#Delta#eta|", "IAA",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
		hfr->Draw();
		//Legend definition
		TLegend *leg = new TLegend(0.45,0.4,0.85,0.78,"","brNDC");
		leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

		latexRun.DrawLatexNDC( 0.25, 0.85 ,strRun);

		leg->AddEntry((TObject*)NULL,hIAADeltaEtaSig[0][ic][iPTT][iPTA]->GetTitle(),"");

		for(int iS=0;iS<Nsets;iS++) {
			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[iS]);
			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->SetMarkerColor(gColors[iS]);
			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->SetLineColor(gColors[iS]);
			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->Draw("p,same");
			leg->AddEntry(hIAADeltaEtaSig[iS][ic][iPTT][iPTA],sLeg[iS],"pl");
		}

		
		leg->Draw();

		//==== Lower pad
		p = fpad[ic]->GetPad(2);
		p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
		TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, -0.1, 5);
		hset( *hfr1, "|#Delta#eta|", Form("Ratio to %s",sLeg[iRef].Data()),1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
		hfr1->Draw();
		for(int i=0;i<Nsets;i++) {
			hRatioEP_IAA[i][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[i]);
			hRatioEP_IAA[i][ic][iPTT][iPTA]->SetMarkerColor(gColors[i]);
			hRatioEP_IAA[i][ic][iPTT][iPTA]->SetLineColor(gColors[i]);
			hRatioEP_IAA[i][ic][iPTT][iPTA]->Draw("p,same");
		}
		gPad->GetCanvas()->SaveAs(Form("figs_iaa/PbPbCompEPIAA_C%02dT%02dA%02d.pdf",ic,iPTT,iPTA));
	}
}



