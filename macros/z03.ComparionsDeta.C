
#include "include/Filipad.h"

void LoadData();
void compare();
void DrawSignal(int iPTT, int iPTA);

double lowx=-0.8;
double highx=0.8;
double ly = -0.1;
double hy = 0.4;
double lowIAA = 0.;
double highIAA = 2.;

TLatex latexRun;
TString strRun = "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV";

const int Nsets = 6;
TString infiles[Nsets] = {
	"sysErrors/Signal_AMPT_LHC13f3cJCIAA_TPC_E00_LHC11a_p4_AOD113_noSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3cJCIAA_TPC_E90_LHC11a_p4_AOD113_noSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3cJCIAA_V0A_E00_LHC11a_p4_AOD113_noSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3cJCIAA_V0A_E90_LHC11a_p4_AOD113_noSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3cJCIAA_V0P_E00_LHC11a_p4_AOD113_noSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3cJCIAA_V0P_E90_LHC11a_p4_AOD113_noSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
};
TFile *fin[Nsets];

int gMarkers[Nsets] = {20,24,21,25,22,26};
int gColors[Nsets]={kBlack, kRed, kBlue, kBlue, kPink, kGray};

const int kMAXD       = 20; //maximal number of pT trigger bins
const int kCENT       = 10; //maximal number of pT trigger bins
enum dataType { AA, pp };


TH1D *hDeltaEtaSig[Nsets][2][kCENT][kMAXD][kMAXD]; // background substracted signal based on fit
TH1D *hIAADeltaEtaSig[Nsets][kCENT][kMAXD][kMAXD]; // background substracted signal IAA

TVector *TriggPtBorders;
TVector *AssocPtBorders;
TVector *CentBinBorders;
int NumCent[2];
int NPTT;
int NPTA;
int iRef=0;

void LoadData() {
	
	for(int i=0;i<Nsets;i++){
		fin[i] = TFile::Open(infiles[i]);
	}

	int irefD = 0;
	TriggPtBorders             = (TVector*) fin[irefD]->Get("TriggPtBorders");
	AssocPtBorders             = (TVector*) fin[irefD]->Get("AssocPtBorders");
	CentBinBorders             = (TVector*) fin[irefD]->Get("CentBinBorders");
	NumCent[AA]    = CentBinBorders->GetNoElements()-1;
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
				} // ipta
			} // iptt 
		} // ic
	} // pp or AA
	} // iset
}


void Compare(){
	LoadData();
	int iPTT=3;
	int iPTA=4;

	DrawSignal(iPTT, iPTA);

}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void DrawSignal(int iPTT, int iPTA) {
	lowx = -0.01;
	for(int ic=0;ic<NumCent[AA];ic++) {
		Filipad *fpad = new Filipad(ic+1, 1.1, 0.5, 100, 100, 0.7,NumCent[AA]);
		fpad->Draw();
		//==== Upper pad
		TPad *p = fpad->GetPad(1); //upper pad
		p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
		TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
		hset( *hfr, "|#Delta#eta|", "1/N_{trigg} dN/d#Delta#eta",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
		hfr->Draw();
		//Legend definition
		TLegend *leg = new TLegend(0.45,0.6,0.85,0.82,"","brNDC");
		leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

		latexRun.DrawLatexNDC( 0.47, 0.85 ,strRun);

		for(int iS=0;iS<Nsets;iS++) {
			hDeltaEtaSig[iS][AA][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[iS]);
			hDeltaEtaSig[iS][AA][ic][iPTT][iPTA]->Draw("p,same");
		}


		leg->AddEntry(hDeltaEtaSig[0][AA][ic][iPTT][iPTA],hDeltaEtaSig[0][AA][ic][iPTT][iPTA]->GetTitle(),"p");
		
		leg->Draw();

		//==== Lower pad
		p = fpad->GetPad(2);
		p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
		TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, lowIAA, highIAA);
		hset( *hfr1, "|#Delta#eta|", "AA/pp",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
		hfr1->Draw();
		for(int iS=0;iS<Nsets;iS++) {
			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[iS]);
			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->Draw("p,same");
		}
		//gPad->GetCanvas()->SaveAs(Form("figs_svn/FigA4_v%d_modelcomparisonBest.eps",i+2));
	}
}
