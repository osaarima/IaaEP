
#include "include/Filipad.h"

void LoadData();
void LoadDataQM();
void compare();
void DrawSignal(int padid, int iPTT, int iPTA);
void DrawIAA(int padID, int iPTT, int iPTA);
void DrawIAAQM(int padID, int iPTT, int iPTA);

double lowx=-0.1;
double highx=0.3;
double ly = -0.05;
double hy = 4.0;
double lowIAA = -0.2;
double highIAA = 5.2;

TLatex latexRun;
TString strRun = "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV";

const int Nsets = 2;
TString infiles[Nsets] = {
//	"sysErrors/Signal_LHC10h_AOD86_MgFpMgFm_5217_JCIAA_TPCOnly_H0_T0_LHC11a_p4_AOD113_noSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
  "sysErrors/IAA_QM12.root",
  "sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_LHC11a_AOD113_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
//	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_pythia8230_pp2.76TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
  
};
TFile *fin[Nsets];

TString sLeg[Nsets] = {
	"LHC10h, TPCOnly",
	"AMPT String melting"
};

int gMarkers[] = {20,24,21,25,23,27,29,30};
int gColors[]={kBlack, kRed, kBlue, kDeepSea, kPink, kGray, kRed, kBlack};

const int kMAXD       = 20; //maximal number of pT trigger bins
const int kCENT       = 10; //maximal number of pT trigger bins
enum dataType { AA, pp };


TH1D *hDeltaEtaSig[Nsets][2][kCENT][kMAXD][kMAXD]; // background substracted signal based on fit
TH1D *hIAADeltaEtaSig[Nsets][kCENT][kMAXD][kMAXD]; // background substracted signal IAA

TH1D *hRatio_DeltaEtaSig[Nsets][kCENT][kMAXD][kMAXD]; //Data TPC V0A V0P to AMPT inclusive
TH1D *hRatio_IAA[Nsets][kCENT][kMAXD][kMAXD];

TGraphAsymmErrors *grFinalIAAStat[Nsets][kCENT][kMAXD][kMAXD];      //For Filip's QM, IAA with statstistical 
TGraphAsymmErrors *grFinalIAASyst[Nsets][kCENT][kMAXD][kMAXD];      //For Filip's QM, IAA with Systematics
TGraphAsymmErrors *grFinalIAAScale[Nsets][kCENT][kMAXD][kMAXD];      //For Filip's QM, Converting box to TGraph
TBox *sysbox[Nsets][kCENT][kMAXD][kMAXD];                           //For Filip's QM, Scaling Box


TVector *TriggPtBorders;
TVector *AssocPtBorders;
TVector *CentBinBorders;
int NumCent[2];
int NPTT;
int NPTA;
int iRef=1; // data -0

//------------------------------------------------------------------------------------------------
void LoadData() {
	
	for(int i=0;i<Nsets;i++){
		fin[i] = TFile::Open(infiles[i]);
	}

	int irefD = 0;
	TriggPtBorders             = (TVector*) fin[irefD]->Get("TriggPtBorders");
	AssocPtBorders             = (TVector*) fin[irefD]->Get("AssocPtBorders");
	CentBinBorders             = (TVector*) fin[irefD]->Get("CentBinBorders");
	NumCent[AA]    = CentBinBorders->GetNoElements()-2;  // 5TeV 1 less cent
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
	// Calculation Various Ratios
	// 1. Out/In ratio for each data setfor(int i=0;i<Nsets;i++){
	for(int ic=0; ic<NumCent[AA]; ic++){
		for(int iptt=0; iptt<NPTT; iptt++){
			for(int ipta=0;ipta<NPTA;ipta++) {
				for(int i=0;i<Nsets;i++){
					hRatio_DeltaEtaSig[i][ic][iptt][ipta] = (TH1D*)hDeltaEtaSig[i][AA][ic][iptt][ipta]->Clone();
					hRatio_DeltaEtaSig[i][ic][iptt][ipta]->Divide(hDeltaEtaSig[iRef][AA][ic][iptt][ipta]);

					hRatio_IAA[i][ic][iptt][ipta] = (TH1D*)hIAADeltaEtaSig[i][ic][iptt][ipta]->Clone();
					hRatio_IAA[i][ic][iptt][ipta]->Divide(hIAADeltaEtaSig[iRef][ic][iptt][ipta]);
				}
			} // ipta
		} // iptt 
	} // ic
}

void LoadDataQM() {
	
	for(int i=0;i<Nsets;i++){
		fin[i] = TFile::Open(infiles[i]);
	}

	int irefD = 1;
	TriggPtBorders             = (TVector*) fin[irefD]->Get("TriggPtBorders");
	AssocPtBorders             = (TVector*) fin[irefD]->Get("AssocPtBorders");
	CentBinBorders             = (TVector*) fin[irefD]->Get("CentBinBorders");
	NumCent[AA]    = CentBinBorders->GetNoElements()-2;  // 5TeV 1 less cent
	NumCent[pp]    = 1; 
	NPTT     = TriggPtBorders->GetNoElements()-1;
	NPTA     = AssocPtBorders->GetNoElements()-1;
	cout <<"PbPb"<<endl;
	cout <<"bins:  c="<<  NumCent[AA] <<" ptt="<< NPTT <<" pta="<< NPTA  << endl; 
	cout <<"+++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	//------------ R e a d    D a t a ------------    
	cout <<"Reading data...."<<endl;
	for(int i=0;i<Nsets;i++){
    if (i != 0) {
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
    }  // Read ONLY Our Analysis Files
	} // iset
	// Calculation Various Ratios
	// 1. Out/In ratio for each data setfor(int i=0;i<Nsets;i++){
	for(int ic=0; ic<NumCent[AA]; ic++){
		for(int iptt=0; iptt<NPTT; iptt++){
			for(int ipta=0;ipta<NPTA;ipta++) {
				for(int i=0;i<Nsets;i++){
          if (i != 0) {
            hRatio_DeltaEtaSig[i][ic][iptt][ipta] = (TH1D*)hDeltaEtaSig[i][AA][ic][iptt][ipta]->Clone();
            hRatio_DeltaEtaSig[i][ic][iptt][ipta]->Divide(hDeltaEtaSig[iRef][AA][ic][iptt][ipta]);

            hRatio_IAA[i][ic][iptt][ipta] = (TH1D*)hIAADeltaEtaSig[i][ic][iptt][ipta]->Clone();
            hRatio_IAA[i][ic][iptt][ipta]->Divide(hIAADeltaEtaSig[iRef][ic][iptt][ipta]);
//            cout << hRatio_IAA[i][ic][iptt][ipta]
          }  // Do it for only Our Analysis Files
				}
			} // ipta
		} // iptt 
	} // ic

  int pttmap[] = {0, 0, 1, 2};        // ptt = 3 -> filip 2 ...
  int ptamap[] = {0, 0, 0, 1, 2, 3, 4};  // pta = 6 -> filip 3 pta = 5 -> filip 3
  int centmap[] = {0, 0, 1, 2, 3, 4};

  double fPttArray[] = {3, 4, 6, 8, 15};
  double fPtaArray[] = {0.5, 1, 2, 3, 4, 6, 8, 10};

  for (int iset = 0; iset < 1; ++iset) {
    for(int ic=0; ic<NumCent[AA]; ic++){
      for(int iptt=0; iptt<NPTT; iptt++){
        for(int ipta=0;ipta<NPTA -1;ipta++) {
          
          if (fPttArray[ pttmap[iptt] ] - fPtaArray[ ptamap[ipta+1] ] < 0.1) continue;

//          if (ipta == 6) continue;
          cout << "Difference : " <<  ic << iptt << ipta << " " << fPttArray[ pttmap[iptt] ] - fPtaArray[ ptamap[ipta+1] ] <<  endl;
//          if (fPttArray[iptt] - fPtaArray[ipta] < -1.1) continue;

            grFinalIAASyst[iset][ic][iptt][ipta] = (TGraphAsymmErrors *)fin[iset]->Get(Form("grFinalIAASyst%02d%02d%02d", centmap[ic], pttmap[iptt], ptamap[ipta]));
            grFinalIAAStat[iset][ic][iptt][ipta] = (TGraphAsymmErrors *)fin[iset]->Get(Form("grFinalIAAStat%02d%02d%02d", centmap[ic], pttmap[iptt], ptamap[ipta]));
            sysbox[iset][ic][iptt][ipta] = (TBox *)fin[iset]->Get(Form("sysbox%02d%02d%02d", centmap[ic], pttmap[iptt], ptamap[ipta]));
            cout << iset << ic << iptt << ipta << " " << endl;
            cout << iset << centmap[ic] << pttmap[iptt] << ptamap[ipta] << " " << sysbox[iset][ic][iptt][ipta] << endl;

            double errorlow = 1 - sysbox[iset][ic][iptt][ipta]->GetY1();
            double errorhigh = 1 - sysbox[iset][ic][iptt][ipta]->GetY2();
            grFinalIAAScale[iset][ic][iptt][ipta] = new TGraphAsymmErrors();
            grFinalIAAScale[iset][ic][iptt][ipta] ->SetName(Form("grFinalIAAScale%02d%02d%02d", centmap[ic], pttmap[iptt], ptamap[ipta]));
            grFinalIAAScale[iset][ic][iptt][ipta]->SetPoint(0, 0.2, 1);
            grFinalIAAScale[iset][ic][iptt][ipta]->SetPointError(0, 1, 1, errorlow, errorhigh);
            grFinalIAAScale[iset][ic][iptt][ipta]->SetFillStyle(3003);
            grFinalIAAScale[iset][ic][iptt][ipta]->SetFillColorAlpha(kBlack, 0.6);

        }
      }
    }
  }


}





//------------------------------------------------------------------------------------------------
void Compare(){
//	LoadData();
  LoadDataQM();
	int ic=0;
	for(int iptt=3; iptt<NPTT; iptt++){
		for(int ipta=1;ipta<NPTA-1;ipta++) {
//			DrawIAA(ic++,iptt, ipta);
      DrawIAAQM(ic++, iptt, ipta);
	 	}
	}

}


//------------------------------------------------------------------------------------------------
void DrawIAA(int padID, int iPTT, int iPTA) {
	Filipad *fpad[NumCent[AA]];
	lowx = -0.01;
	for(int ic=0;ic<NumCent[AA];ic++) {
		fpad[ic] = new Filipad(padID+ic+1, 1.1, 0.5, 100, 100, 0.7,NumCent[AA]);
		fpad[ic]->Draw();
		//==== Upper pad
		TPad *p = fpad[ic]->GetPad(1); //upper pad
		p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
		TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
		hset( *hfr, "|#Delta#eta|", "IAA(|#Delta#eta|)",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
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
		TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, lowIAA, highIAA);
		hset( *hfr1, "|#Delta#eta|", Form("Ratio to %s",sLeg[iRef].Data()),1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
		hfr1->Draw();
		for(int i=0;i<Nsets;i++) {
			hRatio_IAA[i][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[i]);
			hRatio_IAA[i][ic][iPTT][iPTA]->SetMarkerColor(gColors[i]);
			hRatio_IAA[i][ic][iPTT][iPTA]->SetLineColor(gColors[i]);
			hRatio_IAA[i][ic][iPTT][iPTA]->Draw("p,same");
		}
		gPad->GetCanvas()->SaveAs(Form("figs_iaa/IAA_PbPbComp_C%02dT%02dA%02d.pdf",ic,iPTT,iPTA));
	}
	//for(int ic=0;ic<NumCent[AA];ic++) delete fpad[ic];
}

void DrawIAAQM(int padID, int iPTT, int iPTA) {
  
	TCanvas *fpad[NumCent[AA]];
	lowx = -0.01;
	for(int ic=0;ic<NumCent[AA];ic++) {
		fpad[ic] = new TCanvas(Form("fpad%d", padID+ic+1), Form(""), 800, 800);
    fpad[ic]->SetTopMargin(0.05);
    fpad[ic]->SetRightMargin(0.05);
    fpad[ic]->SetLeftMargin(0.15);
    fpad[ic]->SetBottomMargin(0.12);
		fpad[ic]->Draw();
		//==== Upper pad
//		TPad *p = fpad[ic]->GetPad(0); //upper pad
//		p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
    fpad[ic]->SetTickx(); fpad[ic]->SetLogx(0); fpad[ic]->SetLogy(0);
		TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
    hfr->SetStats(0);
    hfr->GetXaxis()->CenterTitle(1);
    hfr->GetYaxis()->CenterTitle(1);

    hfr->GetXaxis()->SetTitleOffset(1.0);
    hfr->GetYaxis()->SetTitleOffset(1.0);

//    hfr->GetXaxis()->SetTitleFont(43);
//    hfr->GetYaxis()->SetTitleFont(43);
    hfr->GetXaxis()->SetTitleSize(0.06);
    hfr->GetYaxis()->SetTitleSize(0.06);

    hfr->GetXaxis()->SetLabelOffset(0.01);
    hfr->GetYaxis()->SetLabelOffset(0.01);

//    hfr->GetXaxis()->SetLabelFont(43);
//    hfr->GetYaxis()->SetLabelFont(43);
    hfr->GetXaxis()->SetLabelSize(0.04);
    hfr->GetYaxis()->SetLabelSize(0.05);

//    hfr->GetXaxis()->SetNdivisions(510);
//    hfr->GetYaxis()->SetNdivisions(505);

    hfr->GetXaxis()->SetTitle("|#Delta#eta|");
    hfr->GetYaxis()->SetTitle("I_AA (|#Delta#eta|)");

		hfr->Draw();
		//Legend definition
		TLegend *leg = new TLegend(0.45,0.4,0.85,0.78,"","brNDC");
		leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

		latexRun.DrawLatexNDC( 0.25, 0.85 ,strRun);

		leg->AddEntry((TObject*)NULL,hIAADeltaEtaSig[1][ic][iPTT][iPTA]->GetTitle(),"");

		for(int iS=1;iS<Nsets;iS++) {
			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[iS]);
			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->SetMarkerColor(gColors[iS]);
			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->SetLineColor(gColors[iS]);
			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->Draw("p,same");
			leg->AddEntry(hIAADeltaEtaSig[iS][ic][iPTT][iPTA],sLeg[iS],"pl");
		}

    grFinalIAAStat[0][ic][iPTT][iPTA]->Draw("p");
    grFinalIAASyst[0][ic][iPTT][iPTA]->Draw("p");
    grFinalIAAScale[0][ic][iPTT][iPTA]->Draw("2");
//    sysbox[0][ic][iPTT][iPTA]->Draw();
		
		leg->Draw();

		fpad[ic]->SaveAs(Form("figs_iaa/HP_IAA_PbPb_C%02dT%02dA%02d.pdf",ic,iPTT,iPTA));
	}
	//for(int ic=0;ic<NumCent[AA];ic++) delete fpad[ic];
}


