
#include "include/Filipad.h"

void LoadData();
void DrawPP();

double lowx=-20.;
double highx=20.;
double ly = 1e-5;
double hy = 100.;
double lowIAA = -0.2;
double highIAA = 5;

/*
void hset(TH1& hid, TString xtit="", TString ytit="",             
        double titoffx = 0.9, double titoffy = 1.2,               
        double titsizex = 0.06, double titsizey = 0.06,           
        double labeloffx = 0.01, double labeloffy = 0.001,        
        double labelsizex = 0.05, double labelsizey = 0.05,       
        int divx = 510, int divy=510)                             
{                                                                 
    hid.SetStats(0);                                              
                                                                  
    hid.GetXaxis()->CenterTitle(1);                               
    hid.GetYaxis()->CenterTitle(1);                               
                                                                  
    hid.GetXaxis()->SetTitleOffset(titoffx);                      
    hid.GetYaxis()->SetTitleOffset(titoffy);                      
                                                                  
    hid.GetXaxis()->SetTitleSize(titsizex);                       
    hid.GetYaxis()->SetTitleSize(titsizey);                       
                                                                  
    hid.GetXaxis()->SetLabelOffset(labeloffx);                    
    hid.GetYaxis()->SetLabelOffset(labeloffy);                    
                                                                  
    hid.GetXaxis()->SetLabelSize(labelsizex);                     
    hid.GetYaxis()->SetLabelSize(labelsizey);                     
                                                                  
    hid.GetXaxis()->SetNdivisions(divx);                          
    hid.GetYaxis()->SetNdivisions(divy);                          
                                                                  
    hid.GetXaxis()->SetTitle(xtit);                               
    hid.GetYaxis()->SetTitle(ytit);                               
}                                                                 
*/


TLatex latexRun;
TString strRun = "PYTHIA #sqrt{#it{s}} = 2.76 TeV";

const int Nsets = 6;
TString infiles[Nsets] = {
	"vtx/zvtx.root",
	"vtx/zvtx.root",
	"vtx/zvtx.root",
	"vtx/zvtx.root",
	"vtx/zvtx.root",
	"vtx/zvtx.root"
};
TFile *fin[Nsets];
TString hnames[Nsets] = {
	"LHC15o_Pb_Pb_5TeV_z",
	"LHC17n_Xe_Xe_544TeV_z",
	"LHC16l_pp_13TeV_z",
	"LHC15n_pp_5TeV_z",
	"LHC11a_pp_2760GeV_z",
	"LHC10c_pp_900GeV_z"
};

TString sLeg[Nsets] = {
	"#sqrt{s_{NN}} = 5.02 TeV, PbPb",
	"#sqrt{s_{NN}} = 5.44 TeV, XeXe",
	"#sqrt{s} = 13 TeV, pp",
	"#sqrt{s} = 5 TeV, pp",
	"#sqrt{s} = 2.76 TeV, pp",
	"#sqrt{s} = 0.9 TeV, pp"
};

int gMarkers[] = {20,24,21,25,23,27,29,30};
int gColors[]={kBlack, kRed, kBlue, kDeepSea, kPink, kRed, kBlack};

enum dataType { AA, pp };
TH1D *hvtx[Nsets];
TH1D *hRatios[Nsets];

int iRef=1; // data:

void LoadData();
void DrawVtx();
void vtx() {
	LoadData();
	DrawVtx();
}

//------------------------------------------------------------------------------------------------
void LoadData() {

	for(int i=0;i<Nsets;i++){
		fin[i] = TFile::Open(infiles[i]);
	}

	int irefD = 0;
	for(int i=0;i<Nsets;i++){
		hvtx[i] = (TH1D *)fin[i]->Get(Form("%s",hnames[i].Data()));
		//hvtx[i]->Scale(1./hvtx[i]->GetEntries());
	} // iset
}

//------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------
void DrawVtx() {
	Filipad *fpad;
	int padID = 1;
	fpad = new Filipad(padID+1, 1.1, 0.5, 100, 100, 0.7,1);
	fpad->Draw();
	//==== Upper pad
	TPad *p = fpad->GetPad(1); //upper pad
	p->SetTickx(); p->SetLogx(0); p->SetLogy(1); p->cd();
	TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
	hset( *hfr, "z_vtx", "A.U",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
	hfr->Draw();
	//Legend definition
	TLegend *leg = new TLegend(0.65,0.5,0.85,0.78,"","brNDC");
	leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

	for(int i=0;i<Nsets;i++){
		hvtx[i]->SetMarkerStyle(gMarkers[i]);
		hvtx[i]->SetMarkerColor(gColors[i]);
		hvtx[i]->SetLineColor(gColors[i]);
		hvtx[i]->Draw("p,same");
		leg->AddEntry(hvtx[i],Form("%s",sLeg[i].Data()),"pl");
	}

	leg->Draw();
	for(int i=0;i<Nsets;i++){
		hRatios[i] = (TH1D*)hvtx[i]->Clone();
		hRatios[i]->Divide(hvtx[iRef]);
	}
	//==== Lower pad
	p = fpad->GetPad(2);
	p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
	TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, lowIAA, highIAA);
	hset( *hfr1, "z_vtx", "Ratio",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
	hfr1->Draw();
	for(int i=0;i<Nsets;i++){
		hRatios[i]->Draw("p,same");
	}
	gPad->GetCanvas()->SaveAs(Form("figs_iaa/vtx.pdf"));
}
