
#include "include/Filipad.h"


TGraphErrors* get_ratio( TGraphErrors * l, TGraphErrors *r );

const int kMAXD       = 20; //maximal number of pT trigger bins
const int kCENT       = 10; //maximal number of pT trigger bins
const int kZvtx       = 15; //maximal number of pT trigger bins

TH1D *hTriggPtBin[2][kCENT][kMAXD]; 
TH1D *hAssocPtBin[2][kCENT][kMAXD][kMAXD]; 

TH1D *hIAAEta[kCENT][kMAXD][kMAXD]; // in eta,phi 

TFile *fin;
// adding DeltaEta histograms after mixed event corrections
// with $\Delta\phi$ < 0.2 ??? check
TH1D *hDeltaEta[2][kCENT][kMAXD][kMAXD]; // filipp summed DeltaEta AA-0 pp-1
// save this into an additional root file, fit and IAA will be calculated in z02.CalIAADeta.C
Bool_t saveDeta = kTRUE;
double lowx=-0.5;
double highx=0.5;
double ly = -0.1;
double hy = 1.5;



void DoAnalysis(TString inFile="sysErrors/_AA_moon1_pp_moon1_Iaa_R0.2_1.0_1.60_Near_Wing0.root",  TString oname=""){

	enum dataType { AA, pp };
	cout <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	cout <<"PbPb : "<< inFile << endl;
	cout <<"Discription :"<< oname << endl;

	fin = TFile::Open(inFile);
	
	TVector *TriggPtBorders;
	TVector *AssocPtBorders;
	TVector *CentBinBorders;

	TriggPtBorders             = (TVector*) fin->Get("TriggPtBorders");
	AssocPtBorders             = (TVector*) fin->Get("AssocPtBorders");
	CentBinBorders             = (TVector*) fin->Get("CentBinBorders");

	int reference = 0; //0=pp, 1=60-90%
	int logy = 1;
	int plotIAA = 1;
	int saveRoot = 1;  
	int sdx = 250;
	double lpta = 0.;
	double hpta = 10.;
	double hIAA = 5.;
	// for merging mixed event

	int NumCent[2]    = { CentBinBorders->GetNoElements()-1, 1}; 
	int NumPtt     = TriggPtBorders->GetNoElements()-1;
	int NumPta     = AssocPtBorders->GetNoElements()-1;
	cout <<"PbPb"<<endl;
	cout <<"bins:  c="<<  NumCent <<" ptt="<< NumPtt <<" pta="<< NumPta  << endl; 
	cout <<"+++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

	const int NC = NumCent[AA];
	const int NPTT = NumPtt;
	const int NPTA = NumPta;
	double MeanPta[NC][NPTT][NPTA];


	//------------ R e a d    D a t a ------------    
	cout <<"Reading data...."<<endl;
	for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp  //(*fCentralityBinBorders)[i+1]
		for(int ic=0; ic<NumCent[idtyp]; ic++){
			for(int iptt=0; iptt<NPTT; iptt++){
				for(int ipta=0;ipta<NPTA;ipta++) {
					hDeltaEta[idtyp][ic][iptt][ipta] = (TH1D *)fin->Get(Form("hDeltaEtaType%02dC%02dT%02dA%02d",idtyp,ic,iptt,ipta));
					hDeltaEta[idtyp][ic][iptt][ipta]->Print();
					hDeltaEta[idtyp][ic][iptt][ipta]->SetXTitle("#Delta#eta");
					hDeltaEta[idtyp][ic][iptt][ipta]->GetXaxis()->CenterTitle(kTRUE);
					hDeltaEta[idtyp][ic][iptt][ipta]->GetXaxis()->SetTitleOffset(2);

				} // ipta
			} // iptt 
		} // ic
	} // pp or AA

	int iPTT=2;
	int iPTA=3;

	for(int ic=0;ic<NC;ic++) {
		Filipad *fpad = new Filipad(ic+1, 1.1, 0.4, 100, 100, 0.7, 5);
		fpad->Draw();
		//==== Upper pad
		TPad *p = fpad->GetPad(1); //upper pad
		p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
		TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
		hset( *hfr, "#Delta#eta", "",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
		hfr->Draw();
		//Legend definition
		TLegend *leg = new TLegend(0.25,0.50,0.85,0.82,"","brNDC");
		leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

		cout <<"OK1"<<endl;

		leg->AddEntry((TObject*)NULL,"xx","");
		hDeltaEta[AA][ic][iPTT][iPTA]->SetMarkerStyle(20);
		hDeltaEta[AA][ic][iPTT][iPTA]->Draw("p,same");
		hDeltaEta[pp][0][iPTT][iPTA]->SetMarkerStyle(24);
		hDeltaEta[pp][0][iPTT][iPTA]->Draw("p,same");

		leg->AddEntry(hDeltaEta[AA][ic][iPTT][iPTA],hDeltaEta[AA][ic][iPTT][iPTA]->GetTitle(),"p");
        leg->AddEntry(hDeltaEta[pp][0][iPTT][iPTA],hDeltaEta[pp][0][iPTT][iPTA]->GetTitle(),"p");
		
		leg->Draw();

		cout <<"OK2"<<endl;
		//==== Lower pad
		p = fpad->GetPad(2);
		p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
		TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, -0.35, 0.35);
		hset( *hfr1, "#Delta#eta", "AA/pp",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
		hfr1->Draw();
		//gPad->GetCanvas()->SaveAs(Form("figs_svn/FigA4_v%d_modelcomparisonBest.eps",i+2));
	}
}


TGraphErrors* get_ratio( TGraphErrors * l, TGraphErrors *r ){
	TGraphErrors * gr_ratio = new TGraphErrors( l->GetN() );
	TGraph ger( r->GetN(), r->GetX(), r->GetEY() );
	for( int i=0; i< l->GetN(); i++ ){
		double x = l->GetX()[i];
		double y1 = l->GetY()[i];
		double ey1 = l->GetEY()[i];
		double y2 = r->Eval(x);
		double ey2 = ger.Eval(x);

		double ratio = y1 / y2; 
		gr_ratio->SetPoint( i,  x, ratio);
		gr_ratio->SetPointError( i,  0, ratio*TMath::Sqrt( ey1*ey1/y1/y1+ey2*ey2/y2/y2));
	}   
	return gr_ratio;
}

void SaveCanvas(TString name, TDirectory * dir=0 ){
	//pPrint("figs/"+name,gPad->GetCanvas()->GetName());
	//ppdf("figs.zMixCut3/"+name,gPad->GetCanvas());
	TDirectory *odir = gDirectory;
	name.ReplaceAll(".", "o");
	if( dir ) dir->cd();
	gPad->GetCanvas()->Write("figs/"+name );
	odir->cd();
}


void RemovePoints(TGraphErrors *ge , double xhigh)
{
	// Remove zero points from TGraphErrors.

	if(!ge){return;}

	Int_t nPoints = ge->GetN();
	Double_t x = 0.; 
	Double_t y = 0.; 
	int p =0; 
	while(p<nPoints) {
		ge->GetPoint(p,x,y);
		if( x < 0.21 || x>xhigh  )
			//if( x < 0.21  )
		{   
			ge->RemovePoint(p);
			//cout<<Form(" WARNING (%s): point %d is < 1.e-10 and it was removed from the plot !!!!",ge->GetName(),p+1)<<endl;
			nPoints = ge->GetN();
		} else {
			p++;
		}   
	} // end of for(Int_t p=0;p<nPoints;p++)

	//cout<<endl;
	return;

} // end of void RemoveZeroPoints(TGraphErrors *ge)
