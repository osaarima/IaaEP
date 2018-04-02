void DoAnalysis(double sgnEta=0.2, double bgRbegin=1.0, double bgRend=1.6, double bckScale=1.00, double Side = 0., TString inFile="", TString ppInFile="", TString oname="");
Double_t Gaus2D(Double_t *x, Double_t *par);
void ApplyWingCorrection(TH2D* H, TH1D *hcorr);
void SubtracFlowBackground(TH2D* H, TH1D *hflow);
TGraphErrors* get_ratio( TGraphErrors * l, TGraphErrors *r );
double IntegralOfSmallR( TH2D *hist, double rs, double rb, double &val, double &err, double &area );
double IntegralOfCone( TH2D *hist, double rs, double r0, double r1, double x0, double y0, double &val, double &err, double &area );
double IntegralOfRec( TH2D *hist, double rs, double r0, double r1, double x0, double y0, double &val, double &err, double &area );
double IntegralOfMoon( TH2D *hist, double rs, double rb, double &val, double &err, double &area );
double IntegralForpp( TH2D *hist, double rs, double rb, double &val, double &err, double &area );
void NormalizeToBinWidth2D ( TH2D *hist );
void SaveCanvas(TString name, TDirectory * dir=0 );
double GetGeoAccCorrFlat(double deltaEta);

const int kMAXD       = 20; //maximal number of pT trigger bins
const int kCENT       = 10; //maximal number of pT trigger bins
const int kZvtx       = 15; //maximal number of pT trigger bins

TH1D *hTriggPtBin[2][kCENT][kMAXD]; 
TH1D *hAssocPtBin[2][kCENT][kMAXD][kMAXD]; 
TH2D *hDphiAssoc2DIAA[2][3][kCENT][kMAXD][kMAXD]; // 0:kReal 1:kMixed 2:signal 
TH1D *hWingCorrection[2][kCENT][kMAXD][kMAXD];
TH1D *hFlowBackground[kCENT][kMAXD][kMAXD];
TH2D *h2DIAA[kCENT][kMAXD][kMAXD]; // in eta,phi 
TH1D *hIAAEta[kCENT][kMAXD][kMAXD]; // in eta,phi 
TH1D *hIAAPhi[kCENT][kMAXD][kMAXD]; // in eta,phi 

// before the z average
TH2D *hDphiAssoc2DIAAVtxAA[3][kZvtx][kCENT][kMAXD][kMAXD]; // 0:kReal 1:kMixed 2:signal 
TH2D *hDphiAssoc2DIAAVtxPP[3][kZvtx][kCENT][kMAXD][kMAXD]; // 0:kReal 1:kMixed 2:signal 
TH1D *hTriggPtBinVtx[2][kCENT][kZvtx][kMAXD]; 
TFile *fin[2];
TFile *fmix;
TF2 *g2D;
// adding DeltaEta histograms after mixed event corrections
// with $\Delta\phi$ < 0.2 ??? check
TH1D *hDeltaEta[2][kCENT][kMAXD][kMAXD]; // filipp summed DeltaEta AA-1 pp-0
// save this into an additional root file, fit and IAA will be calculated in z02.CalIAADeta.C
Bool_t saveDeta = kTRUE;


void run(){

	const int NAA = 1;
	TString fileAA[NAA] = {
		"legotrain_JCIaa/data/JCIaa_legotrain_PbPb_CF-5059_20180329-1139_runlist_3-LHC10h_AOD86_MgFpMgFm.root"
	};
	TString commentAA[NAA] = {
		"AA_moon1"
	};

	const int NPP = 1;
	TString filePP[NPP] = {
		"legotrain_JCIaa/data/JCIaa_legotrain_CF_pp-1677_20180326-1905-2760GeV_LHC11a_p4_AOD113_noSDD.root"
	};
	TString commentPP[NPP] = {
		"pp_moon1"
	};

	// Moon
	const int NR = 1;
	double dR[NR] = {0.2};
	double BgRbegin[1] = {1.0};
	//double BgRbegin[6] = {1.0,1.1,1.2,1.3,1.4};
	int NBG=1;
	for(int iA=0;iA<NAA;iA++) {
		for(int iP=0;iP<NPP;iP++) {
			for(int iR=0;iR<NR;iR++){
				for( int iB=0;iB<NBG;iB++){
					DoAnalysis( dR[iR], BgRbegin[iB], 1.6, 1, 0, fileAA[iA],filePP[iP],commentAA[iA]+"_"+commentPP[iP] );
				}
			}
		}
	}
}


void DoAnalysis(double sgnEta=0.2, double bgRbegin=1.0, double bgRend=1.6, double bckScale=1.00, double Side = 0., TString inFile="", TString ppInFile="", TString oname=""){

	int applyWingCorrection = 0;
	int doMixMerge = 1;
	int takeMixExt = 0;
	int correctMix = 1;
	double mixptt = 4.0, mixpta = 3.0;
	// this is just a cross check

	enum fillType { kReal, kMixed, kSignal };
	enum dataType { AA, pp };
	fillType fTyp;
	cout <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	cout <<"PbPb : "<< inFile << endl;
	cout <<"pp :"<< ppInFile << endl;
	cout <<"Discription :"<< oname << endl;

	fin[AA] = TFile::Open(inFile);
	fin[pp] = TFile::Open(ppInFile);
	//fmix = new TFile("mergedmixDJ_Cut1.root","read");

	TString TopDir[2] = {"JCIAA_TPCOnly_H0_T0","JCIAA_GlobalSDD_H0_T0"};

	double bgnEta[2] = {bgRbegin,bgRend}; //  ietaBckThr=5 R=1.0

	TVector *TriggPtBorders[2];
	TVector *AssocPtBorders[2];
	TVector *CentBinBorders[2];
	TVector *EtaGapThresholds[2];
	TVector *zVertBins[2];

	for(int idtyp=0; idtyp<2; idtyp++){ 
		TriggPtBorders[idtyp]              = (TVector*) fin[idtyp]->Get(TopDir[idtyp]+"/JCard/TriggPtBorders");
		AssocPtBorders[idtyp]              = (TVector*) fin[idtyp]->Get(TopDir[idtyp]+"/JCard/AssocPtBorders");
		CentBinBorders[idtyp]              = (TVector*) fin[idtyp]->Get(TopDir[idtyp]+"/JCard/CentBinBorders");
		EtaGapThresholds[idtyp]			   = (TVector*) fin[idtyp]->Get(TopDir[idtyp]+"/JCard/EtaGapThresholds");
		zVertBins[idtyp]				   = (TVector*) fin[idtyp]->Get(TopDir[idtyp]+"/JCard/zVertBins");
	}
	// Temporary
	int nzvtx[2] = { 1, zVertBins[pp]->GetNoElements()-1}; 


	cout <<"++++++++++++++++  Settings +++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	cout <<"R_signal <  "<< sgnEta <<", "<< bgnEta[0] <<"<R_bck<"<< bgnEta[1] << endl;
	cout <<"Wing Correction = "<< applyWingCorrection << endl;
	cout <<"Merging pt bins of mixed event (ptt>"<<mixptt<<",pta>"<<mixpta<<") = "<< doMixMerge << endl;
	cout <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;


	enum subType { kNormal };
	subType method = kNormal;

	const char *sidelabel[2] = {"Near","Far"};
	//=========================================
	int icbin[2]    = {0, 4};
	int reference = 0; //0=pp, 1=60-90%
	int subtract = 1;
	int logy = 1;
	int plotPythia = 1;
	int plotIAA = 1;
	int saveRoot = 1;  
	TString prefix = "sysErrors/";
	double ymin1 = 1e-3; //1e-3; 
	double ymax1 = subtract==0 ? 1e4 : 20; //40;
	double ymax2 = subtract==0 ? 35 : 6;  //35 150
	int sdx = 250;
	double lpta = 0.;
	double hpta = 10.;
	double hIAA = 5.;
	// for merging mixed event

	//double v3PhaseCorr = 0.948; // Needs to be finished
	double v3PhaseCorr = 1.0; //0.95; // Needs to be finished
	int collors[5] ={920, 622, 406, 590, 606};
	//=========================================

	g2D = new TF2("g2D",Gaus2D,-1.6,1.6,-0.4,0.4,5);
	g2D->SetParNames("Const","#Delta#eta_{0}","#sigma_{#Delta#eta}","#Delta#phi_{0}","#sigma_{#Delta#phi}");
	g2D->SetParameters(5,0,0.3,0,0.3);

	double ntrigg = 0;

	int NumEtaGaps = EtaGapThresholds[pp]->GetNoElements()-1; 
	cout <<"pp"<<endl;
	cout <<"bins:  "<<" eta="<< NumEtaGaps <<" zvtx="<<nzvtx[pp]<< endl; 

	int NumCent[2]    = { CentBinBorders[AA]->GetNoElements()-1, 1}; 
	NumEtaGaps = EtaGapThresholds[AA]->GetNoElements()-1; 
	int NumPtt     = TriggPtBorders[AA]->GetNoElements()-1;
	int NumPta     = AssocPtBorders[AA]->GetNoElements()-1;
	cout <<"PbPb"<<endl;
	cout <<"bins:  c="<<  NumCent[0] <<" eta="<< NumEtaGaps <<" ptt="<< NumPtt <<" pta="<< NumPta  <<" zvtx="<<nzvtx[AA]<< endl; 
	cout <<"+++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

	const int NC = NumCent[0];
	const int NPTT = NumPtt;
	const int NPTA = NumPta;
	double MeanPta[2][NC][NPTT][NPTA];

	double InclYieldAA[NC][NPTT][NPTA];
	double InclYieldpp[NPTT][NPTA];
	double eInclYieldAA[NC][NPTT][NPTA];
	double eInclYieldpp[NPTT][NPTA];

	double BackgroundYieldAA[NC][NPTT][NPTA];
	double BackgroundYieldpp[NPTT][NPTA];
	double eBackgroundYieldAA[NC][NPTT][NPTA];
	double eBackgroundYieldpp[NPTT][NPTA];

	double SignalYieldAA[NC][NPTT][NPTA];
	double SignalYieldpp[NPTT][NPTA];
	double eSignalYieldAA[NC][NPTT][NPTA];
	double eSignalYieldpp[NPTT][NPTA];

	double InclYieldAAFar[NC][NPTT][NPTA];
	double InclYieldppFar[NPTT][NPTA];
	double eInclYieldAAFar[NC][NPTT][NPTA];
	double eInclYieldppFar[NPTT][NPTA];

	double BackgroundYieldAAFar[NC][NPTT][NPTA];
	double BackgroundYieldppFar[NPTT][NPTA];
	double eBackgroundYieldAAFar[NC][NPTT][NPTA];
	double eBackgroundYieldppFar[NPTT][NPTA];

	double SignalYieldAAFar[NC][NPTT][NPTA];
	double SignalYieldppFar[NPTT][NPTA];
	double eSignalYieldAAFar[NC][NPTT][NPTA];
	double eSignalYieldppFar[NPTT][NPTA];
	// === normalization =================
	// Near: signal in first Rgap bin
	// Far:  Signal in second phi-pi bin
	// ===================================
	double MixedEventStatAA[NC][NPTT][NPTA];
	double MixedEventStatPP[NPTT][NPTA];

	//------------ R e a d    D a t a ------------    
	cout <<"Reading data...."<<endl;
	for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp  //(*fCentralityBinBorders)[i+1]
		cout << idtyp <<"\t"<< nzvtx[idtyp] << endl;
		for(int ic=0; ic<NumCent[idtyp]; ic++){
			for(int iptt=0; iptt<NumPtt; iptt++){
				if(idtyp==AA) for(int iz=0; iz<nzvtx[idtyp]; iz++) hTriggPtBinVtx[idtyp][ic][iz][iptt] = (TH1D*) fin[idtyp]->Get(Form("%s/AliJHistos/hTriggPtBin/hTriggPtBinC%02dV%02dT%02d",TopDir[idtyp].Data(),ic, 05, iptt));
				if(idtyp==pp) for(int iz=0; iz<nzvtx[idtyp]; iz++) hTriggPtBinVtx[idtyp][ic][iz][iptt] = (TH1D*) fin[idtyp]->Get(Form("%s/AliJHistos/hTriggPtBin/hTriggPtBinC%02dV%02dT%02d",TopDir[idtyp].Data(),ic, iz, iptt));
				for(int ipta=0;ipta<NumPta;ipta++) {
					//hAssocPtBin[idtyp][ic][iptt][ipta] = (TH1D*) fin[idtyp]->Get(Form("%s/AliJHistos/hAssocPtBin/hAssocPtBin%02d%02d%02d", TopDir[idtyp].Data(),ic, iptt,ipta));//distribution of trigger pT
					//MeanPta[idtyp][ic][iptt][ipta] = hAssocPtBin[idtyp][ic][iptt][ipta]->GetMean();
					MeanPta[idtyp][ic][iptt][ipta] = (*AssocPtBorders[AA])[ipta+1] + (*AssocPtBorders[AA])[ipta+2]/2.;
					for(int ityp=0; ityp<2; ityp++){
						if(idtyp==AA) {
							for(int iz=0; iz<nzvtx[idtyp]; iz++){
								//cout << Form("%s/hDphiDetaPta/hDphiDetaPtaD%02dC%02dV%02dT%02dA%02d",TopDir[idtyp].Data(), ityp, ic, 5, iptt, ipta) << endl;
								hDphiAssoc2DIAAVtxAA[ityp][iz][ic][iptt][ipta] = (TH2D *)fin[idtyp]->Get(Form("%s/AliJHistos/hDphiDetaPta/hDphiDetaPtaD%02dC%02dV%02dT%02dA%02d",TopDir[idtyp].Data(), ityp, ic, 5, iptt, ipta));
								hDphiAssoc2DIAAVtxAA[ityp][iz][ic][iptt][ipta]->SetXTitle("#Delta#eta");
								hDphiAssoc2DIAAVtxAA[ityp][iz][ic][iptt][ipta]->GetXaxis()->CenterTitle(kTRUE);
								hDphiAssoc2DIAAVtxAA[ityp][iz][ic][iptt][ipta]->GetXaxis()->SetTitleOffset(2);
								hDphiAssoc2DIAAVtxAA[ityp][iz][ic][iptt][ipta]->SetYTitle("#Delta#phi/#pi");
								hDphiAssoc2DIAAVtxAA[ityp][iz][ic][iptt][ipta]->GetYaxis()->CenterTitle(kTRUE);
								hDphiAssoc2DIAAVtxAA[ityp][iz][ic][iptt][ipta]->GetYaxis()->SetTitleOffset(2);
								//hDphiAssoc2DIAAVtxAA[ityp][iz][ic][iptt][ipta]->SetZTitle("#frac{1}{N_{trigg}} #frac{1}{d#Delta#eta#Delta#phi}");
								hDphiAssoc2DIAAVtxAA[ityp][iz][ic][iptt][ipta]->GetZaxis()->CenterTitle(kTRUE);
								hDphiAssoc2DIAAVtxAA[ityp][iz][ic][iptt][ipta]->GetZaxis()->SetTitleOffset(2);
								hDphiAssoc2DIAAVtxAA[ityp][iz][ic][iptt][ipta]->Rebin2D(3,3);
							}
						}
						if(idtyp==pp) {
							for(int iz=0; iz<nzvtx[idtyp]; iz++){
								hDphiAssoc2DIAAVtxPP[ityp][iz][ic][iptt][ipta] = (TH2D *)fin[idtyp]->Get(Form("%s/AliJHistos/hDphiDetaPta/hDphiDetaPtaD%02dC%02dV%02dT%02dA%02d", TopDir[idtyp].Data(), ityp, ic, iz, iptt, ipta));
								hDphiAssoc2DIAAVtxPP[ityp][iz][ic][iptt][ipta]->SetXTitle("#Delta#eta");
								hDphiAssoc2DIAAVtxPP[ityp][iz][ic][iptt][ipta]->GetXaxis()->CenterTitle(kTRUE);
								hDphiAssoc2DIAAVtxPP[ityp][iz][ic][iptt][ipta]->GetXaxis()->SetTitleOffset(2);
								hDphiAssoc2DIAAVtxPP[ityp][iz][ic][iptt][ipta]->SetYTitle("#Delta#phi/#pi");
								hDphiAssoc2DIAAVtxPP[ityp][iz][ic][iptt][ipta]->GetYaxis()->CenterTitle(kTRUE);
								hDphiAssoc2DIAAVtxPP[ityp][iz][ic][iptt][ipta]->GetYaxis()->SetTitleOffset(2);
								//hDphiAssoc2DIAAVtxPP[ityp][iz][ic][iptt][ipta]->SetZTitle("#frac{1}{N_{trigg}} #frac{1}{d#Delta#eta#Delta#phi}");
								hDphiAssoc2DIAAVtxPP[ityp][iz][ic][iptt][ipta]->GetZaxis()->CenterTitle(kTRUE);
								hDphiAssoc2DIAAVtxPP[ityp][iz][ic][iptt][ipta]->GetZaxis()->SetTitleOffset(2);
								hDphiAssoc2DIAAVtxPP[ityp][iz][ic][iptt][ipta]->Rebin2D(3,3);
							}
						}
					} // type
				} // ipta
			} // iptt 
		} // ic
	} // pp or AA

	int imixpta = 0;
	int imixptt = 0;
	if(doMixMerge) {
		cout <<"Merging event mixing bin...."<<endl;
		int imixpta = 0;
		int imixptt = 0;
		for(int iptt=0; iptt<NumPtt; iptt++){
			if((*TriggPtBorders[AA])[iptt+1]>=mixptt ) { imixptt=iptt;break; }
		}
		for(int ipta=0;ipta<NumPta;ipta++) {
			if((*AssocPtBorders[AA])[ipta+1]>=mixpta) { imixpta=ipta;break; }
		}

		cout << " merge after iptt="<< imixptt<<"\t ipta="<< imixpta << endl; 
		for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp
			for(int ic=0; ic<NumCent[idtyp]; ic++){
				for(int iz=0; iz<nzvtx[idtyp]; iz++){
					if(idtyp==AA) {
						for(int iiptt=imixptt+1; iiptt<NumPtt; iiptt++) {
							for(int iipta=imixpta+1;iipta<NumPta;iipta++) {	
								hDphiAssoc2DIAAVtxAA[kMixed][iz][ic][imixptt][imixpta]->Add(hDphiAssoc2DIAAVtxAA[kMixed][iz][ic][iiptt][iipta]);
							}
						}
					} //AA 
					if(idtyp==pp) {
						for(int iiptt=imixptt+1; iiptt<NumPtt; iiptt++) {
							for(int iipta=imixpta+1;iipta<NumPta;iipta++) {	
								hDphiAssoc2DIAAVtxPP[kMixed][iz][ic][imixptt][imixpta]->Add(hDphiAssoc2DIAAVtxPP[kMixed][iz][ic][iiptt][iipta]);
							}
						}
					} // pp 
				} 
			} 
		} 

		// Copying the higher bins from merged one.
		if(takeMixExt) {
			for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp
				for(int ic=0; ic<NumCent[idtyp]; ic++){
					for(int iptt=imixptt; iptt<NumPtt; iptt++){
						for(int ipta=imixpta;ipta<NumPta;ipta++) {
							for(int iz=0; iz<nzvtx[idtyp]; iz++){
								if(idtyp==AA)	hDphiAssoc2DIAAVtxAA[kMixed][iz][ic][iptt][ipta] = (TH2D*)fmix->Get(Form("hDphiAssoc2DIAAVtxAAMixZ%02dC%02d_mergedmixAA",iz,ic)); 
								if(idtyp==pp)	hDphiAssoc2DIAAVtxPP[kMixed][iz][ic][iptt][ipta] = (TH2D*)fmix->Get(Form("hDphiAssoc2DIAAVtxAAMixZ%02dC%02d_mergedmixPP",iz,ic)); 
							}
						}
					}
				}
			} 
		} else {	
			for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp
				for(int ic=0; ic<NumCent[idtyp]; ic++){
					for(int iz=0; iz<nzvtx[idtyp]; iz++){
						for(int iptt=imixptt; iptt<NumPtt; iptt++){
							for(int ipta=imixpta;ipta<NumPta;ipta++) {
								if(idtyp==AA)	hDphiAssoc2DIAAVtxAA[kMixed][iz][ic][iptt][ipta] = (TH2D*)hDphiAssoc2DIAAVtxAA[kMixed][iz][ic][imixptt][imixpta]->Clone();
								if(idtyp==pp)	hDphiAssoc2DIAAVtxPP[kMixed][iz][ic][iptt][ipta] = (TH2D*)hDphiAssoc2DIAAVtxPP[kMixed][iz][ic][imixptt][imixpta]->Clone();
							} 
						} 
					}
				}
			} 
		}
		if(0) {
			// Wrrite down the merged bin
			TFile *fmixout = new TFile("mergedmixDJ_Cut1.root","recreate");
			for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp
				for(int ic=0; ic<NumCent[idtyp]; ic++){
					for(int iz=0; iz<nzvtx[idtyp]; iz++){
						fmixout->cd();
						if(idtyp==AA) hDphiAssoc2DIAAVtxAA[kMixed][iz][ic][imixptt][imixpta]->Write(Form("hDphiAssoc2DIAAVtxAAMixZ%02dC%02d_mergedmixAA",iz,ic));
						if(idtyp==pp) hDphiAssoc2DIAAVtxPP[kMixed][iz][ic][imixptt][imixpta]->Write(Form("hDphiAssoc2DIAAVtxAAMixZ%02dC%02d_mergedmixPP",iz,ic));
					} 
				}
			} 
		}

		// for drawing of merged mixed event
	} // doMixMerge

	// make integrated one here
	cout <<"Merging mixed event z bins into one bin for drawing "<<endl;
	for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp
		for(int ic=0; ic<NumCent[idtyp]; ic++){
			for(int iptt=0; iptt<NumPtt; iptt++){
				for(int ipta=0;ipta<NumPta;ipta++) {
					if(idtyp==AA)	hDphiAssoc2DIAA[idtyp][kMixed][ic][iptt][ipta] = (TH2D*)hDphiAssoc2DIAAVtxAA[kMixed][0][ic][iptt][ipta]->Clone();
					if(idtyp==pp)	hDphiAssoc2DIAA[idtyp][kMixed][ic][iptt][ipta] = (TH2D*)hDphiAssoc2DIAAVtxPP[kMixed][0][ic][iptt][ipta]->Clone();
					for(int iz=1; iz<nzvtx[idtyp]; iz++){
						if(idtyp==AA) hDphiAssoc2DIAA[idtyp][kMixed][ic][iptt][ipta]->Add(hDphiAssoc2DIAAVtxAA[kMixed][iz][ic][iptt][ipta]);
						if(idtyp==pp) hDphiAssoc2DIAA[idtyp][kMixed][ic][iptt][ipta]->Add(hDphiAssoc2DIAAVtxPP[kMixed][iz][ic][iptt][ipta]);
					} 
				} 
			}
		}
	} 

	cout <<"Mixed event correction"<<endl;
	//------------ Mixed event correction ------------    
	for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp
		for(int ic=0; ic<NumCent[idtyp]; ic++){
			for(int iptt=0; iptt<NumPtt; iptt++){
				for(int ipta=0;ipta<NumPta;ipta++) {
					if(idtyp==AA) {
						double nmixed = 0.;
						for(int iz=0; iz<nzvtx[idtyp]; iz++){
							hDphiAssoc2DIAAVtxAA[kSignal][iz][ic][iptt][ipta] = (TH2D*)hDphiAssoc2DIAAVtxAA[kReal][iz][ic][iptt][ipta]->Clone();
							double norm  = hDphiAssoc2DIAAVtxAA[kMixed][iz][ic][iptt][ipta]->Integral(); // should be before binwidth co
							nmixed+= norm;
							NormalizeToBinWidth2D ( hDphiAssoc2DIAAVtxAA[kMixed][iz][ic][iptt][ipta] );
							hDphiAssoc2DIAAVtxAA[kMixed][iz][ic][iptt][ipta]->Scale(2*1.6/norm);
							if(correctMix) hDphiAssoc2DIAAVtxAA[kSignal][iz][ic][iptt][ipta]->Divide(hDphiAssoc2DIAAVtxAA[kMixed][iz][ic][iptt][ipta]);
						} // z bin
						MixedEventStatAA[ic][iptt][ipta] = nmixed;
					} // AA
					if(idtyp==pp) {
						double nmixed = 0.;
						for(int iz=0; iz<nzvtx[idtyp]; iz++){
							hDphiAssoc2DIAAVtxPP[kSignal][iz][ic][iptt][ipta] = (TH2D*)hDphiAssoc2DIAAVtxPP[kReal][iz][ic][iptt][ipta]->Clone();
							double norm  = hDphiAssoc2DIAAVtxPP[kMixed][iz][ic][iptt][ipta]->Integral(); // should be before binwidth co
							nmixed+= norm;
							NormalizeToBinWidth2D ( hDphiAssoc2DIAAVtxPP[kMixed][iz][ic][iptt][ipta] );
							hDphiAssoc2DIAAVtxPP[kMixed][iz][ic][iptt][ipta]->Scale(2*1.6/norm);
							if(correctMix) hDphiAssoc2DIAAVtxPP[kSignal][iz][ic][iptt][ipta]->Divide(hDphiAssoc2DIAAVtxPP[kMixed][iz][ic][iptt][ipta]);
						} // z bin
						MixedEventStatPP[iptt][ipta] = nmixed;
					} // PP
				} // ipta
			} // iptt 
		} // ic
	} // pp or AA

	cout <<"++++++++++++++++++++ Mixed event statistics +++++++++++++++++++++++"<< endl;
	cout << " collision cent ptt  pta   entries "<<endl;
	for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp
		for(int ic=0; ic<NumCent[idtyp]; ic++){
			for(int iptt=0; iptt<NumPtt; iptt++){
				for(int ipta=0;ipta<NumPta;ipta++) {
					if((*AssocPtBorders[AA])[ipta+2]<=(*TriggPtBorders[AA])[iptt+1]) {
						if(doMixMerge) { if(iptt==imixptt&&ipta==imixpta) cout <<"-----------------------------------------------------"<<endl;}
						if(idtyp==AA) printf("%5s  %5.0f  %5.1f  %5.1f  %e\n","AA",(*CentBinBorders[AA])[ic+1],(*TriggPtBorders[AA])[iptt+1],(*AssocPtBorders[AA])[ipta+1],MixedEventStatAA[ic][iptt][ipta]);
						if(idtyp==pp) printf("%5s  %5.0f  %5.1f  %5.1f  %e\n","PP",(*CentBinBorders[AA])[ic+1],(*TriggPtBorders[AA])[iptt+1],(*AssocPtBorders[AA])[ipta+1],MixedEventStatPP[iptt][ipta]);
					}
				}
			}
		}
	}
	//--------------------------------------------------------------
	// calculating a weighted average
	//--------------------------------------------------------------
	cout <<"Calculating the weighted average.. over zvertex bins"<<endl;
	// Normalized per trigger yields
	// AA
	for(int idtyp=0; idtyp<1; idtyp++){ // 0 = AA, 1 = pp
		for(int ic=0; ic<NumCent[idtyp]; ic++){
			for(int iptt=0; iptt<NumPtt; iptt++){
				hTriggPtBin[idtyp][ic][iptt] = (TH1D*)hTriggPtBinVtx[idtyp][ic][0][iptt]->Clone();
				for(int iz=1; iz<nzvtx[idtyp]; iz++) hTriggPtBin[idtyp][ic][iptt]->Add(hTriggPtBinVtx[idtyp][ic][iz][iptt]);
				for(int ipta=0;ipta<NumPta;ipta++) {
					if(idtyp==AA) {
						hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta] = (TH2D*)hDphiAssoc2DIAAVtxAA[kSignal][0][ic][iptt][ipta]->Clone();
						for(int iz=1; iz<nzvtx[idtyp]; iz++){
							hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta]->Add(hDphiAssoc2DIAAVtxAA[kSignal][iz][ic][iptt][ipta]);
						}
						double ntriggall = hTriggPtBin[idtyp][ic][iptt]->Integral();
						cout <<"Number of Trigger particles "<< ic <<"\t"<< (*TriggPtBorders[AA])[iptt+1] <<"<ptt<"<< (*TriggPtBorders[AA])[iptt+2] <<"\t"<< ntriggall << endl;
						hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta]->Scale(1./ntriggall);
						// Check the wing correction
						int phil = hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta]->GetYaxis()->FindBin(1-0.5);
						int phih = hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta]->GetYaxis()->FindBin(1+0.5);
						TString hname = Form("hWingCorrection%02d%02d%02d%02d",idtyp,ic,iptt,ipta);
						hWingCorrection[idtyp][ic][iptt][ipta] = (TH1D*) hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta]->ProjectionX(hname.Data(),phil,phih);
						double WingNorm = hWingCorrection[idtyp][ic][iptt][ipta]->Integral();
						double bw = hWingCorrection[idtyp][ic][iptt][ipta]->GetBinWidth(1);
						hWingCorrection[idtyp][ic][iptt][ipta]->Scale(2*1.6/bw/WingNorm);
						if(applyWingCorrection && idtyp==AA) ApplyWingCorrection(hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta], hWingCorrection[idtyp][ic][iptt][ipta]);
					}	
				}
			}
		}
	}
	// PP
	for(int idtyp=1; idtyp<2; idtyp++){ // 0 = AA, 1 = pp
		for(int ic=0; ic<NumCent[idtyp]; ic++){
			for(int iptt=0; iptt<NumPtt; iptt++){
				hTriggPtBin[idtyp][ic][iptt] = (TH1D*)hTriggPtBinVtx[idtyp][ic][0][iptt]->Clone();
				for(int iz=1; iz<nzvtx[idtyp]; iz++) hTriggPtBin[idtyp][ic][iptt]->Add(hTriggPtBinVtx[idtyp][ic][iz][iptt]);
				for(int ipta=0;ipta<NumPta;ipta++) {
					if(idtyp==pp) {
						hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta] = (TH2D*)hDphiAssoc2DIAAVtxPP[kSignal][0][ic][iptt][ipta]->Clone();
						for(int iz=1; iz<nzvtx[idtyp]; iz++){
							hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta]->Add(hDphiAssoc2DIAAVtxPP[kSignal][iz][ic][iptt][ipta]);
						}
						double ntriggall = hTriggPtBin[idtyp][ic][iptt]->Integral();
						hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta]->Scale(1./ntriggall);
						// Check the wing correction
						int phil = hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta]->GetYaxis()->FindBin(1-0.5);
						int phih = hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta]->GetYaxis()->FindBin(1+0.5);
						TString hname = Form("hWingCorrection%02d%02d%02d%02d",idtyp,ic,iptt,ipta);
						hWingCorrection[idtyp][ic][iptt][ipta] = (TH1D*) hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta]->ProjectionX(hname.Data(),phil,phih);
						double WingNorm = hWingCorrection[idtyp][ic][iptt][ipta]->Integral();
						double bw = hWingCorrection[idtyp][ic][iptt][ipta]->GetBinWidth(1);
						hWingCorrection[idtyp][ic][iptt][ipta]->Scale(2*1.6/bw/WingNorm);
					}
				}	
			}
		}
	}

	cout <<"Calculating Signal yields.. for integrated IAA.."<<endl;
	// Now calcuate yields
	for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp
		for(int ic=0; ic<NumCent[idtyp]; ic++){
			for(int iptt=0; iptt<NumPtt; iptt++){
				ntrigg = hTriggPtBin[idtyp][ic][iptt]->Integral();
				for(int ipta=0;ipta<NumPta;ipta++) {
					int phil = hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta]->GetYaxis()->FindBin(1-0.2);
					int phih = hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta]->GetYaxis()->FindBin(1+0.2);
					// Signal
					double  IncYield, eIncYield, IncYield_Area;
					IntegralOfCone( hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta], sgnEta,0., sgnEta, 0, Side , IncYield, eIncYield, IncYield_Area);
					// Background
					double  BckYield, eBckYield, BckYield_Area  ;
					//if(idtyp==pp) IntegralForpp( hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta], sgnEta, 1.0, BckYield, eBckYield, BckYield_Area );
					//hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta]->Fit("g2D","V");
					if (sgnEta<=0.3) {
						IntegralOfSmallR( hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta], sgnEta, bgnEta[0], BckYield, eBckYield, BckYield_Area ); // if(R<0.3     ) circle not moon
					} else {
						IntegralOfMoon( hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta], sgnEta, bgnEta[0], BckYield, eBckYield, BckYield_Area );
					}
					//if(ic==4) IntegralForpp( hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta], sgnEta, bgnEta[0], BckYield, eBckYield, BckYield_Area );

					double dpt = ((*AssocPtBorders[AA])[ipta+2]-(*AssocPtBorders[AA])[ipta+1]);
					double nearMixedEventRatio = 1.;

					// NORM **WARNING** => absolute norm can give 4% difference because of binning in case R=0.5
					double NearEtaBinNorm = IncYield_Area/BckYield_Area;
					double NormBckYield = BckYield * NearEtaBinNorm * nearMixedEventRatio * bckScale;
					double SignalYield = (IncYield - NormBckYield)/dpt;
					double eSignalYield = TMath::Sqrt(eIncYield*eIncYield + eBckYield*eBckYield)/dpt;
					if(idtyp==AA) { 
						InclYieldAA[ic][iptt][ipta] = IncYield; eInclYieldAA[ic][iptt][ipta] = eIncYield ; 
						BackgroundYieldAA[ic][iptt][ipta] = NormBckYield; eBackgroundYieldAA[ic][iptt][ipta] = eBckYield * NearEtaBinNorm * bckScale ; 
						SignalYieldAA[ic][iptt][ipta] = SignalYield; eSignalYieldAA[ic][iptt][ipta] = eSignalYield ; 
					}
					if(idtyp==pp) { 
						InclYieldpp[iptt][ipta] = IncYield; eInclYieldpp[iptt][ipta] = eIncYield;
						BackgroundYieldpp[iptt][ipta] = NormBckYield; eBackgroundYieldpp[iptt][ipta] = eBckYield * NearEtaBinNorm * bckScale ; 
						SignalYieldpp[iptt][ipta] = SignalYield; eSignalYieldpp[iptt][ipta] = eSignalYield;
					}
					if(idtyp==AA){
						h2DIAA[ic][iptt][ipta] = (TH2D*)hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta]->Clone();
						h2DIAA[ic][iptt][ipta]->Divide(hDphiAssoc2DIAA[pp][kSignal][0][iptt][ipta]);
						int etal = h2DIAA[ic][iptt][ipta]->GetXaxis()->FindBin(-0.4);
						int etah = h2DIAA[ic][iptt][ipta]->GetXaxis()->FindBin(0.4);
						int phil = h2DIAA[ic][iptt][ipta]->GetYaxis()->FindBin(-0.4);
						int phih = h2DIAA[ic][iptt][ipta]->GetYaxis()->FindBin(0.4);
						hIAAEta[ic][iptt][ipta] = (TH1D*)h2DIAA[ic][iptt][ipta]->ProjectionX(Form("hIAAEtaC%02dT%02dA%02d",ic,iptt,ipta),phil,phih); // phi near
						hIAAPhi[ic][iptt][ipta] = (TH1D*)h2DIAA[ic][iptt][ipta]->ProjectionY(Form("hIAAPhiC%02dT%02dA%02d",ic,iptt,ipta),etal,etah); // eta nea,
					}
				} // pta
			} // ptt 
		} // cent
	} // type 
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	cout <<"Calculating DeltaEta, just projection "<<endl;
	for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp
		for(int ic=0; ic<NumCent[idtyp]; ic++){
			for(int iptt=0; iptt<NumPtt; iptt++){
				ntrigg = hTriggPtBin[idtyp][ic][iptt]->Integral();
				for(int ipta=0;ipta<NumPta;ipta++) {
					int phil = hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta]->GetYaxis()->FindBin(1-0.2);
					int phih = hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta]->GetYaxis()->FindBin(1+0.2);
					// Signal
					double dpt = ((*AssocPtBorders[AA])[ipta+2]-(*AssocPtBorders[AA])[ipta+1]);
					hDeltaEta[idtyp][ic][iptt][ipta] = (TH1D*)hDphiAssoc2DIAA[idtyp][kSignal][ic][iptt][ipta]->ProjectionX(Form("hDeltaEtaType%02dC%02dT%02dA%02d",idtyp,ic,iptt,ipta),phil,phih); // phi near
				} // pta
			} // ptt 
		} // cent
	} // type 


	// Make TGraphs , yield and IAA
	cout <<"Making graphs of yields and IAA..."<<endl;
	TGraphErrors *grSignalYieldpp[NPTT];
	TGraphErrors *grSignalYieldAA[NC][NPTT];

	TGraphErrors *grInclYieldpp[NPTT];
	TGraphErrors *grInclYieldAA[NC][NPTT];

	TGraphErrors *grBckYieldpp[NPTT];
	TGraphErrors *grBckYieldAA[NC][NPTT];

	TGraphErrors *grIAA[NC][NPTT];

	for(int iptt=0;iptt<NPTT;iptt++) {
		grSignalYieldpp[iptt] = new TGraphErrors(NPTA, MeanPta[0][0][iptt], SignalYieldpp[iptt],0,eSignalYieldpp[iptt]);
		grInclYieldpp[iptt] = new TGraphErrors(NPTA, MeanPta[0][0][iptt], InclYieldpp[iptt],0,eInclYieldpp[iptt]);
		grBckYieldpp[iptt] = new TGraphErrors(NPTA, MeanPta[0][0][iptt], BackgroundYieldpp[iptt],0,eBackgroundYieldpp[iptt]);
		for(int ic=0;ic<NC;ic++) {
			// ?? meanpt is zero for perp
			grInclYieldAA[ic][iptt] = new TGraphErrors(NPTA, MeanPta[1][0][iptt], InclYieldAA[ic][iptt],0,eInclYieldAA[ic][iptt]);
			grBckYieldAA[ic][iptt] = new TGraphErrors(NPTA, MeanPta[1][0][iptt], BackgroundYieldAA[ic][iptt],0,eBackgroundYieldAA[ic][iptt]);
			grSignalYieldAA[ic][iptt] = new TGraphErrors(NPTA, MeanPta[1][0][iptt], SignalYieldAA[ic][iptt],0,eSignalYieldAA[ic][iptt]);
			grIAA[ic][iptt] = get_ratio(grSignalYieldAA[ic][iptt],grSignalYieldpp[iptt]);
		}
	}
	// Write down Iaa and yields into a root file
	if(saveRoot) {
		cout <<"Writing the results into a file..."<< endl;
		TFile *fout = new TFile(Form("%s_%s_Iaa_R%.1f_%.1f_%.2f_%s_Wing%d.root",prefix.Data(),oname.Data(),sgnEta,bgnEta[0],bgnEta[1],sidelabel[int(Side)],applyWingCorrection),"recreate");
		fout->cd();
		for(int iptt=0;iptt<NPTT;iptt++) {
			grInclYieldpp[iptt]->Write(Form("grInclYieldpp_%02d",iptt));
			grBckYieldpp[iptt]->Write(Form("grBckYieldpp_%02d",iptt));
			grSignalYieldpp[iptt]->Write(Form("grYieldpp_%02d",iptt));
			for(int ic=0;ic<NC;ic++) {
				grInclYieldAA[ic][iptt]->Write(Form("grInclYieldAA_%02d%02d",ic,iptt));
				grBckYieldAA[ic][iptt]->Write(Form("grBckYieldAA_%02d%02d",ic,iptt));
				grSignalYieldAA[ic][iptt]->Write(Form("grYieldAA_%02d%02d",ic,iptt));
				grIAA[ic][iptt]->Write(Form("grIAA_%02d%02d",ic,iptt));
			}
		}
		// Deltaeta
		for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp
			for(int ic=0; ic<NumCent[idtyp]; ic++){
				for(int iptt=0; iptt<NumPtt; iptt++){
					for(int ipta=0;ipta<NumPta;ipta++) {
						hDeltaEta[idtyp][ic][iptt][ipta]->Write();
					} // pta
				} // ptt 
			} // cent
		} // type 
		fout->cd();
		//WriteJCard(fin[AA],fout);
		fout->Close();
	}

	cout <<"Closing files.."<<endl;
	// When a file is closed, all histograms in memory associated with this file are automatically deleted.
	fin[pp]->Close();
	fin[AA]->Close();
	//fmix->Close();
}



void ApplyWingCorrection(TH2D* H, TH1D *hcorr){ 
	// normalize each bin with the bin width
	int nbx = H->GetNbinsX();
	int nby = H->GetNbinsY();

	for(int i=1;i<=nbx;i++){ //loop over all bins
		for(int j=1;j<=nby;j++){ //loop over all bins
			double bc = H->GetBinContent(i,j);
			double be = H->GetBinError(i,j);
			double corr = hcorr->GetBinContent(i);
			H->SetBinContent(i,j,bc/corr);
			H->SetBinError(i,j,be/corr);
		}
	}    
}

void NormalizeToBinWidth2D(TH2D* H){ 
	// normalize each bin with the bin width
	int nbx = H->GetNbinsX();
	int nby = H->GetNbinsY();

	for(int i=1;i<=nbx;i++){ //loop over all bins
		for(int j=1;j<=nby;j++){ //loop over all bins
			double bc = H->GetBinContent(i,j);
			double be = H->GetBinError(i,j);
			double bwx = H->GetXaxis()->GetBinWidth(i);
			double bwy = H->GetYaxis()->GetBinWidth(j);

			H->SetBinContent(i,j,bc/bwx/bwy);
			H->SetBinError(i,j,be/bwx/bwy);
		}
	}    
}

double IntegralOfCone( TH2D *hist, double rs, double r0, double r1, double x0, double y0, double &val, double &err, double &area ){
	double rs2 = rs*rs;
	double r02 = r0*r0;
	double r12 = r1*r1;
	double sum = 0;
	double error = 0;
	int nbin=0;
	double dx = hist->GetXaxis()->GetBinWidth(1);
	double dy = hist->GetYaxis()->GetBinWidth(1);
	for( int ix=1;ix <= hist->GetNbinsX();ix++ ){
		double x = hist->GetXaxis()->GetBinCenter(ix)-x0;
		for( int iy=1;iy <= hist->GetNbinsY();iy++ ){
			double y = hist->GetYaxis()->GetBinCenter(iy)*TMath::Pi()-y0*TMath::Pi();
			if(y*y>rs2) continue;
			if( x*x+y*y>r02 && x*x+y*y<r12 ){
				sum += hist->GetBinContent(ix,iy); 
				error += TMath::Power(hist->GetBinError(ix,iy),2); 
				nbin++;
			}
		}
	}
	val = sum;
	err = TMath::Sqrt(error);
	area = nbin*dx*dy*TMath::Pi();
	return sum;
}
// moon for background
double IntegralOfMoon( TH2D *hist, double rs, double rb, double &val, double &err, double &area ){
	double rs2 = rs*rs;
	double rb2 = rb*rb;
	double sum = 0;
	double error = 0;
	int nbin=0;
	double dx = hist->GetXaxis()->GetBinWidth(1); // eta
	double dy = hist->GetYaxis()->GetBinWidth(1); // phi
	for( int ix=1;ix <= hist->GetNbinsX();ix++ ){
		double x = hist->GetXaxis()->GetBinCenter(ix);
		for( int iy=1;iy <= hist->GetNbinsY();iy++ ){
			double y = hist->GetYaxis()->GetBinCenter(iy)*TMath::Pi();
			double flowL = TMath::Sqrt(rb2-y*y); 
			double flowH = flowL + TMath::Sqrt(rs2-y*y); 
			if( TMath::Abs(x)>flowL && TMath::Abs(x)<flowH ){
				sum += hist->GetBinContent(ix,iy); 
				error += TMath::Power(hist->GetBinError(ix,iy),2); 
				nbin++;
			}
		}
	}
	val = sum;
	err = TMath::Sqrt(error);
	area = nbin*dx*dy*TMath::Pi();
	return sum;
}

// moon for background
double IntegralOfSmallR( TH2D *hist, double rs, double rb, double &val, double &err, double &area ){
	double rs2 = rs*rs;
	double rb2 = rb*rb;
	double sum = 0;
	double error = 0;
	int nbin=0;
	double dx = hist->GetXaxis()->GetBinWidth(1); // eta
	double dy = hist->GetYaxis()->GetBinWidth(1); // phi
	if(rs <=  0.3 ) {
		double xcent = rb+rs;
		cout << "rs : rb : xcent = "<< rs <<"\t"<< rb <<"\t"<< xcent <<  endl;
		for( int ix=1;ix <= hist->GetNbinsX();ix++ ){
			double x = hist->GetXaxis()->GetBinCenter(ix);
			for( int iy=1;iy <= hist->GetNbinsY();iy++ ){
				double y = hist->GetYaxis()->GetBinCenter(iy)*TMath::Pi();
				double smallR2 =.2;
				if(x<0) {
					smallR2 = (x+xcent)*(x+xcent) + y*y;
				} else {
					smallR2 = (x-xcent)*(x-xcent) + y*y;
				}
				if( smallR2 < rs2 ){
					sum += hist->GetBinContent(ix,iy); 
					error += TMath::Power(hist->GetBinError(ix,iy),2); 
					nbin++;
				}
			}
		}
	} else {
		for( int ix=1;ix <= hist->GetNbinsX();ix++ ){
			double x = hist->GetXaxis()->GetBinCenter(ix);
			for( int iy=1;iy <= hist->GetNbinsY();iy++ ){
				double y = hist->GetYaxis()->GetBinCenter(iy)*TMath::Pi();
				double flowL = TMath::Sqrt(rb2-y*y); 
				double flowH = flowL + TMath::Sqrt(rs2-y*y); 
				if( TMath::Abs(x)>flowL && TMath::Abs(x)<flowH ){
					sum += hist->GetBinContent(ix,iy); 
					error += TMath::Power(hist->GetBinError(ix,iy),2); 
					nbin++;
				}
			}
		}
	}
	val = sum;
	err = TMath::Sqrt(error);
	area = nbin*dx*dy*TMath::Pi();
	return sum;
}
// for pp
double IntegralForpp( TH2D *hist, double rs, double rb, double &val, double &err, double &area ){
	double rs2 = rs*rs;
	double rb2 = rb*rb;
	double sum = 0;
	double error = 0;
	int nbin=0;
	double dx = hist->GetXaxis()->GetBinWidth(1); // eta
	double dy = hist->GetYaxis()->GetBinWidth(1); // phi
	for( int ix=1;ix <= hist->GetNbinsX();ix++ ){
		double x = hist->GetXaxis()->GetBinCenter(ix);
		for( int iy=1;iy <= hist->GetNbinsY();iy++ ){
			double y = hist->GetYaxis()->GetBinCenter(iy)*TMath::Pi();
			double R2 = x*x + y*y;
			if( R2 > rb2 && R2<1.6*1.6 ){
				sum += hist->GetBinContent(ix,iy); 
				error += TMath::Power(hist->GetBinError(ix,iy),2); 
				nbin++;
			}
		}
	}
	val = sum;
	err = TMath::Sqrt(error);
	area = nbin*dx*dy*TMath::Pi();
	return sum;
}


double IntegralOfRec( TH2D *hist, double rs, double r0, double r1, double x0, double y0, double &val, double &err, double &area ){
	double rs2 = rs*rs;
	double r02 = r0*r0;
	double r12 = r1*r1;
	double sum = 0;
	double error = 0;
	int nbin=0;
	double dx = hist->GetXaxis()->GetBinWidth(1);
	double dy = hist->GetYaxis()->GetBinWidth(1);
	for( int ix=1;ix <= hist->GetNbinsX();ix++ ){
		double x = hist->GetXaxis()->GetBinCenter(ix)-x0;
		for( int iy=1;iy <= hist->GetNbinsY();iy++ ){
			double y = hist->GetYaxis()->GetBinCenter(iy)*TMath::Pi()-y0*TMath::Pi();
			if(y*y>rs2) continue;
			if( x*x>r02 && x*x<r12 ){
				sum += hist->GetBinContent(ix,iy); 
				error += TMath::Power(hist->GetBinError(ix,iy),2); 
				nbin++;
			}
		}
	}
	val = sum;
	err = TMath::Sqrt(error);
	area = nbin*dx*dy*TMath::Pi();
	return sum;
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


Double_t Gaus2D(Double_t *x, Double_t *par) {
	if(par[2] > 0 && par[4] > 0) {
		double rx = ( x[0] - par[1] )/par[2];
		double ry = ( x[1] - par[3] )/par[4];
		return par[0]*TMath::Exp(-(rx*rx+ry*ry)/2.);
	} else {
		return 0.;
	}
}

// Flow subtraction
void SubtracFlowBackground(TH2D* H, TH1D *hflow){ 
	// normalize each bin with the bin width
	int nbx = H->GetNbinsX(); // eta
	int nby = H->GetNbinsY(); // phi

	for(int i=1;i<=nbx;i++){ //loop over all bins 
		for(int j=1;j<=nby;j++){ //loop over all bins 
			double bc = H->GetBinContent(i,j);
			double be = H->GetBinError(i,j);
			double flow = hflow->GetBinContent(j);
			double eflow = hflow->GetBinError(j);
			double jet = bc - flow;
			double ejet = TMath::Sqrt(be*be+eflow*eflow); 
			H->SetBinContent(i,j,jet);
			H->SetBinError(i,j,ejet);
		}    
	}    
}



double GetGeoAccCorrFlat(double deltaEta){
	//FK// calculate acceptance correction on pseudorapidity triangle
	double absDEta = fabs(deltaEta);
	double fmaxEtaRange = 0.8;
	double denominator = 1 - absDEta/(2*fmaxEtaRange);
	//double denominator = 1 - (absDEta - ftriggFiducCut)/(2*fmaxEtaRange-2*ftriggFiducCut);//When Fid>0 max_Deta je mensi nez 2*EtaMax
	//double denominator = 1 - (absDEta - ftriggFiducCut)/(2*fmaxEtaRange-ftriggFiducCut);
	if(denominator > 1e-6)
		return 1.0/denominator;
	else
		return 0;
}
