
#include <stdio.h>
#include <vector>
#include <TRandom.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TComplex.h>

typedef unsigned int uint;

//#define NH 3
#define NC 6

//https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/jonderwa/2014-Apr-29-analysis_note-analysis_note_event_plane_calibration.pdf
//https://www.hepdata.net/record/78365
//https://indico.cern.ch/event/703569/contributions/2886293/attachments/1597688/2531598/2018.02.08-Slupecki-ToyFlow.pdf
//https://aliceinfo.cern.ch/ArtSubmission/sites/aliceinfo.cern.ch.ArtSubmission/files/draft/cholm/2017-Jun-15-paper_draft-cds-rb-adraft-20170109-1.pdf

//#define ETADST_N 34
#define ETADST_N 36
static double etadst[ETADST_N] = {
	-3.8, //underflow
	-3.375,-3.125,-2.875,-2.625,-2.375,-2.125,-1.875,-1.625,
    -1.375,-1.125,-0.875,-0.625,-0.375,-0.125, 0.125, 0.375,
     0.625, 0.875, 1.125, 1.375, 1.625, 1.875, 2.125, 2.375,
     2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.125, 4.375,
     4.625, 4.875,
	 5.2, //overflow
};

static double etanch[NC+2][ETADST_N] = {
	{1640,1643,1670,1718,1787,1835,1912,1968,2001,2021,2017,1995,1970,1943,1929,1929,1943,1970,1995,2017,2021,2001,1968,1912,1835,1787,1718,1670,1643,1563,1474,1370,1324,1281,1244,1240},
	{1360,1364,1391,1424,1474,1507,1569,1644,1679,1682,1672,1646,1621,1597,1583,1583,1597,1621,1646,1672,1682,1679,1644,1569,1507,1474,1424,1391,1364,1292,1218,1132,1093,1062,1032,1030},
	{1000,1038,1061,1080,1114,1136,1178,1229,1253,1256,1247,1229,1210,1191,1181,1181,1191,1210,1229,1247,1256,1253,1229,1178,1136,1114,1080,1061,1038,977,921.3,857.7,829.6,807.4,787,780},
	{700,714,726,738,759,772,797,827,842,844,838,826,811.9,799.2,792.4,792.4,799.2,811.9,826,838,844,842,827,797,772,759,738,726,714,665,625.4,582.6,565.5,551.4,538,530},
	{460,475,482.7,489.7,502.6,510.6,522,539.9,549,549.3,545.5,537.5,527.6,519.3,514.7,514.7,519.3,527.6,537.5,545.5,549.3,549,539.9,522,510.6,502.6,489.7,482.7,475,440,413.6,386.7,375.6,368,359.9,350},
	{300,302,306.3,310.1,317.9,322.3,327.6,335.1,340,340.2,337.7,332.5,326.3,320.7,317.5,317.5,320.7,326.3,332.5,337.7,340.2,340,335.1,327.6,322.3,317.9,310.1,306.3,302,277.5,261.3,244.7,238.4,233.8,229.4,225},
	{177,178,179.9,181.7,186,188.2,189.8,193.5,196.4,196.5,194.8,191.4,187.5,184.3,182.5,182.5,184.3,187.5,191.4,194.8,196.5,196.4,193.5,189.8,188.2,186,181.7,179.9,178,163.2,153.4,143.8,140.3,138.7,136,135},
	{93,94.9,96.1,96.8,98.3,98.8,99.1,101.2,102.7,103.1,102,100.3,98,96.1,95.2,95.2,96.1,98,100.3,102,103.1,102.7,101.2,99.1,98.8,98.3,96.8,96.1,94.9,86.8,81.9,77.3,75.8,75.1,73.8,73}
};

static double CentBins[NC+3] = {0,5,10,20,30,40,50, 60, 70};

enum DETECTOR{
	D_TPC, //TPC full coverage
	D_TPC_ETAA, //TPC with eta gap
	D_TPC_ETAC,
	D_V0A,
	D_V0C,
	D_V0P, //V0+
	D_COUNT
};
static double cov[D_COUNT][2] = {
	{-1.5,1.5},
	{-1.5,-0.4},
	{0.4,1.5},
	{2.8,5.1},
	{-3.7,-1.7},
	{2.19,5.08}
};

enum RESOLUTION{
	R_V0A,
	R_V0C,
	R_V0P,
	R_COUNT
};
static const char *presn[] = {"V0A","V0C","V0P"};

int checkplane(double evp, double phi) {
  double pi = TMath::Pi();
  double diff = phi- evp;
  if (diff < 0) diff += 2*pi;

//  if (diff < 1./6*pi && diff >= 0) return 0;
//  else if (diff < 2./6*pi && diff >= 1./6*pi ) return 1;
//  else if (diff < 4./6*pi && diff >= 2./6*pi ) return 2;
//  else if (diff < 5./6*pi && diff >= 4./6*pi ) return 3;
//  else if (diff < 7./6*pi && diff >= 5./6*pi ) return 4;
//  else if (diff < 8./6*pi && diff >= 7./6*pi ) return 5;
//  else if (diff < 10./6*pi && diff >= 8./6*pi ) return 6;
//  else if (diff < 11./6*pi && diff >= 10./6*pi ) return 7;
//  else if (diff < 12./6*pi && diff >= 11./6*pi ) return 0;

//  if (diff < 1./8*pi && diff >= 0) return 0;
//  else if (diff < 3./8*pi && diff >= 1./8*pi ) return 1;
//  else if (diff < 5./8*pi && diff >= 3./8*pi ) return 2;
//  else if (diff < 7./8*pi && diff >= 5./8*pi ) return 3;
//  else if (diff < 9./8*pi && diff >= 7./8*pi ) return 4;
//  else if (diff < 11./8*pi && diff >= 9./8*pi ) return 5;
//  else if (diff < 13./8*pi && diff >= 11./8*pi ) return 6;
//  else if (diff < 15./8*pi && diff >= 13./8*pi ) return 7;
//  else if (diff < 16./8*pi && diff >= 15./8*pi ) return 0;

  if (diff < 1./12*pi && diff >= 0) return 0;
  else if (diff < 4./12*pi && diff >= 2./12*pi ) return 1;
  else if (diff < 7./12*pi && diff >= 5./12*pi ) return 2;
  else if (diff < 10./12*pi && diff >= 8./12*pi ) return 3;
  else if (diff < 13./12*pi && diff >= 11./12*pi ) return 4;
  else if (diff < 16./12*pi && diff >= 14./12*pi ) return 5;
  else if (diff < 19./12*pi && diff >= 17./12*pi ) return 6;
  else if (diff < 22./12*pi && diff >= 20./12*pi ) return 7;
  else if (diff < 24./12*pi && diff >= 23./12*pi ) return 0;

  return -9;
}

double CheckDetectorPhi(double phi) {

  double pi = TMath::Pi();
  double angle[8] = {-3*pi/4, -1*pi/2, -1*pi/4, 0, pi/4, pi/2, 3*pi/4, pi};
  double medianangle[8] = {-7*pi/8, -5*pi/8, -3*pi/8, -1*pi/8, pi/8, 3*pi/8, 5*pi/8, 7*pi/8};
  int i = 0;
  for ( ; i < 8; i++) {
    if (phi < angle[i]) break;
  }
  return medianangle[i];
}



//void resolution(){
int main(int argc, char **pargv){
	uint seed = argc > 1?atol(pargv[1]):1000;
	uint evtc = argc > 2?atol(pargv[2]):1000;
	printf("seed:\t%u\nevents:\t%u\n",seed,evtc);
	
	//Read flow coefficients ----------------------------
	const char *pgn[3] = {"Hist1D_y1","Hist1D_y1_e1","Hist1D_y1_e2"};
	const char *pglabel[4] = {"Hist1D_y%u","Hist1D_y%u_e1","Hist1D_y%u_e2plus", "Hist1D_y%u_e2minus"};
	TFile *pff = new TFile("anizo-flow.root","read");

  TFile *highf = new TFile("276_HighPtFlow.root");

#define N_VN 2
	TGraphErrors *pgr_v[N_VN];
	TGraphErrors *high_v[N_VN+1];

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

  for (uint i = 0; i < N_VN+1; ++i) {
      high_v[i] =  new TGraphErrors();
    for (uint c = 0; c < NC; ++c) {
      TH1D *temp_hist[5];
      for (int ilbl = 0; ilbl < 3; ilbl++) {
        temp_hist[ilbl] = (TH1D*)highf->Get(Form("Table %u/%s", ilbl+1, Form(pglabel[ilbl], c+1)));
      }  // y1 : v2, y2 : v2{4}, y3 : v3;
      
      high_v[i]->SetPoint(high_v[i]->GetN(), (CentBins[c] + CentBins[c+1])/2, temp_hist[0]->GetBinContent( temp_hist[0]->FindBin(8.1)) );
      high_v[i]->SetPointError(high_v[i]->GetN() - 1 , 0, temp_hist[1]->FindBin(8.1));
    }
  }

	pff->Close();
	delete pff;
	highf->Close();
	delete highf;

	double pi = TMath::Pi();
	TF1 *pdf = new TF1("df",
		"[0]*(1+2*[1]*cos(x-[5])+2*[2]*cos(2*(x-[6]))+2*[3]*cos(3*(x-[7]))+2*[4]*cos(4*(x-[8])))",-pi,pi);
	pdf->SetParameter(0,100.0); // v0
	pdf->SetParameter(1,0.0); // v1
	TF1 *pdf_high = new TF1("dfhigh",
		"[0]*(1+2*[1]*cos(x-[5])+2*[2]*cos(2*(x-[6]))+2*[3]*cos(3*(x-[7]))+2*[4]*cos(4*(x-[8])))",-pi,pi);
	pdf_high->SetParameter(0,100.0); // v0
	pdf_high->SetParameter(1,0.0); // v1

	TRandom *prng = new TRandom(seed);

	//correlation histograms
	TH1D *pah[R_COUNT][NC], *pbh[R_COUNT][NC], *pch[R_COUNT][NC];
  TH1D *evph[R_COUNT][NC];
  TH2D *contami2d[R_COUNT][NC];
  TH2D *highcontami2d[R_COUNT][NC];
  TH2D *evpcorr2d[R_COUNT][NC];
  TH2D *evpcorrvsdet2d[NC];
  TH1D *evpdifference[R_COUNT][NC];
	for(uint i = 0; i < R_COUNT; ++i){
		for(uint j = 0; j < NC; ++j){
			pah[i][j] = new TH1D(Form("h_res_%s_a%02u",presn[i],j),"h_res",1024,-1.5,1.5);
			pbh[i][j] = new TH1D(Form("h_res_%s_b%02u",presn[i],j),"h_res",1024,-1.5,1.5);
			pch[i][j] = new TH1D(Form("h_res_%s_c%02u",presn[i],j),"h_res",1024,-1.5,1.5);
      contami2d[i][j] = new TH2D(Form("h_contami2d_%s_c%02u", presn[i], j), Form("h_contami2d_%s_%.0f-%.0fC", presn[i],CentBins[j], CentBins[j+1] ), 8, -0.5, 7.5, 8, -0.5, 7.5);
      highcontami2d[i][j] = new TH2D(Form("h_highcontami2d_%s_c%02u", presn[i], j), Form("h_highcontami2d_%s_%.0f-%.0fC", presn[i],CentBins[j], CentBins[j+1] ), 8, -0.5, 7.5, 8, -0.5, 7.5);
      evph[i][j] = new TH1D(Form("h_evp_%s_%02u", presn[i], j), "h_evp", 64, -1.6, 1.6);
      evpcorr2d[i][j] = new TH2D(Form("h_evpcorr2d_%s_c%02u", presn[i], j), Form("h_evpcorr2d_%s_%.0f-%.0fC", presn[i],CentBins[j], CentBins[j+1] ), 128, -4, 4, 128, -4, 4);
      evpcorr2d[i][j]->GetXaxis()->SetTitle("Evp from True");
      evpcorr2d[i][j]->GetYaxis()->SetTitle("Evp from Reco");
      
      evpdifference[i][j] = new TH1D(Form("h_evpdiff_%s_c%02u", presn[i], j), Form("Event Plane Difference %s, %.0f-%.0fC", presn[i], CentBins[j], CentBins[j+1]), 128, -4, 4); 
      evpdifference[i][j]->GetXaxis()->SetTitle("Event plane Difference");
      evpdifference[i][j]->GetYaxis()->SetTitle("N Events");

      if (i == 0) {
        evpcorrvsdet2d[j] = new TH2D(Form("h_evpcorrdet2d_c%02u", j), Form("h_evpcorr2d_%.0f-%.0fC", CentBins[j], CentBins[j+1] ), 64, -1.6, 1.6, 64, -1.6, 1.6);
        evpcorrvsdet2d[j]->GetXaxis()->SetTitle("Evp from V0A");
        evpcorrvsdet2d[j]->GetYaxis()->SetTitle("Evp from V0C");
      }
      
		}
	}

	//Create multiplicity(eta) distribution & integrate ---------------
	/*double etadst[ETADST_N];
	for(uint i = 0; i < ETADST_N; ++i)
		etadst[i] = -3.5+(double)i*0.25+0.125;*/
	
	TGraph *pgr_nch[D_COUNT];
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

//#define MAX_N 10000
	//double eta[MAX_N], phi[MAX_N];
	//double phi[MAX_N];

	for(uint evt = 0; evt < evtc; ++evt){
		//Event generation ----------------------------

		double cent = prng->Uniform(0,50.0);
		pdf->SetParameter(2,pgr_v[0]->Eval(cent) * 0.9); //v2
		pdf->SetParameter(3,pgr_v[1]->Eval(cent)); //v3
		pdf->SetParameter(4,0.01);//pgr_v[2]->Eval(cent)); //v4

		pdf->SetParameter(5,prng->Uniform(-pi,pi));  //EP for v1
		pdf->SetParameter(6,prng->Uniform(-pi/2.0,pi/2.0)); //EP for v2
		pdf->SetParameter(7,prng->Uniform(-pi/3.0,pi/3.0)); //EP for v3
		pdf->SetParameter(8,prng->Uniform(-pi/4.0,pi/4.0)); //EP for v4

		pdf_high->SetParameter(2,high_v[0]->Eval(cent) * 0.9); //v2
		pdf_high->SetParameter(3,high_v[2]->Eval(cent)); //v3
		pdf_high->SetParameter(4,0.01);//pgr_v[2]->Eval(cent)); //v4

		pdf_high->SetParameter(5,pdf->GetParameter(5));  //EP for v1
		pdf_high->SetParameter(6,pdf->GetParameter(6)); //EP for v2
		pdf_high->SetParameter(7,pdf->GetParameter(7)); //EP for v3
		pdf_high->SetParameter(8,pdf->GetParameter(8)); //EP for v4


		uint cid = 0;
		for(; cid < NC; ++cid)
			if(cent < CentBins[cid+1])
				break;

		//Q-vectors -----------------------------------
    double trueevp = pdf->GetParameter(6);

    std::vector <double> trackphi[D_COUNT];
		TComplex Qsd[D_COUNT];
		for(uint s = 0; s < D_COUNT; ++s){
			uint ntracks = (uint)pgr_nch[s]->Eval(cent) * 0.9;

			TComplex Qa2 = TComplex(0,0);
			for(uint i = 0; i < ntracks; ++i){
				double tphi = pdf->GetRandom();
//				phi = TMath::Floor(8.0*(phi+pi)/(2.0*pi))*(2.0*pi)/8-pi;
        double phi = CheckDetectorPhi(tphi);
        trackphi[s].push_back(phi);
				Qa2 += TComplex(TMath::Cos(2.0*phi),TMath::Sin(2.0*phi));
			}
      double thphi = pdf_high->GetRandom();
      double hphi = CheckDetectorPhi(thphi);
      trackphi[s].push_back(hphi);
			Qa2 += TComplex(TMath::Cos(2.0*hphi),TMath::Sin(2.0*hphi));
			
			Qa2 /= (double) (ntracks+1);
			Qsd[s] = Qa2/TComplex::Abs(Qa2);
		}

    //Calculate Evp

    for (uint s = 0; s < R_COUNT; ++s) {
      double recoevp = TMath::ATan2(Qsd[s].Im(), Qsd[s].Re())/2;
      double evpdiff = trueevp - recoevp;
      
      evph[s][cid]->Fill(recoevp);
      evpcorr2d[s][cid]->Fill(trueevp, recoevp);
      evpdifference[s][cid]->Fill(evpdiff);
      
      for (uint i = 0; i < trackphi[s].size(); ++i) {
        
//        contami2d[s][cid]->Fill(checkplane(trueevp, trackphi[i]), 0);
        contami2d[s][cid]->Fill(checkplane(trueevp, trackphi[s][i]), checkplane(recoevp, trackphi[s][i]));
      }  // For Full contamination
      highcontami2d[s][cid]->Fill(checkplane(trueevp, trackphi[s].back()), checkplane(recoevp, trackphi[s].back()));
    }

//    evpcorrvsdet2d[cid]->Fill(TMath::ATan2(Qsd[0].Im(), Qsd[0].Re()), TMath::ATan2(Qsd[1].Im(), Qsd[1].Re()));


		//Calculate the resolution components
		TComplex ab, ac, bc;
		/*ab = Qsd[D_V0A]*TComplex::Conjugate(Qsd[D_TPC]);
		ac = Qsd[D_V0A]*TComplex::Conjugate(Qsd[D_V0C]);
		bc = Qsd[D_TPC]*TComplex::Conjugate(Qsd[D_V0C]);
		pah[R_V0A][cid]->Fill(ab.Re());
		pbh[R_V0A][cid]->Fill(ac.Re());
		pch[R_V0A][cid]->Fill(bc.Re());

		ab = Qsd[D_V0C]*TComplex::Conjugate(Qsd[D_V0A]);
		ac = Qsd[D_V0C]*TComplex::Conjugate(Qsd[D_TPC]);
		bc = Qsd[D_V0A]*TComplex::Conjugate(Qsd[D_TPC]);
		pah[R_V0C][cid]->Fill(ab.Re());
		pbh[R_V0C][cid]->Fill(ac.Re());
		pch[R_V0C][cid]->Fill(bc.Re());*/
		ab = Qsd[D_V0A]*TComplex::Conjugate(Qsd[D_TPC_ETAA]);
		ac = Qsd[D_V0A]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
		bc = Qsd[D_TPC_ETAA]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
		pah[R_V0A][cid]->Fill(ab.Re());
		pbh[R_V0A][cid]->Fill(ac.Re());
		pch[R_V0A][cid]->Fill(bc.Re());

		ab = Qsd[D_V0C]*TComplex::Conjugate(Qsd[D_TPC_ETAA]);
		ac = Qsd[D_V0C]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
		//bc = Qsd[D_TPC_ETAA]*TComplex::Conjugate(Qsd[D_TPC]);
		pah[R_V0C][cid]->Fill(ab.Re());
		pbh[R_V0C][cid]->Fill(ac.Re());
		pch[R_V0C][cid]->Fill(bc.Re());

		ab = Qsd[D_V0P]*TComplex::Conjugate(Qsd[D_TPC_ETAA]);
		ac = Qsd[D_V0P]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
		bc = Qsd[D_TPC_ETAA]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
		pah[R_V0P][cid]->Fill(ab.Re());
		pbh[R_V0P][cid]->Fill(ac.Re());
		pch[R_V0P][cid]->Fill(bc.Re());

	}

	for(uint i = 0; i < D_COUNT; ++i)
		delete pgr_nch[i];
	for(uint i = 0; i < N_VN; ++i) {
		delete pgr_v[i];
    delete high_v[i];
  }

	delete prng;
	delete pdf;
  delete pdf_high;

	TFile *pfo = new TFile(argc > 3?pargv[3]:"results.root","recreate");
	pfo->cd();


	for(uint i = 0; i < R_COUNT; ++i){
		for(uint j = 0; j < NC; ++j){
			//double a[2] = {pah[i]->GetMean(),pah[i]->GetMeanError()};
			//double b[2] = {pbh[i]->GetMean(),pah[i]->GetMeanError()};
			//double c[2] = {pch[i]->GetMean(),pah[i]->GetMeanError()};
			//
			//double R = TMath::Sqrt(a[0]*b[0]/c[0]);
			//double e = TMath::Sqrt(b[0]*a[1]*a[1]/(a[0]*c[0])
			//	+a[0]*b[1]*b[1]/(b[0]*c[0])+a[0]*b[0]*c[1]/(c[0]*c[0]*c[0]));
			//printf("R2(%u) = %lf pm %lf\n",i,R,e);

			pah[i][j]->Write(Form("h_%s_a%02u",presn[i],j));
			pbh[i][j]->Write(Form("h_%s_b%02u",presn[i],j));
			pch[i][j]->Write(Form("h_%s_c%02u",presn[i],j));

      contami2d[i][j]->GetXaxis()->SetTitle("Evp from True");
      contami2d[i][j]->GetYaxis()->SetTitle("Evp from Reco");
      contami2d[i][j]->Write(Form("h_contami2d_%s_%02u", presn[i], j));
      highcontami2d[i][j]->GetXaxis()->SetTitle("Evp from True");
      highcontami2d[i][j]->GetYaxis()->SetTitle("Evp from Reco");
      highcontami2d[i][j]->Write(Form("h_highcontami2d_%s_%02u", presn[i], j));

      evpcorr2d[i][j]->Write(Form("h_evpcorr2d_%s_%02u", presn[i], j));
      evph[i][j]->Write(Form("h_evp_%s_%02u", presn[i], j));
      evpdifference[i][j]->Write(Form("h_evpdiff_%s_%02u", presn[i], j));
 
      if (i == 0) {
      evpcorrvsdet2d[j]->Write(Form("h_evpcorrvsdet2d_%02u", j));
      delete evpcorrvsdet2d[j];
      }

			delete pah[i][j];
			delete pbh[i][j];
			delete pch[i][j];
      delete contami2d[i][j];
      delete highcontami2d[i][j];
      delete evpcorr2d[i][j];
      delete evph[i][j];
      delete evpdifference[i][j];
    }
	}

	pfo->Close();
	delete pfo;

	return 0;
}

