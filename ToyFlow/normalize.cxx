
#define NC 6             //Changed for high pT test
static const char *presn[] = {"V0A","V0C","V0P"};

enum RESOLUTION{
	R_V0A,
	R_V0C,
	R_V0P,
	R_COUNT
};

static double CentBins[NC+3] = {0,5,10,20,30,40,50,60,70};
normalize() {

  double pi = TMath::Pi();
  string name;
  TFile *f1 = new TFile("lowerv2_35.root");                     // Resolution & smearing histos
  TFile *f2 = new TFile("normal_lowerv2_35.root", "recreate");  // Output after normalization
 // TFile *f1 = new TFile("lowerv2_45_true.root");                     // Resolution & smearing histos
 // TFile *f2 = new TFile("normal_lowerv2_45_true.root", "recreate");  // Output after normalization
  TH2D *rawfile[R_COUNT][NC];

  TH1D *pah[R_COUNT][NC], *pbh[R_COUNT][NC], *pch[R_COUNT][NC];

  TH1D *evpcorr_1d[R_COUNT][NC][8];

  TH1D *resolution[R_COUNT];

  TH2D *evpcorr2d[R_COUNT][NC];

  TH2D *normalized[R_COUNT][NC];
  TH1D *evpdifference[R_COUNT][NC];
  TH2D *rawfilehigh[R_COUNT][NC];
  TH2D *evpcorr2dhigh[R_COUNT][NC];
  TH2D *normalizedhigh[R_COUNT][NC];
  for (int s = 0; s < R_COUNT; s++) {
    resolution[s] = new TH1D(Form("resolution_%s", presn[s]), Form("resolution_%s", presn[s]), NC, CentBins);
    for (int c = 0; c < NC; c++) {
      rawfile[s][c] = (TH2D*)f1->Get(Form("h_contami2d_%s_%02u", presn[s], c));
      normalized[s][c] = new TH2D(Form("h_contami2d_normal_%s_%02u", presn[s], c), Form("%s, Scaled with 1/Reconstructed in same Reconstructed Sector, %.0f-%.0fC", presn[s], CentBins[c], CentBins[c+1]), rawfile[s][c]->GetNbinsX(), rawfile[s][c]->GetXaxis()->GetBinLowEdge(1), rawfile[s][c]->GetXaxis()->GetBinLowEdge(rawfile[s][c]->GetNbinsX() + 1), rawfile[s][c]->GetNbinsY(), rawfile[s][c]->GetYaxis()->GetBinLowEdge(1), rawfile[s][c]->GetYaxis()->GetBinLowEdge(rawfile[s][c]->GetNbinsY() + 1) );
      normalized[s][c]->GetXaxis()->SetTitle("Sector from true Event-plane");
      normalized[s][c]->GetYaxis()->SetTitle(Form("Sector from Reco Event-plane %s", presn[s]));
			pah[s][c]=(TH1D*)f1->Get(Form("h_%s_a%02u",presn[s],c));
			pbh[s][c]=(TH1D*)f1->Get(Form("h_%s_b%02u",presn[s],c));
			pch[s][c]=(TH1D*)f1->Get(Form("h_%s_c%02u",presn[s],c));

      evpdifference[s][c] = (TH1D*) f1->Get(Form("h_evpdiff_%s_%02d", presn[s], c));
      evpdifference[s][c]->GetXaxis()->SetTitle("Event Plane Difference #Psi_{true}-#Psi{reco}");
      evpdifference[s][c]->GetYaxis()->SetTitleOffset(1.3);
      evpdifference[s][c]->SetTitle("");
      
      evpcorr2d[s][c] = (TH2D*)f1->Get(Form("h_evpcorr2d_%s_%02u", presn[s], c));
      resolution[s]->SetBinContent(c+1, TMath::Sqrt(pah[s][c]->GetMean()*pbh[s][c]->GetMean()/pch[s][c]->GetMean()));
      evpcorr2d[s][c]->GetXaxis()->SetRangeUser(-0.5*pi, 0.5*pi);
      evpcorr2d[s][c]->GetYaxis()->SetRangeUser(-0.5*pi, 0.5*pi);

      for (int j = 0; j < rawfile[s][c]->GetNbinsY() + 1; j++) {
        double normalizationFactor  = 0;
        for (int i = 0; i < rawfile[s][c]->GetNbinsX() + 1; i++) {
          normalizationFactor += rawfile[s][c]->GetBinContent(i, j);
        }
        for (int i = 0; i < rawfile[s][c]->GetNbinsX() + 1; i++) {

          normalized[s][c]->SetBinContent(i, j,rawfile[s][c]->GetBinContent(i, j)/normalizationFactor );   // Normalization process Error is not added yet
        }
      }

      for (int ibin = 0; ibin < 8; ibin++) {
        evpcorr_1d[s][c][ibin] = normalized[s][c]->ProjectionX(Form("h_contami1d_%s_%02u%02u", presn[s], c, ibin), ibin+1, ibin+1);
        evpcorr_1d[s][c][ibin]->SetTitle(Form("Event plane Reconstructed to Sector %d, %.0f-%.0fC", ibin, CentBins[c], CentBins[c+1]));
      }



      rawfilehigh[s][c] = (TH2D*)f1->Get(Form("h_highcontami2d_%s_%02u", presn[s], c));
      normalizedhigh[s][c] = new TH2D(Form("h_highcontami2d_normal_%s_%02u", presn[s], c), Form("%s, Scaled with 1/Reconstructed in same Reconstructed Sector, %.0f-%.0fC", presn[s], CentBins[c], CentBins[c+1]), rawfilehigh[s][c]->GetNbinsX(), rawfilehigh[s][c]->GetXaxis()->GetBinLowEdge(1), rawfilehigh[s][c]->GetXaxis()->GetBinLowEdge(rawfilehigh[s][c]->GetNbinsX() + 1), rawfilehigh[s][c]->GetNbinsY(), rawfilehigh[s][c]->GetYaxis()->GetBinLowEdge(1), rawfilehigh[s][c]->GetYaxis()->GetBinLowEdge(rawfilehigh[s][c]->GetNbinsY() + 1) );
      normalizedhigh[s][c]->GetXaxis()->SetTitle("Sector from true Event-plane");
      normalizedhigh[s][c]->GetYaxis()->SetTitle(Form("Sector from Reco Event-plane %s", presn[s]));
      
      for (int j = 0; j < rawfilehigh[s][c]->GetNbinsY() + 1; j++) {
        double normalizationFactor  = 0;
        for (int i = 0; i < rawfilehigh[s][c]->GetNbinsX() + 1; i++) {
          normalizationFactor += rawfilehigh[s][c]->GetBinContent(i, j);
        }
        for (int i = 0; i < rawfilehigh[s][c]->GetNbinsX() + 1; i++) {

          normalizedhigh[s][c]->SetBinContent(i, j,rawfilehigh[s][c]->GetBinContent(i, j)/normalizationFactor );
        }
      }
    

    }
  }

  f2->cd();
  TCanvas *c1 = new TCanvas();
  /* Drawing and Writing to file */
  for (int s = 0; s < R_COUNT; s++) {
    resolution[s]->Write();

    resolution[s]->Draw();
    c1->SaveAs(Form("plots/Resolution%s.pdf", presn[s]));
    for (int c = 0; c < NC; c++) {
      normalized[s][c]->SetStats(0);
      normalized[s][c]->Scale(100.);
      normalized[s][c]->Write();
      normalized[s][c]->Draw("text");
      c1->SaveAs(Form("plots/%s_%d.pdf", presn[s], c));
      
      normalizedhigh[s][c]->SetStats(0);
      normalizedhigh[s][c]->Scale(100.);
      normalizedhigh[s][c]->Write();
      normalizedhigh[s][c]->Draw("colz");
      c1->SaveAs(Form("plots/high%s_%d.pdf", presn[s], c));

      evpdifference[s][c]->SetStats(0);
      
        TLatex l2;
        l2.SetTextSize(0.04);
      evpdifference[s][c]->Write();
      evpdifference[s][c]->SetMaximum(evpdifference[s][c]->GetMaximum() * 1.2);
      evpdifference[s][c]->Draw();
        l2.DrawLatexNDC(0.12, 0.85, Form("Toy MC, #sqrt{s_{NN}}=5.02TeV, Centrality %.0f-%.0f%%", CentBins[c], CentBins[c+1]));
        l2.DrawLatexNDC(0.12, 0.8, Form("Event plane reconstructed with V0C Detector"));
      c1->SaveAs(Form("plots/evpdifference%s_%d.pdf", presn[s], c));

      evpcorr2d[s][c]->SetStats(0);
      evpcorr2d[s][c]->Draw("colz");
      c1->SaveAs(Form("plots/%s_evpcorr2d_%d.pdf", presn[s], c));

      for (int ibin = 0; ibin < 8; ibin++) {
        evpcorr_1d[s][c][ibin]->SetStats(0);
        evpcorr_1d[s][c][ibin]->Scale(100);
        evpcorr_1d[s][c][ibin]->GetYaxis()->SetTitle("Fraction of reconstructed from true plane (%)");
        evpcorr_1d[s][c][ibin]->GetYaxis()->SetRangeUser(0, 70);
        evpcorr_1d[s][c][ibin]->SetTitle("");
        TLatex l1;
        evpcorr_1d[s][c][ibin]->Draw();
        l1.SetTextSize(0.04);
        l1.DrawLatex(0, 65, Form("Toy MC, #sqrt{s_{NN}}=5.02TeV, Centrality %.0f-%.0f%%", CentBins[c], CentBins[c+1]));
        l1.DrawLatex(0, 61.5, Form("Event plane reconstructed with V0C Detector"));
        l1.DrawLatex(0, 58, Form("Reconstructed to Sector %d", ibin));
        c1->SaveAs(Form("plots/%s_evpcorr1d_%d%d.pdf", presn[s], c, ibin+1));
      }
    }
  }

}
