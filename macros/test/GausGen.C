TH1D *Flip(TH1D* hin, int i);

void GausGen() {

	// create a random number generator
	gRandom = new TRandom3();
	// create a hgausogram 
	TH1D * hgaus = new TH1D("hgaus", "Gauss", 100, -1., 1.0);
	// fill in the hgausogram
	for (int i = 0; i < 10000; ++i)
		hgaus->Fill(gRandom->Gaus(0, .2));

	TCanvas * c1= new TCanvas("c1", "random",5,5,800,600);
	//hgaus->Scale(1.,"width");
	hgaus->Draw("e");


	TH1D *hgausflip = Flip((TH1D*)hgaus, 1);
	hgausflip->SetLineColor(2);
	hgausflip->Draw("same");
}

//------------------------------------------------------------------------------------------------
TH1D *Flip(TH1D* hin, int idtyp){
	int nb  = hin->GetNbinsX();
	double scale = 2.0; // MUST Have dN/deta -> dN/d|eta|
	double max = hin->GetBinLowEdge(nb+1);
	TString hname = hin->GetName();
	TString newName = Form("%s_flip%d",hname.Data(),idtyp);

	TH1D *hFlip = new TH1D(newName.Data(), newName.Data(), (int) nb/2., 0, max);
	hFlip->SetTitle(hin->GetTitle());
	int zero = hin->FindBin(0.00001);
	for(int ib=zero; ib<=nb; ib++){
		double valPos = hin->GetBinContent(ib);
		double errPos = hin->GetBinError(ib);
		double valNeg = hin->GetBinContent(nb - ib+1);
		double errNeg = hin->GetBinError(nb - ib+1);

		hFlip->SetBinContent(ib-zero+1, (valPos+valNeg)/scale);
		hFlip->SetBinError(ib-zero+1, sqrt( errPos*errPos + errNeg * errNeg )/scale);
	}

	return hFlip;
}
