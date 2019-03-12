#include  "JToyFlowHistos.h"


//______________________________________________________________________________
JToyFlowHistos::JToyFlowHistos():
	fUseDirectory(true)
{   // constructor

	
	fTopDirectory = gDirectory;

}

//______________________________________________________________________________
JToyFlowHistos::JToyFlowHistos(const JToyFlowHistos& obj){
	// copy constructor
}

//______________________________________________________________________________
JToyFlowHistos& JToyFlowHistos::operator=(const JToyFlowHistos& obj){
	// copy constructor
	return *this;
}

//______________________________________________________________________________
void JToyFlowHistos::CreateHistos()
{
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
}






















