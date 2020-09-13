const int nFourierTerm = 4;
  
double Fourier(double *x, double *p){
 double f=0;
 f += p[0];
 for(int i=1;i<nFourierTerm;i++){ f += p[i]*cos(x[0]*(double)i); }
 return f;
}


const int nMult = 4;
const int nPtTrig = 4;

int MultMin[nMult] = {0 ,40 ,90 ,110};
int MultMax[nMult] = {40,90,110,1000};
float pTMin[nPtTrig] = { 0.1, 1.0, 2.0, 3.0 };
float pTMax[nPtTrig] = { 1.0, 2.0, 3.0, 4.0 };



void GetPlotFig2(){

 TFile* fin = new TFile("./data/output_set05_grp000_ana0.root","read");

 TH1D* hNTrig[25];
 for(int i=0;i<25;i++) hNTrig[i] = (TH1D*)fin->Get(Form("hNTrig_%d",i));

 TH1D* hNTrigMerged[nMult];
 for(int i=0;i<nMult;i++){
	hNTrigMerged[i] = (TH1D*)hNTrig[0]->Clone();
	hNTrigMerged[i]->Reset();
	for(int j=0;j<25;j++){
		if( j*10 + 5 > MultMin[i] && j*10 + 5 < MultMax[i] ){
			hNTrigMerged[i]->Add( hNTrig[j],1 );
		}
	}
 }


 TH2D* hCorr_Fig2[nMult][nPtTrig];
 TH2D* hCorrMixing_Fig2[nMult][nPtTrig];
 TH2D* hC[nMult][nPtTrig];

 double NTrig[nMult][nPtTrig];

 TH1D* hFPhi[nMult][nPtTrig];
 TF1* fFourier[nMult][nPtTrig];

 TH1D* hFPhiZYAM[nMult][nPtTrig];
 double ZYAMPhi[nMult][nPtTrig];
 double ZYAM[nMult][nPtTrig];

 for(int i=0;i<nMult;i++){
	for(int j=0;j<nPtTrig;j++){
		hCorr_Fig2[i][j] = (TH2D*)fin->Get(Form("hCorr_Fig2_%d_%d",i,j));
		hCorrMixing_Fig2[i][j] = (TH2D*)fin->Get(Form("hCorrMixing_Fig2_%d_%d",i,j));

		hC[i][j] = (TH2D*)hCorr_Fig2[i][j]->Clone();

		hCorrMixing_Fig2[i][j]->Scale( 1.0/hCorrMixing_Fig2[i][j]->GetMaximum() );
		hC[i][j]->Divide( hCorrMixing_Fig2[i][j] );

                hC[i][j]->RebinY(2);
                hC[i][j]->Scale( 1, "width" );

		NTrig[i][j] = 0.0;
		for(int k=0;k<hNTrigMerged[i]->GetNbinsX();k++){
			if( hNTrigMerged[i]->GetBinCenter(k+1) > pTMin[j]  && hNTrigMerged[i]->GetBinCenter(k+1) < pTMax[j] ){
				NTrig[i][j] += hNTrigMerged[i]->GetBinContent(k+1);
			}
		}
	

		hC[i][j]->Scale( 1.0/NTrig[i][j] );

		hC[i][j]->GetXaxis()->SetRangeUser(-4.0,4.0);


		hFPhi[i][j] = (TH1D*)hC[i][j]->ProjectionY(Form("hFPhi_%d_%d",i,j),hC[i][j]->GetXaxis()->FindBin(-4.0),hC[i][j]->GetXaxis()->FindBin(-2.0),"e" );
		hFPhi[i][j]->Add( (TH1D*)hC[i][j]->ProjectionY(Form("hFPhi_%d_%d",i,j),hC[i][j]->GetXaxis()->FindBin(2.0),hC[i][j]->GetXaxis()->FindBin(4.0),"e" ), 1.0 );
		cout << hC[i][j]->GetXaxis()->FindBin(-4.0) << ", " << hC[i][j]->GetXaxis()->FindBin(-2.0) << endl;
		hFPhi[i][j]->Scale( hC[i][j]->GetXaxis()->GetBinWidth(1) );
		hFPhi[i][j]->Scale( 1.0/2.0 );

		fFourier[i][j] = new TF1("f1",Fourier,-10.0,10.0,nFourierTerm);
		for(int n=0;n<nFourierTerm;n++) fFourier[i][j]->SetParLimits(n, -10, 10);

		hFPhiZYAM[i][j] = (TH1D*)hFPhi[i][j]->Clone();
		hFPhi[i][j]->Fit( fFourier[i][j], "" , "" );		

		ZYAMPhi[i][j] = fFourier[i][j]->GetMinimumX(-1.5,1.5);
		ZYAM[i][j] = fFourier[i][j]->GetMinimum(-1.5,1.5);

		for(int k=0;k<hFPhiZYAM[i][j]->GetNbinsX();k++){
			hFPhiZYAM[i][j]->AddBinContent( k+1, -ZYAM[i][j] );
		}
	}
 }

 TFile* fhep_cms = new TFile("/Users/junleekim/RIDGE/Macro_v2/HEPData/CMSResults.root","read");
 TGraphAsymmErrors* ghep_cms[nMult][nPtTrig];
 for(int i=0;i<nMult;i++){
        for(int j=0;j<nPtTrig;j++){
		ghep_cms[i][j] = (TGraphAsymmErrors*)fhep_cms->Get(Form("Table %d/Graph1D_y1",i*8+j*2+1));
	}
 }

 TCanvas* c = new TCanvas("c","c",800,600);
 gStyle->SetOptStat(0);
 gPad->SetLeftMargin(0.13);
 gPad->SetBottomMargin(0.13);
 for(int i=0;i<nMult;i++){
        for(int j=0;j<nPtTrig;j++){
		ghep_cms[i][j]->SetTitle(Form("%s, %.1lf<p_{#font[22]{T}}<%.1lf, %d<N_{#font[22]{ch}}<%d",ghep_cms[i][j]->GetTitle(),pTMin[j],pTMax[j],MultMin[i],MultMax[i]));
		ghep_cms[i][j]->Draw("AP");
		hFPhiZYAM[i][j]->Draw("same");
		c->SaveAs(Form("figs/fig2_dphi_%d_%d.pdf",i,j));
	}
 }
}
