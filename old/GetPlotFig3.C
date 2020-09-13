const int nFourierTerm = 4;
  
double Fourier(double *x, double *p){
 double f=0;
 f += p[0];
 for(int i=1;i<nFourierTerm;i++){ f += p[i]*cos(x[0]*(double)i); }
 return f;
}


const int nMult = 1;
const int nPtTrig = 8;

float pTMin[nPtTrig] = {0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0 };
float pTMax[nPtTrig] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0 };
int MultMin[nMult] = {110};
int MultMax[nMult] = {1000};


void GetPlotFig3(){

 TFile* fin = new TFile("output_set05_grp000_ana1.root","read");

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


 TH2D* hCorr_Fig3[nMult][nPtTrig];
 TH2D* hCorrMixing_Fig3[nMult][nPtTrig];
 TH2D* hC[nMult][nPtTrig];

 double NTrig[nMult][nPtTrig];

 TH1D* hFPhi[nMult][nPtTrig];
 TF1* fFourier[nMult][nPtTrig];

 TH1D* hFPhiZYAM[nMult][nPtTrig];
 double ZYAMPhi[nMult][nPtTrig];
 double ZYAM[nMult][nPtTrig];

 double YieldCntl[nMult][nPtTrig];
 double YieldStat[nMult][nPtTrig];

 for(int i=0;i<nMult;i++){
	for(int j=0;j<nPtTrig;j++){
		hCorr_Fig3[i][j] = (TH2D*)fin->Get(Form("hCorr_Fig3_%d",j));
		hCorrMixing_Fig3[i][j] = (TH2D*)fin->Get(Form("hCorrMixing_Fig3_%d",j));

		hC[i][j] = (TH2D*)hCorr_Fig3[i][j]->Clone();

		hCorrMixing_Fig3[i][j]->Scale( 1.0/hCorrMixing_Fig3[i][j]->GetMaximum() );
		hC[i][j]->Divide( hCorrMixing_Fig3[i][j] );

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

		YieldCntl[i][j] = 0.0;
		YieldStat[i][j] = 0.0;
		for(int k=0;k<hFPhiZYAM[i][j]->GetNbinsX();k++){
			if( fabs( hFPhiZYAM[i][j]->GetBinCenter(k+1) ) > fabs( ZYAMPhi[i][j] ) ) continue;
			YieldCntl[i][j] += hFPhiZYAM[i][j]->GetBinContent(k+1) * hFPhiZYAM[i][j]->GetBinWidth(k+1);
			YieldStat[i][j] += pow( hFPhiZYAM[i][j]->GetBinError(k+1)*hFPhiZYAM[i][j]->GetBinWidth(k+1),2 );
		}
		YieldStat[i][j] = sqrt( YieldStat[i][j] );
	}
 }

 double ptcntl[nPtTrig];
 double ptstat[nPtTrig];

 for(int i=0;i<nPtTrig;i++){
	ptcntl[i] = (pTMin[i]+pTMax[i])/2.0;
	ptstat[i] = (pTMax[i]-pTMin[i])/2.0;
 }

 TGraphErrors* g = new TGraphErrors(nPtTrig,ptcntl,YieldCntl[0],ptstat,YieldStat[0]);
 g->SetLineColor(kRed);

 TFile* fhep_cms = new TFile("/Users/junleekim/RIDGE/Macro_v2/HEPData/CMSResults.root","read"); //33, 35
 TGraphAsymmErrors* ghep_cms = (TGraphAsymmErrors*)fhep_cms->Get("Table 33/Graph1D_y1");

 g->Draw("AP");
 ghep_cms->Draw("P");

 TCanvas* c = new TCanvas("c","c",600,800);
 c->Divide(3,3);
 for(int i=0;i<nPtTrig;i++){
	c->cd(i+1);
	hFPhiZYAM[0][i]->Draw();
 }

}
