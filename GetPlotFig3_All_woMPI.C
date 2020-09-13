void SetStyle(TH2D* h, int option=0){
  
 h->SetTitle("");
 h->GetXaxis()->SetTitle("#Delta#eta");
 h->GetYaxis()->SetTitle("#Delta#varphi");
 if( option == 0 ) h->GetZaxis()->SetTitle("(1/N_{#font[22]{trig}}) (d^{2}N_{#font[22]{pair}}/d#Delta#eta d#Delta#varphi)");
 else if( option == 1 ) h->GetZaxis()->SetTitle("S(#Delta#eta, #Delta#varphi)");
 else if( option == 2 ) h->GetZaxis()->SetTitle("B(#Delta#eta, #Delta#varphi) / B(0, 0)");
 h->GetXaxis()->CenterTitle();
 h->GetYaxis()->CenterTitle();
 
 h->GetXaxis()->SetTitleOffset(1.2);
 h->GetYaxis()->SetTitleOffset(1.2);
 h->GetZaxis()->SetTitleOffset(1.2);
 
 h->GetXaxis()->SetTitleSize( h->GetXaxis()->GetTitleSize()*1.5 );
 h->GetYaxis()->SetTitleSize( h->GetYaxis()->GetTitleSize()*1.5 );
 h->GetZaxis()->SetTitleSize( h->GetZaxis()->GetTitleSize()*1.5 );

// h->GetXaxis()->SetRangeUser(-2.0,2.0);
// h->GetYaxis()->SetRangeUser(-1.0,1.0);
 
 if( option==0 ) h->SetMaximum( h->GetBinContent( h->GetBin( h->GetXaxis()->FindBin(0.0), h->GetYaxis()->FindBin(1.0*3.141592) ) )*1.03 );

}

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
int MultMin[nMult] = {0};
int MultMax[nMult] = {1000};

const int nBin_for_ntrig = 50;

void GetPlotFig3_All_woMPI(){

 const int nfiles = 7;

 const int set_num[nfiles] = {
        5,
	5, 5, 5,
	5, 5, 5 };
 const int grp_num[nfiles] = {
	16,
	18, 17, 31,
	1,  3,  5 };
 TH1D* hNTrig[nfiles][nBin_for_ntrig];
 TH1D* hNTrigMerged[nfiles][nMult];

 TH2D* hCorr_Fig2[nfiles][nMult][nPtTrig];
 TH2D* hCorrMixing_Fig2[nfiles][nMult][nPtTrig];
 TH2D* hC[nfiles][nMult][nPtTrig];

 double NTrig[nfiles][nMult][nPtTrig];

 TH1D* hFPhi[nfiles][nMult][nPtTrig];
 TF1* fFourier[nfiles][nMult][nPtTrig];

 TH1D* hFPhiZYAM[nfiles][nMult][nPtTrig];
 double ZYAMPhi[nfiles][nMult][nPtTrig];
 double ZYAM[nfiles][nMult][nPtTrig];

 double YieldCntl[nfiles][nMult][nPtTrig];
 double YieldStat[nfiles][nMult][nPtTrig];

 TFile* fin;

 for(int f=0;f<nfiles;f++){
        if( f < 4 ) fin = new TFile(Form("./data/output_set%02d_grp%03d_ana4.root",set_num[f],grp_num[f]),"read");
	else if( f > 3 ) fin = new TFile(Form("./data/output_set%02d_grp%03d_ana5.root",set_num[f],grp_num[f]),"read");

        for(int i=0;i<nBin_for_ntrig;i++){
                hNTrig[f][i] = (TH1D*)fin->Get(Form("hNTrig_%d",i));
                hNTrig[f][i]->SetName(Form("hNTrig_%d_%d",f,i));
        }
        for(int i=0;i<nMult;i++){
                hNTrigMerged[f][i] = (TH1D*)hNTrig[0][0]->Clone();
                hNTrigMerged[f][i]->SetName(Form("hNTrigNorm_%d_%d",f,i));
                hNTrigMerged[f][i]->Reset();
                for(int j=0;j<nBin_for_ntrig;j++){
			if( j*5 + 2.5 > MultMin[i] && j*5 + 2.5 < MultMax[i] ){
//                        if( j*10 + 5 > MultMin[i] && j*10 + 5 < MultMax[i] ){
                                hNTrigMerged[f][i]->Add( hNTrig[f][j],1 );
                        }
                }
        }

        for(int i=0;i<nMult;i++){
                for(int j=0;j<nPtTrig;j++){
                        hCorr_Fig2[f][i][j] = (TH2D*)fin->Get(Form("hCorr_Fig3_%d",j));
                        hCorrMixing_Fig2[f][i][j] = (TH2D*)fin->Get(Form("hCorrMixing_Fig3_%d",j));
                        hC[f][i][j] = (TH2D*)hCorr_Fig2[f][i][j]->Clone();

                        hCorr_Fig2[f][i][j]->SetName(Form("hCorr_Fig3_%d_%d",f,j));
                        hCorrMixing_Fig2[f][i][j]->SetName(Form("hCorrMixing_Fig3_%d_%d",f,j));
                        hC[f][i][j]->SetName(Form("hC_%d_%d",f,j));


                        hCorrMixing_Fig2[f][i][j]->Scale( 1.0/hCorrMixing_Fig2[f][i][j]->GetMaximum() );
                        hC[f][i][j]->Divide( hCorrMixing_Fig2[f][i][j] );

                        hC[f][i][j]->RebinY(2);
                        hC[f][i][j]->Scale( 1, "width" );

                        NTrig[f][i][j] = 0.0;

                        for(int k=0;k<hNTrigMerged[f][i]->GetNbinsX();k++){
                                if( hNTrigMerged[f][i]->GetBinCenter(k+1) > pTMin[j]  && hNTrigMerged[f][i]->GetBinCenter(k+1) < pTMax[j] ){
                                        NTrig[f][i][j] += hNTrigMerged[f][i]->GetBinContent(k+1);
                                }
                        }

                        hC[f][i][j]->Scale( 1.0/NTrig[f][i][j] );
                        hC[f][i][j]->GetXaxis()->SetRangeUser(-4.0,4.0);

                        hFPhi[f][i][j] = (TH1D*)hC[f][i][j]->ProjectionY(Form("hFPhi_%d_%d_%d",f,i,j),hC[f][i][j]->GetXaxis()->FindBin(-4.0),hC[f][i][j]->GetXaxis()->FindBin(-2.0),"e" );
                        hFPhi[f][i][j]->Add( (TH1D*)hC[f][i][j]->ProjectionY(Form("hFPhi_%d_%d_%d",f,i,j),hC[f][i][j]->GetXaxis()->FindBin(2.0),hC[f][i][j]->GetXaxis()->FindBin(4.0),"e" ), 1.0 );
                        hFPhi[f][i][j]->Scale( hC[f][i][j]->GetXaxis()->GetBinWidth(1) );
                        hFPhi[f][i][j]->Scale( 1.0/2.0 );

                        fFourier[f][i][j] = new TF1("f1",Fourier,-10.0,10.0,nFourierTerm);
                        for(int n=0;n<nFourierTerm;n++) fFourier[f][i][j]->SetParLimits(n, -5, 5);

                        hFPhiZYAM[f][i][j] = (TH1D*)hFPhi[f][i][j]->Clone();
                        hFPhi[f][i][j]->Fit( fFourier[f][i][j], "" , "" );
                        ZYAMPhi[f][i][j] = fFourier[f][i][j]->GetMinimumX(-1.5,1.5);
                        ZYAM[f][i][j] = fFourier[f][i][j]->GetMinimum(-1.5,1.5);

                        for(int k=0;k<hFPhiZYAM[f][i][j]->GetNbinsX();k++){
                                hFPhiZYAM[f][i][j]->AddBinContent( k+1, -ZYAM[f][i][j] );
                        }

	                YieldCntl[f][i][j] = 0.0;
	                YieldStat[f][i][j] = 0.0;
	                for(int k=0;k<hFPhiZYAM[f][i][j]->GetNbinsX();k++){
	                        if( fabs( hFPhiZYAM[f][i][j]->GetBinCenter(k+1) ) > fabs( ZYAMPhi[f][i][j] ) ) continue;
	                        YieldCntl[f][i][j] += hFPhiZYAM[f][i][j]->GetBinContent(k+1) * hFPhiZYAM[f][i][j]->GetBinWidth(k+1);
	                        YieldStat[f][i][j] += pow( hFPhiZYAM[f][i][j]->GetBinError(k+1)*hFPhiZYAM[f][i][j]->GetBinWidth(k+1),2 );
	                }
	                YieldStat[f][i][j] = sqrt( YieldStat[f][i][j] );
                }
        }
 }


 double ptcntl[nPtTrig];
 double ptstat[nPtTrig];

 for(int i=0;i<nPtTrig;i++){
	ptcntl[i] = (pTMin[i]+pTMax[i])/2.0;
	ptstat[i] = (pTMax[i]-pTMin[i])/2.0;
 }

 TGraphErrors* g[nfiles];
 for(int f=0;f<nfiles;f++){
	g[f] = new TGraphErrors(nPtTrig,ptcntl,YieldCntl[f][0],ptstat,YieldStat[f][0]);
 }


 char ModnameFanc[nfiles][1000] = {
	"PYTHIA8 default, no MPI",
	"String Shoving, g=3, no MPI",
//	"String Shoving, g=4, no MPI",
	"String Shoving, g=5, no MPI",
//	"String Shoving, g=6, no MPI",
	"String Shoving, g=7, no MPI",

	"String Shoving, g=3, LM",
//        "String Shoving, g=4, LM",
        "String Shoving, g=5, LM",
//        "String Shoving, g=6, LM",
        "String Shoving, g=7, LM"
};

 TFile* fhep_cms = new TFile("data/CMSResults.root","read"); //33, 35
 TGraphAsymmErrors* ghep_cms = (TGraphAsymmErrors*)fhep_cms->Get("Table 33/Graph1D_y1");
 ghep_cms->GetXaxis()->SetRangeUser(0,6);
 ghep_cms->SetTitle("");

 ghep_cms->GetXaxis()->SetTitle("p_{T,trig(assoc)} (GeV/c)");
 ghep_cms->GetYaxis()->SetTitle("Ridge Yield");

 ghep_cms->GetXaxis()->SetTitleSize(0.06);
 ghep_cms->GetYaxis()->SetTitleSize(0.05);

 TLegend* leg = new TLegend(0.55,0.5,0.89,0.89);
 leg->SetLineWidth(0.0);
 leg->SetFillColorAlpha(0,0);
// leg->SetNColumns(2);

 TCanvas* c = new TCanvas("c","c",800,600);
 c->SetLeftMargin(0.14);
 c->SetBottomMargin(0.14);


 ghep_cms->SetMaximum(0.03);
 ghep_cms->SetMarkerStyle(20);
 ghep_cms->SetMarkerColor(1);
 ghep_cms->SetLineColor(1);
// ghep_cms->Draw("AP");

 g[0]->GetXaxis()->SetRangeUser(0,4.0);
 g[0]->SetTitle("");

 g[0]->GetXaxis()->SetTitle("p_{T,trig(assoc)} (GeV/c)");
 g[0]->GetYaxis()->SetTitle("Ridge Yield");

 g[0]->GetXaxis()->SetTitleSize(0.06);
 g[0]->GetYaxis()->SetTitleSize(0.05);
 g[0]->SetMaximum(0.007);
 g[0]->SetMinimum(-0.0002);
 g[0]->Draw("AL3");

// leg->AddEntry( ghep_cms, "CMS", "p");
 leg->SetHeader( "pp #sqrt{s}=13 TeV" ); 
 for(int f=0;f<nfiles;f++){
	g[f]->SetLineWidth(10.0);
	if( f == 0 ){
		g[f]->SetLineColor(f+2);
		g[f]->SetFillColor(f+2);
	}
	else{
		g[f]->SetLineColor(f%3+3);
		g[f]->SetFillColor(f%3+3);
	}
	if( f > 3 ){
		g[f]->SetLineStyle(10);
	}
	g[f]->Draw("L3");
	leg->AddEntry( g[f], ModnameFanc[f],"l");
 }
 leg->Draw();
 c->SaveAs("figs/fig3/All_pt_wompi.pdf");
 leg->Clear();


/*

 ghep_cms->Draw("AP");
 leg->AddEntry( ghep_cms, "CMS", "p");
 for(int f=1;f<6;f++){
        g[f]->SetMarkerColor(f+2);
//        g[f]->SetMarkerStyle(f+21);
        g[f]->SetLineColor(f+2);
        g[f]->Draw("P");
        leg->AddEntry( g[f], Modname[f],"p");
 }
 leg->Draw();
// c->SaveAs("figs/fig3/All_pt_loose_gamp.pdf");
 leg->Clear();

 leg->SetNColumns(1);
 ghep_cms->Draw("AP");
 leg->AddEntry( (TObject*)0, "#scale[1.3]{pp, #sqrt{s}=13 TeV}", "" );
 leg->AddEntry( ghep_cms, "CMS, N_{ch}>110", "p");
 for(int f=0;f<nfiles;f++){
	if( f!=3 && f!=6 && f!=7 ) continue;
        g[f]->SetMarkerColor(f+2);
        g[f]->SetLineColor(f+2);
	g[f]->SetFillColorAlpha(f+2,0.5);
        g[f]->Draw("3");
        leg->AddEntry( g[f], ModnameFanc[f],"f");
 }
 leg->Draw();
// c->SaveAs("figs/fig3/All_pt_loose_ptcut.pdf");
 leg->Clear();

// g[0]->SetFillStyle(1);
// g[3]->SetFillStyle(1);
 g[0]->SetLineColor(2);
 g[3]->SetLineColor(4);

 g[0]->SetFillColor(2);
 g[3]->SetFillColor(4);

 TLegend* leg2 = new TLegend(0.35,0.56,0.85,0.89);
 leg2->SetLineWidth(0.0);
 leg2->SetFillColorAlpha(0,0); 
 ghep_cms->Draw("AP");
 leg2->AddEntry( (TObject*)0, "pp, #sqrt{s}=13 TeV", "" );
 leg2->AddEntry( ghep_cms, "CMS, N_{ch}>110", "p" );
 leg2->AddEntry( g[0], "PYTHIA8 default, N_{ch}>90", "f" );
 leg2->AddEntry( g[3], "PYTHIA8 String Shoving, N_{ch}>90", "f" );
 g[0]->Draw("3C"); g[3]->Draw("3");
 leg2->Draw();

// c->SaveAs("figs/fig3/Reprodpt.pdf");
*/

 TCanvas* c1 = new TCanvas("c1","c1",800,600);
 gPad->SetLeftMargin(0.14);
 gStyle->SetOptStat(0);
 SetStyle(hC[0][0][2]); hC[0][0][2]->Draw("surf1");
 c1->SaveAs("corr_pd.pdf");
 SetStyle(hC[5][0][3]); hC[5][0][3]->Draw("surf1");
 c1->SaveAs("corr_nompi.pdf");

}
