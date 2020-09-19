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

 if( option==0 ) h->SetMaximum( h->GetBinContent( h->GetBin( h->GetXaxis()->FindBin(0.0), h->GetYaxis()->FindBin(1.0*3.141592) ) )*1.03 );

}

const int nFourierTerm = 4;
  
double Fourier(double *x, double *p){
 double f=0;
 f += p[0];
 for(int i=1;i<nFourierTerm;i++){ f += p[i]*cos(x[0]*(double)i); }
 return f;
}


const int nMult = 4;
const int nPtTrig = 4;

int MultMin[nMult] = {0 ,35 ,90 ,110};
int MultMax[nMult] = {35,90,110,1000};
float pTMin[nPtTrig] = { 0.1, 1.0, 2.0, 3.0 };
float pTMax[nPtTrig] = { 1.0, 2.0, 3.0, 4.0 };

const int nBin_for_ntrig = 50;

void GetPlotFig2_ForPaper(){

 const int nfiles = 8;

 const int set_num[nfiles] = {
	5, 5, 5, 5, 5,
	5, 5, 5 };
 const int grp_num[nfiles] = {
	0, 1, 2, 3, 4,
	5, 8, 13 };

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


 TFile* fin;

 for(int f=0;f<nfiles;f++){
// for(int f=0;f<1;f++){
	fin = new TFile(Form("./data/output_set%02d_grp%03d_ana0.root",set_num[f],grp_num[f]),"read");
	cout << Form("./data/output_set%02d_grp%03d_ana0.root",set_num[f],grp_num[f]) << endl;

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
//			if( j*10 + 5 > MultMin[i] && j*10 + 5 < MultMax[i] ){
				hNTrigMerged[f][i]->Add( hNTrig[f][j],1 );
			}
		}
	}

	for(int i=0;i<nMult;i++){
		for(int j=0;j<nPtTrig;j++){
			hCorr_Fig2[f][i][j] = (TH2D*)fin->Get(Form("hCorr_Fig2_%d_%d",i,j));
			hCorrMixing_Fig2[f][i][j] = (TH2D*)fin->Get(Form("hCorrMixing_Fig2_%d_%d",i,j));
			hC[f][i][j] = (TH2D*)hCorr_Fig2[f][i][j]->Clone();

			hCorr_Fig2[f][i][j]->SetName(Form("hCorr_Fig2_%d_%d_%d",f,i,j));
			hCorrMixing_Fig2[f][i][j]->SetName(Form("hCorrMixing_Fig2_%d_%d_%d",f,i,j));
			hC[f][i][j]->SetName(Form("hC_%d_%d_%d",f,i,j));


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
		}
	}
 } 


 TFile* fhep_cms = new TFile("data/CMSResults.root","read");
 TGraphAsymmErrors* ghep_cms[nMult][nPtTrig];
 for(int i=0;i<nMult;i++){
        for(int j=0;j<nPtTrig;j++){
		ghep_cms[i][j] = (TGraphAsymmErrors*)fhep_cms->Get(Form("Table %d/Graph1D_y1",i*8+j*2+1));
	}
 }

 TCanvas* c = new TCanvas("c","c",800,600);
 c->SetLeftMargin(0.14);
 c->SetBottomMargin(0.14);
 gStyle->SetOptStat(0);

 TLegend* leg = new TLegend(0.2,0.6,0.5,0.85);
 leg->SetLineWidth(0.0);
 leg->SetFillColorAlpha(0,0);
 leg->SetNColumns(2);

 double gAmp[nfiles] = {
	0, 3, 4, 5, 6, 
	7, 5,    5 };
 double pT[nfiles] = {
	0, 2, 2, 2, 2,
	2, 2.25, 2.5 };

 char Modname[nfiles][1000] = {
	"g=0, pTcut=0",
	"g=3, pTcut=2",
	"g=4, pTcut=2",
	"g=5, pTcut=2",
	"g=6, pTcut=2",

	"g=7, pTcut=2",
	"g=5, pTcut=2.25",
	"g=5, pTcut=2.5" };

 for(int i=0;i<nMult;i++){
	for(int j=0;j<nPtTrig;j++){
		leg->AddEntry( (TObject*)0, hCorr_Fig2[0][i][j]->GetTitle(), "" );
		ghep_cms[i][j]->SetMarkerStyle(20);
		ghep_cms[i][j]->SetMarkerColor(1);
		ghep_cms[i][j]->SetLineColor(1);
		ghep_cms[i][j]->Draw("AP");
		leg->AddEntry( ghep_cms[i][j], "CMS", "p");
		for(int f=0;f<nfiles;f++){
			if( !hFPhiZYAM[f][i][j] ) continue;
			hFPhiZYAM[f][i][j]->SetLineColor(f+2);
			hFPhiZYAM[f][i][j]->SetMarkerColor(f+2);
//			hFPhiZYAM[f][i][j]->SetMarkerStyle(21+f);
			hFPhiZYAM[f][i][j]->Draw("same");
			leg->AddEntry( hFPhiZYAM[f][i][j], Modname[f], "p");
		}
		leg->Draw();
//		c->SaveAs(Form("figs/fig2/fig2_dphi_all_%d_%d.pdf",i,j));
		leg->Clear();


                leg->AddEntry( (TObject*)0, hCorr_Fig2[0][i][j]->GetTitle(), "" );
                ghep_cms[i][j]->SetMarkerStyle(20);
                ghep_cms[i][j]->SetMarkerColor(1);
                ghep_cms[i][j]->SetLineColor(1);
                ghep_cms[i][j]->Draw("AP");
                leg->AddEntry( ghep_cms[i][j], "CMS", "p");
		for(int f=1;f<6;f++){
                        if( !hFPhiZYAM[f][i][j] ) continue;
                        hFPhiZYAM[f][i][j]->SetLineColor(f+2);
                        hFPhiZYAM[f][i][j]->SetMarkerColor(f+2);
//                        hFPhiZYAM[f][i][j]->SetMarkerStyle(21+f);
                        hFPhiZYAM[f][i][j]->Draw("same");
                        leg->AddEntry( hFPhiZYAM[f][i][j], Modname[f], "p");
		}
                leg->Draw();
//                c->SaveAs(Form("figs/fig2/fig2_dphi_gamp_%d_%d.pdf",i,j));
                leg->Clear();


                leg->AddEntry( (TObject*)0, hCorr_Fig2[0][i][j]->GetTitle(), "" );
                ghep_cms[i][j]->SetMarkerStyle(20);
                ghep_cms[i][j]->SetMarkerColor(1);
                ghep_cms[i][j]->SetLineColor(1);
                ghep_cms[i][j]->Draw("AP");
                leg->AddEntry( ghep_cms[i][j], "CMS", "p");
		for(int f=0;f<nfiles;f++){
			if( f!=3 && f!=6 && f!=7 ) continue;
                        if( !hFPhiZYAM[f][i][j] ) continue;
                        hFPhiZYAM[f][i][j]->SetLineColor(f+2);
                        hFPhiZYAM[f][i][j]->SetMarkerColor(f+2);
//                        hFPhiZYAM[f][i][j]->SetMarkerStyle(21);
                        hFPhiZYAM[f][i][j]->Draw("same");
                        leg->AddEntry( hFPhiZYAM[f][i][j], Modname[f], "p");
                }
                leg->Draw();
//                c->SaveAs(Form("figs/fig2/fig2_dphi_pt_%d_%d.pdf",i,j));
                leg->Clear();
		for(int f=0;f<nfiles;f++){
			SetStyle( hC[f][i][j] );
			hC[f][i][j]->Draw("surf1");
//			c->SaveAs(Form("figs/fig2/fig2_corr_%d_%d_%d.pdf",f,i,j));
		}
	}
 }

 SetStyle( hCorr_Fig2[3][1][1],1 );
 SetStyle( hCorrMixing_Fig2[3][1][1],2 );

 hCorr_Fig2[3][1][1]->Draw("surf1"); c->SaveAs(Form("figs/fig2/fig2_same_stringshoving.pdf"));
 hCorr_Fig2[0][1][1]->Draw("surf1"); c->SaveAs(Form("figs/fig2/fig2_same_pythiadefault.pdf"));

 hCorrMixing_Fig2[3][1][1]->Draw("surf1"); c->SaveAs(Form("figs/fig2/fig2_mixing_stringshoving.pdf"));
 hCorrMixing_Fig2[0][1][1]->Draw("surf1"); c->SaveAs(Form("figs/fig2/fig2_mixing_pythiadefault.pdf"));

 hFPhi[0][1][1]->Draw(); c->SaveAs(Form("figs/fig2/fig2_dphi_stringshoving.pdf"));
 hFPhiZYAM[0][1][1]->Draw(); c->SaveAs(Form("figs/fig2/fig2_dphzyami_stringshoving.pdf"));


 TCanvas* c1 = new TCanvas("c1","c1",800,600);
 gStyle->SetOptStat(0);
 gPad->SetLeftMargin(0.13);
 gPad->SetBottomMargin(0.13);
 c1->Divide(4,4,0,0);
 for(int i=0;i<nMult;i++){
        for(int j=0;j<nPtTrig;j++){
		c1->cd(1+j+4*i);
                ghep_cms[i][j]->SetTitle(Form("%s, %.1lf<p_{#font[22]{T}}<%.1lf, %d<N_{#font[22]{ch}}<%d",ghep_cms[i][j]->GetTitle(),pTMin[j],pTMax[j],MultMin[i],MultMax[i]));
                ghep_cms[i][j]->Draw("AP");
                hFPhiZYAM[0][i][j]->Draw("same");
//		hFPhiZYAM[3][i][j]->SetMarkerStyle(23);
		hFPhiZYAM[3][i][j]->SetMarkerColor(3);
		hFPhiZYAM[3][i][j]->SetLineColor(3);
		hFPhiZYAM[3][i][j]->Draw("same");
        }
 }
 c1->SaveAs("figs/fig2/allrange.pdf");

 TLegend* leg2 = new TLegend(0.176,0.475,0.978,0.813);
 leg2->SetLineWidth(0.0);
 leg2->SetFillColorAlpha(0,0);

 c->cd();
 for(int i=0;i<nMult;i++){
        for(int j=0;j<nPtTrig;j++){
                ghep_cms[i][j]->SetTitle(Form("%s, %.1lf<p_{#font[22]{T}}<%.1lf, %d<N_{#font[22]{ch}}<%d",ghep_cms[i][j]->GetTitle(),pTMin[j],pTMax[j],MultMin[i],MultMax[i]));
//		leg2->SetHeader( ghep_cms[i][j]->GetTitle() );
		leg2->AddEntry( (TObject*)0, "pp, #sqrt{s}=13 TeV", "" );
		leg2->AddEntry( (TObject*)0, Form("%.1lf<p_{#font[22]{T}}<%.1lf",pTMin[j],pTMax[j]), "" );
		leg2->AddEntry( (TObject*)0, Form("%d<N_{#font[22]{ch}}<%d",MultMin[i],MultMax[i]), "" );
		leg2->AddEntry( ghep_cms[i][j], "CMS", "p");
		ghep_cms[i][j]->SetTitle("");
		ghep_cms[i][j]->GetXaxis()->SetTitle("#Delta#varphi (rad)");
		ghep_cms[i][j]->GetYaxis()->SetTitle("(1/N_{trig}) (dN_{pair}/d#Delta#varphi) - C_{zyam}");
		ghep_cms[i][j]->GetXaxis()->SetTitleSize(0.06);
		ghep_cms[i][j]->GetYaxis()->SetTitleSize(0.05);
		ghep_cms[i][j]->SetMaximum( 0.085 );
		ghep_cms[i][j]->SetMinimum(-0.005 );
                ghep_cms[i][j]->GetXaxis()->SetNdivisions(505);
                ghep_cms[i][j]->GetYaxis()->SetNdivisions(505);
                ghep_cms[i][j]->GetXaxis()->CenterTitle();
                ghep_cms[i][j]->GetYaxis()->CenterTitle();
		ghep_cms[i][j]->GetXaxis()->SetLabelSize(0.05);
		ghep_cms[i][j]->GetYaxis()->SetLabelSize(0.05);
                ghep_cms[i][j]->Draw("AP");
		hFPhiZYAM[0][i][j]->SetFillColorAlpha(kRed,0.7);
                hFPhiZYAM[0][i][j]->Draw("e3,same");
//                hFPhiZYAM[3][i][j]->SetMarkerStyle(22);
                hFPhiZYAM[3][i][j]->SetMarkerColor(4);
                hFPhiZYAM[3][i][j]->SetLineColor(4);
		hFPhiZYAM[3][i][j]->SetFillColorAlpha(4,0.7);
                hFPhiZYAM[3][i][j]->Draw("e3,same");
		leg2->AddEntry( hFPhiZYAM[0][i][j], "PYTHIA8 default", "f");
		leg2->AddEntry( hFPhiZYAM[3][i][j], "PYTHIA8 String Shoving", "f");
		leg2->Draw();
		c->SaveAs(Form("figs/fig2/dphi_%d_%d.pdf",i,j));
		leg2->Clear();
        }
 }

 TCanvas* call = new TCanvas("call","call",1100,700);
 call->Divide(3,2,0,0);
 
 TH1D* hFPhiZYAM_ENLARGED_UNC[2][2][3];

 TLatex* latex = new TLatex();
 latex->SetTextFont(22);
 latex->SetTextSize(0.07);
 for(int i=0;i<2;i++){
	for(int j=0;j<3;j++){
		call->cd(j+1+i*3);
		if( j==0 ) gPad->SetLeftMargin(0.16);
		if( i==1 ) gPad->SetBottomMargin(0.15);
//		gPad->SetNdivision(505);
		ghep_cms[i*3][j+1]->Draw("AP");
		hFPhiZYAM[0][i*3][j+1]->Draw("e3,same");
		hFPhiZYAM[3][i*3][j+1]->Draw("e3,same");
		hFPhiZYAM_ENLARGED_UNC[0][i][j] = (TH1D*)hFPhiZYAM[0][i*3][j+1]->Clone();
		hFPhiZYAM_ENLARGED_UNC[1][i][j] = (TH1D*)hFPhiZYAM[3][i*3][j+1]->Clone();
		for(int p=0;p<hFPhiZYAM_ENLARGED_UNC[0][i][j]->GetNbinsX();p++){
			hFPhiZYAM_ENLARGED_UNC[0][i][j]->SetBinError( p+1, hFPhiZYAM_ENLARGED_UNC[0][i][j]->GetBinError(p+1)*10 );
			hFPhiZYAM_ENLARGED_UNC[1][i][j]->SetBinError( p+1, hFPhiZYAM_ENLARGED_UNC[1][i][j]->GetBinError(p+1)*10 );
		}
		if(i==0){
			hFPhiZYAM_ENLARGED_UNC[0][i][j]->Draw("e3,same");
			hFPhiZYAM_ENLARGED_UNC[1][i][j]->Draw("e3,same");
		}
		if( i==1 ) latex->SetTextSize(0.06);
		latex->DrawLatex(-1.3,0.075,Form("%d<#font[12]{N}_{ch}<%d, %.0f<#font[12]{p}_{T,trig(assoc)}<%.0f GeV/c",
			MultMin[i*3],MultMax[i*3],pTMin[j+1],pTMax[j+1]));
	}
 }
 TLegend* leg3 = new TLegend(0.176,0.475,0.978,0.813);
 leg3->SetLineWidth(0.0);
 leg3->SetFillColorAlpha(0,0);

 call->cd(1);
 leg3->AddEntry( (TObject*)0, "pp, #sqrt{s}=13 TeV", "" );
 leg3->AddEntry( ghep_cms[0][0], "CMS", "p");
 leg3->AddEntry( hFPhiZYAM_ENLARGED_UNC[0][0][0], "PYTHIA8 default", "l");
 leg3->AddEntry( hFPhiZYAM_ENLARGED_UNC[1][0][0], "PYTHIA8 String Shoving", "l");
 leg3->Draw(); 



 call->SaveAs("figs_paper/dphi_all.pdf");

}
