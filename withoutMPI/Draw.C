{
 TFile* fin = new TFile("AY-pT-CMS-kinematics.root","read");
 TH1D* hlow = (TH1D*)fin->Get("hYA_pT_w_MPI_LM");
 TH1D* hhigh = (TH1D*)fin->Get("hYA_pT_w_MPI_HM");
 TH1D* hNoMPI = (TH1D*)fin->Get("hYA_pT_wo_MPI");

 TCanvas* c = new TCanvas("c","c",800,600);
 c->SetLeftMargin(0.14);
 c->SetBottomMargin(0.14);

 hhigh->GetYaxis()->SetTitle("Ridge Yield");
 hhigh->GetXaxis()->SetTitle("p_{#font[22]{T,trig(assoc)}} (GeV/c)");
 hhigh->GetYaxis()->SetTitleSize(0.05);
 hhigh->GetXaxis()->SetTitleSize(0.06);


 TLegend* leg = new TLegend(0.42,0.49,0.86,0.78);
 leg->SetLineWidth(0.0);
 leg->SetFillColorAlpha(0,0);

 leg->SetHeader("Pythia8 pp, #sqrt{s}=13 TeV");
 leg->AddEntry( hhigh, "String Shoving, 90<N_{ch}<150","p");
 leg->AddEntry( hlow, "String Shoving, 0<N_{ch}<20","p");
 leg->AddEntry( hNoMPI, "String Shoving, no MPI","p");


 hhigh->Draw();
 hlow->Draw("same");
 hNoMPI->Draw("same");

 leg->Draw();

 c->SaveAs("mpi1.pdf");

}
