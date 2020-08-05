{
 TFile* fin = new TFile("outfile_PYTHIA8_pp13TeV_template_fit.root","read");

 TGraphErrors* gdata_v22_pT = (TGraphErrors*)fin->Get("gdata_v22_pT");
 TGraphErrors* gdata_v22_mult = (TGraphErrors*)fin->Get("gdata_v22_mult");

 TGraphErrors* gv22_pT_PYTHIA = (TGraphErrors*)fin->Get("gv22_pT_PYTHIA");
 TGraphErrors* gv22_pT_PYTHIA_shoving = (TGraphErrors*)fin->Get("gv22_pT_PYTHIA_shoving");
 TGraphErrors* gv22_mult_PYTHIA = (TGraphErrors*)fin->Get("gv22_mult_PYTHIA");
 TGraphErrors* gv22_mult_PYTHIA_shoving = (TGraphErrors*)fin->Get("gv22_mult_PYTHIA_shoving");



 TCanvas* c = new TCanvas("c","c",800,600);
 c->SetLeftMargin(0.14);
 c->SetBottomMargin(0.14);

 gdata_v22_pT->GetYaxis()->SetTitle("v_{2,2}(p_{T}^{a},p_{T}^{b})");
 gdata_v22_pT->GetXaxis()->SetTitle("p_{T}^{a} (GeV/c)");
 gdata_v22_pT->GetXaxis()->SetTitleSize(0.06);
 gdata_v22_pT->GetYaxis()->SetTitleSize(0.05);

 TLegend* leg = new TLegend(0.15,0.65,0.6,0.89);
 leg->SetLineWidth(0.0);
 leg->SetFillColorAlpha(0,0);

 gdata_v22_pT->SetMaximum(0.01);
 gdata_v22_pT->SetMinimum(-0.0003);
 gdata_v22_pT->Draw("AP");
 gv22_pT_PYTHIA->Draw("L");
 gv22_pT_PYTHIA_shoving->Draw("L");


 leg->AddEntry( (TObject*)0, "pp, #sqrt{s}=13 TeV","");
 leg->AddEntry( (TObject*)0, "0.5<p_{T}^{b}<5.0, N_{ch}, ATLAS Template fit","");
 leg->AddEntry( gdata_v22_pT, "ATLAS", "p");
 leg->AddEntry( gv22_pT_PYTHIA, "PYTHIA8 default", "l");
 leg->AddEntry( gv22_pT_PYTHIA_shoving, "PYTHIA8 String Shoving","l");
 leg->Draw();


 c->SaveAs("atlaspt.pdf");

 gdata_v22_mult->GetYaxis()->SetTitle("v_{2,2}");
 gdata_v22_mult->GetXaxis()->SetTitle("N_{ch}");
 gdata_v22_mult->GetXaxis()->SetTitleSize(0.06);
 gdata_v22_mult->GetYaxis()->SetTitleSize(0.05);
 gdata_v22_mult->SetMaximum(0.01);
 gdata_v22_mult->SetMinimum(-0.0011);

 gdata_v22_mult->Draw("AP");
 gv22_mult_PYTHIA->Draw("L");
 gv22_mult_PYTHIA_shoving->Draw("L");

 TLegend* leg2 = new TLegend(0.15,0.65,0.6,0.89);
 leg2->SetLineWidth(0.0);
 leg2->SetFillColorAlpha(0,0);
 leg2->AddEntry( (TObject*)0, "pp, #sqrt{s}=13 TeV","");
 leg2->AddEntry( (TObject*)0, "0.5<p_{T}^{a,b}<5.0, ATLAS Template fit","");
 leg2->AddEntry( gdata_v22_mult, "ATLAS", "p");
 leg2->AddEntry( gv22_mult_PYTHIA, "PYTHIA8 default", "l");
 leg2->AddEntry( gv22_mult_PYTHIA_shoving, "PYTHIA8 String Shoving","l");
 leg2->Draw();

 c->SaveAs("atlasmult.pdf");

}
