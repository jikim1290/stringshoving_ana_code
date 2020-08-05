#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>
#include <TClonesArray.h>
#include <TStopwatch.h>
#include <TAxis.h>
#include <THnSparse.h>
#include <THashList.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TRandom3.h>
#include <vector>
#include <deque>
#include <TSystem.h>

using namespace std;

TFile* fin;
TTree* tin;

int np;

float p_pt[10000];
float p_eta[10000];
float p_phi[10000];

TH1D* hMult = new TH1D("hMult","hMult",1000,0,1000);

struct Part{
	float pt;
	float eta;
	float phi;
};

typedef std::vector<Part> TRACKS;
typedef std::deque<TRACKS> EVENTPOOL;
typedef std::vector<EVENTPOOL> MIXINGPOOL;

int nMultBin;
int nPtTrigBin;
int nPtAsscoBin;

TRACKS *tracks;
MIXINGPOOL mixingpool;
EVENTPOOL *eventpool;

TAxis multaxis;
TAxis ptaxis;

std::vector<Part> GoodTracks;

double fMult;

int fbookingsize = 5;

double ptrange[9] = {
	0.1, 0.5, 1.0, 1.5, 2.0,
	2.5, 3.0, 4.0, 6.0 };

TH1D* hNTrig[50];


TH2D* hCorr_Fig2_CMS[4][4];
TH2D* hCorrMixing_Fig2_CMS[4][4];

TH2D* hCorr_Fig3_CMS[8];
TH2D* hCorrMixing_Fig3_CMS[8];

TH2D* hCorr_Fig4_CMS[20];
TH2D* hCorrMixing_Fig4_CMS[20];
//****************************************//
// functions
void ReadTree(TFile* fin){
 tin = (TTree*)fin->Get("T");

 tin->SetBranchAddress("np",&np);
 tin->SetBranchAddress("p_pt",p_pt);
 tin->SetBranchAddress("p_eta",p_eta);
 tin->SetBranchAddress("p_phi",p_phi);

}

int GetMult(int IsCMS){
 int mult = 0;
 if( IsCMS ){
	for(int i=0;i<np;i++){
		if( fabs( p_eta[i] ) < 2.4 && p_pt[i] > 0.4 ) mult++;
	}
 }
 return mult;
}
void GenerateHists_CMS(int options){
 if( options == 0 ){
	const int nMult = 4;
	const int nPtTrig = 4;
	const int nPtAssco = 4;

	nMultBin = nMult;
	nPtTrigBin = nPtTrig;
	nPtAsscoBin = nPtAssco;

	int MultMin[nMult] = {0 ,35 ,90 ,110};
	int MultMax[nMult] = {35,90,110,1000};
	float pTMin[nPtTrig] = { 0.1, 1.0, 2.0, 3.0 };
	float pTMax[nPtTrig] = { 1.0, 2.0, 3.0, 4.0 };

	double multbin[5] = {0,35,90,110,1000};
	double ptbin[5] = {0.1, 1.0, 2.0, 3.0, 4.0};

	multaxis = TAxis(nMult,multbin);
	ptaxis = TAxis(nPtTrig,ptbin);

	mixingpool.resize( nMultBin );

	float pTRange[nPtTrig+1];
	for(int i=0;i<nPtTrig;i++){ pTRange[i] = pTMin[i]; }
	pTRange[nPtTrig] = pTMax[nPtTrig-1];

	for(int i=0;i<nMult;i++){
		for(int j=0;j<nPtTrig;j++){
			hCorr_Fig2_CMS[i][j] = new TH2D(Form("hCorr_Fig2_%d_%d",i,j),Form("hCorr_Fig2, %d<N_{ch}<%d, %.1f<p_{T}<%.1f",MultMin[i],MultMax[i],pTMin[j],pTMax[j]),200,-5.0,5.0,50,-TMath::Pi()*0.5,TMath::Pi()*1.5);
			hCorrMixing_Fig2_CMS[i][j] = new TH2D(Form("hCorrMixing_Fig2_%d_%d",i,j),Form("hCorrMixing_Fig2, %d<N_{ch}<%d, %.1f<p_{T}<%.1f",MultMin[i],MultMax[i],pTMin[j],pTMax[j]),200,-5.0,5.0,50,-TMath::Pi()*0.5,TMath::Pi()*1.5);
		}
	}
 } else if( options == 1 ){
	const int nMult = 1;
	const int nPtTrig = 8;
	const int nPtAssco = 8;

        nMultBin = nMult;
        nPtTrigBin = nPtTrig;
        nPtAsscoBin = nPtAssco;

	float pTMin[nPtTrig] = {0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0 };
	float pTMax[nPtTrig] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0 };

	double multbin[2] = {110,1000};
	double ptbin[9] = {0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0};

	multaxis = TAxis(1,multbin);
	ptaxis = TAxis(nPtTrig,ptbin);

	mixingpool.resize( nMultBin );

	float pTRange[nPtTrig+1];
	for(int i=0;i<nPtTrig;i++){ pTRange[i] = pTMin[i]; }
	pTRange[nPtTrig] = pTMax[nPtTrig-1];

	for(int i=0;i<nPtTrig;i++){
		hCorr_Fig3_CMS[i] = new TH2D(Form("hCorr_Fig3_%d",i),Form("hCorr_Fig3, N_{ch}>105, %.lf<p_{T}<%.1f",pTMin[i],pTMax[i]),200,-5.0,5.0,50,-TMath::Pi()*0.5,TMath::Pi()*1.5);
		hCorrMixing_Fig3_CMS[i] = new TH2D(Form("hCorrMixing_Fig3_%d",i),Form("hCorrMixing_Fig3, N_{ch}>105, %.lf<p_{T}<%.1f",pTMin[i],pTMax[i]),200,-5.0,5.0,50,-TMath::Pi()*0.5,TMath::Pi()*1.5);
	}
 } else if( options == 2 ){
	const int nMult = 20;
	const int nPtTrig = 1;
	const int nPtAssco = 1;

        nMultBin = nMult;
        nPtTrigBin = nPtTrig;
        nPtAsscoBin = nPtAssco;

	multaxis = TAxis(20, 0, 200);
	ptaxis = TAxis(1, 1.0, 2.0);

	mixingpool.resize( nMultBin );

	for(int i=0;i<nMult;i++){
		hCorr_Fig4_CMS[i] = new TH2D(Form("hCorr_Fig4_%d",i),Form("hCorr_Fig4, %d<N_{ch}<%d, 1<p_{T}<2",i*10,i*10+10),200,-5.0,5.0,50,-TMath::Pi()*0.5,TMath::Pi()*1.5);
		hCorrMixing_Fig4_CMS[i] = new TH2D(Form("hCorrMixing_Fig4_%d",i),Form("hCorrMixing_Fig4, %d<N_{ch}<%d, 1<p_{T}<2",i*10,i*10+10),200,-5.0,5.0,50,-TMath::Pi()*0.5,TMath::Pi()*1.5);
	}
 } else if( options > 2.5 ){
        const int nMult = 1;
        const int nPtTrig = 8;
        const int nPtAssco = 8;

        nMultBin = nMult;
        nPtTrigBin = nPtTrig;
        nPtAsscoBin = nPtAssco;

        float pTMin[nPtTrig] = {0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0 };
        float pTMax[nPtTrig] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0 };

        double multbin[2];
	if( options == 3 ){
		multbin[0] = 90;
		multbin[1] = 1000;
	}
	else if( options == 4 ){
		multbin[0] = 0;
		multbin[1] = 1000;
	}
	else if( options == 5 ){
		multbin[0] = 0;
		multbin[1] = 20;
	}

        double ptbin[9] = {0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0};

        multaxis = TAxis(1,multbin);
        ptaxis = TAxis(nPtTrig,ptbin);

        mixingpool.resize( nMultBin );

        float pTRange[nPtTrig+1];
        for(int i=0;i<nPtTrig;i++){ pTRange[i] = pTMin[i]; }
        pTRange[nPtTrig] = pTMax[nPtTrig-1];

        for(int i=0;i<nPtTrig;i++){
                hCorr_Fig3_CMS[i] = new TH2D(Form("hCorr_Fig3_%d",i),Form("hCorr_Fig3, N_{ch}>105, %.lf<p_{T}<%.1f",pTMin[i],pTMax[i]),200,-5.0,5.0,50,-TMath::Pi()*0.5,TMath::Pi()*1.5);
                hCorrMixing_Fig3_CMS[i] = new TH2D(Form("hCorrMixing_Fig3_%d",i),Form("hCorrMixing_Fig3, N_{ch}>105, %.lf<p_{T}<%.1f",pTMin[i],pTMax[i]),200,-5.0,5.0,50,-TMath::Pi()*0.5,TMath::Pi()*1.5);
        }
 }


}

bool ObtainGoodTracks(TTree* tin){
 TRACKS *etl;
 EVENTPOOL *ep;

 if( multaxis.FindBin( fMult ) >= 1 && multaxis.FindBin( fMult ) <= multaxis.GetNbins() ){
	ep = &mixingpool[multaxis.FindBin( fMult )-1];
	ep -> push_back( TRACKS() ); //
	etl = &(ep->back());
 }

 GoodTracks.clear();
 Part part;
 for(int i=0;i<np;i++){
	if( fabs( p_eta[i] ) > 2.4 ) continue;
	if( p_pt[i] < 0.1 ) continue;

	part.pt = p_pt[i];
	part.eta = p_eta[i];
	part.phi = p_phi[i];

	GoodTracks.push_back( part );
	etl->push_back( part );

	hNTrig[ (int)fMult / 5 ]->Fill( p_pt[i] );
 }
 if( !GoodTracks.size() ) ep->pop_back();
 if( ep->size() > fbookingsize ){
	ep->pop_front();
 }

 return GoodTracks.size();
}

void GetCorrelations(int options){

 Part pa;
 Part pb;

 for(int i=0;i<GoodTracks.size()-1;i++){
	pa = (Part)GoodTracks.at(i);
//	if( !pa ) continue;
	for(int j=i+1;j<GoodTracks.size();j++){
		pb = (Part)GoodTracks.at(j);
//		if( !pb ) continue;

		double DeltaEta = pa.eta - pb.eta;
		double DeltaPhi = pa.phi - pb.phi;

		if( pa.pt < pb.pt ){
			DeltaEta *= -1.0;
			DeltaPhi *= -1.0;
		}
		DeltaPhi = TVector2::Phi_0_2pi(DeltaPhi);
		if( DeltaPhi > 1.5*TMath::Pi() ) DeltaPhi -= 2.0*TMath::Pi();

		if( options==0 ){
			if( ptaxis.FindBin( pa.pt ) != ptaxis.FindBin( pb.pt ) ) continue;
			if( ptaxis.FindBin( pa.pt ) < 1 || ptaxis.FindBin( pa.pt ) > ptaxis.GetNbins() ) continue;
			hCorr_Fig2_CMS[multaxis.FindBin( fMult )-1][ptaxis.FindBin( pa.pt )-1]->Fill( DeltaEta, DeltaPhi );
		} else if( options==1 || options==3 || options==4 || options==5 ){
			if( ptaxis.FindBin( pa.pt ) != ptaxis.FindBin( pb.pt ) ) continue;
			if( ptaxis.FindBin( pa.pt ) < 1 || ptaxis.FindBin( pa.pt ) > ptaxis.GetNbins() ) continue;
			if( multaxis.FindBin( fMult ) != 1 ) continue;
			hCorr_Fig3_CMS[ptaxis.FindBin( pa.pt )-1]->Fill( DeltaEta, DeltaPhi );
		} else if( options==2 ){
			if( ptaxis.FindBin( pa.pt ) != ptaxis.FindBin( pb.pt ) ) continue;
			if( ptaxis.FindBin( pa.pt ) < 1 || ptaxis.FindBin( pa.pt ) > ptaxis.GetNbins() ) continue;
			if( ptaxis.FindBin( pa.pt ) != 1 ) continue;
			hCorr_Fig4_CMS[multaxis.FindBin( fMult )-1]->Fill( DeltaEta, DeltaPhi );
		}
	}
 }

 TRACKS trackpool;

 int epsize = 1;
 if( multaxis.FindBin( fMult ) >= 1 && multaxis.FindBin( fMult ) <= multaxis.GetNbins() ){
	EVENTPOOL &ep = mixingpool[multaxis.FindBin( fMult )-1];
	epsize = ep.size();

	if (ep.size() < fbookingsize  ) return;
	int n = 0;
	for (auto pool: ep){
		if (n == (ep.size() -1 )) continue;
		for (auto part: pool) trackpool.push_back( part );
		n++;
	}
 }

 Part pMixing;
 for(int i=0;i<GoodTracks.size();i++){
	pa = (Part)GoodTracks.at(i);
	for(int j=0;j<trackpool.size();j++){
		pMixing = (Part)trackpool.at(j);

		double DeltaEta = pa.eta - pMixing.eta;
		double DeltaPhi = pa.phi - pMixing.phi;

		if( pa.pt < pMixing.pt ){
			DeltaEta *= -1.0;
			DeltaPhi *= -1.0;
		}
		DeltaPhi = TVector2::Phi_0_2pi(DeltaPhi);
		if( DeltaPhi > 1.5*TMath::Pi() ) DeltaPhi -= 2.0*TMath::Pi();

		if( options==0 ){
			if( ptaxis.FindBin( pa.pt ) != ptaxis.FindBin( pMixing.pt ) ) continue;
			if( ptaxis.FindBin( pa.pt ) < 1 || ptaxis.FindBin( pa.pt ) > ptaxis.GetNbins() ) continue;
			hCorrMixing_Fig2_CMS[multaxis.FindBin( fMult )-1][ptaxis.FindBin( pa.pt )-1]->Fill( DeltaEta, DeltaPhi );
		} else if( options==1 || options==3 || options==4 || options==5 ){
			if( ptaxis.FindBin( pa.pt ) != ptaxis.FindBin( pMixing.pt ) ) continue;
			if( ptaxis.FindBin( pa.pt ) < 1 || ptaxis.FindBin( pa.pt ) > ptaxis.GetNbins() ) continue;
                        if( multaxis.FindBin( fMult ) != 1 ) continue;
                        hCorrMixing_Fig3_CMS[ptaxis.FindBin( pa.pt )-1]->Fill( DeltaEta, DeltaPhi );
                } else if( options==2 ){
			if( ptaxis.FindBin( pa.pt ) != ptaxis.FindBin( pMixing.pt ) ) continue;
                        if( ptaxis.FindBin( pa.pt ) < 1 || ptaxis.FindBin( pa.pt ) > ptaxis.GetNbins() ) continue;
                        if( ptaxis.FindBin( pa.pt ) != 1 ) continue;
                        hCorrMixing_Fig4_CMS[multaxis.FindBin( fMult )-1]->Fill( DeltaEta, DeltaPhi );
                }
	}
 }
}

//****************************************//
// main
//void ana(){
int main(int argc, char** argv){

 int nruns = atoi(argv[1]);

 int irun = atoi(argv[2]) * nruns;
 int frun = irun + nruns;
 int anaind = atoi(argv[3]);
 int opt1 = atoi(argv[4]);
 int opt2 = atoi(argv[5]);

 TString fname;
 for(int i=0;i<50;i++) hNTrig[i] = new TH1D(Form("hNTrig_%d",i),Form("hNTrig_%d, %d<N_{ch}<%d",i,i*5,i*5+5),100,0,10);

 GenerateHists_CMS(anaind);
 for(int f=irun;f<frun;f++){
	fname = Form("/alice/data/shlim/pp13TeV_set%02d_grp%03d/outfile_pp13TeV_set%02d_grp%03d_%05d.root",opt1,opt2,opt1,opt2,f);
	if( !gSystem->IsFileInIncludePath( fname.Data() ) ) continue;
	fin = new TFile( fname.Data(), "read");
	ReadTree(fin);
	for(int i=0;i<tin->GetEntries();i++){
		tin->GetEntry(i);
		fMult = (double)GetMult(1);
		hMult->Fill( fMult );
		if( fMult > 100 ) cout << fMult << ", " << multaxis.FindBin( fMult ) << endl;
		if( !(multaxis.FindBin( fMult ) >= 1 && multaxis.FindBin( fMult ) <= multaxis.GetNbins()) ) continue;
		if( ObtainGoodTracks(tin) ){
			GetCorrelations(anaind);
		}
	}
 }
 TFile* fout = new TFile(Form("output_set%02d_grp%03d_%d_%d_ana%d.root",opt1,opt2,irun,frun,anaind),"recreate");
 hMult->Write();
 for(int i=0;i<50;i++) hNTrig[i]->Write();
 for(int i=0;i<4;i++){
	for(int j=0;j<4;j++){
		if( hCorr_Fig2_CMS[i][j]) hCorr_Fig2_CMS[i][j]->Write();
		if( hCorrMixing_Fig2_CMS[i][j] )hCorrMixing_Fig2_CMS[i][j]->Write();
	}
 }
 for(int i=0;i<8;i++){
	if( hCorr_Fig3_CMS[i] ) hCorr_Fig3_CMS[i]->Write();
	if( hCorrMixing_Fig3_CMS[i] )hCorrMixing_Fig3_CMS[i]->Write();
 }
 for(int i=0;i<20;i++){
	if( hCorr_Fig4_CMS[i] ) hCorr_Fig4_CMS[i]->Write();
	if( hCorrMixing_Fig4_CMS[i] )hCorrMixing_Fig4_CMS[i]->Write();
 }

 fout->Close();
}
