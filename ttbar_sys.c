#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TAttFill.h"
#include "THStack.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "string"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TMath.h"
#include "TStyle.h"

#include "TRatioPlot.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TEntryList.h"
#include "TInterpreter.h"

#include "myplotfun.h"

void ttbar_sys(Double_t tagJ1, Double_t tagJ2, Double_t Mjj_cut_value, int region_index, int hp_in, int lp_in, int sel_in=1, string r_name="Syst/Norminal", string chain="s")
{
	Int_t nsteps = 40;
	Double_t basef = 500.0;
	Double_t widthM = 50.0;
//	Double_t base_deta = 1.0;
//	Double_t base_deta = 0.5;
//      Double_t width_deta = 0.05;
	Double_t base_deta = 0.0;
        Double_t width_deta = 0.1;
	Double_t sigma_M = 0.01;
	Double_t sigma_deta = 0.1;
	Double_t sb2_M = sigma_M*sigma_M;
	Double_t sb2_M2 = 0.05*0.05;
	Double_t sb2_M3 = 0.1*0.1;
        Double_t sb2_d = sigma_deta*sigma_deta;
	Double_t fatJ_Min = 60.0;
	Double_t fatJ_Max = 120.0;
	Double_t fatJ_pt_Min = 200.0;
	
	Double_t totalSize[6][40];
	Double_t event_d_eta[4][40];
	Double_t Bg_M[40];
	Double_t p1_M[40];
	Double_t p2_M[40];
	Double_t sig_M[40];
	Double_t p1_M2[40];
        Double_t p2_M2[40];
        Double_t sig_M2[40];
	Double_t p1_M3[40];
        Double_t p2_M3[40];
        Double_t sig_M3[40];
	Double_t cutM[40];

	Double_t Bg_d_eta[40];
	Double_t sig_d_eta[40];
	Double_t cut_d_eta[40];
	Double_t sig[40];

	Int_t counter_tt[6] = {0,0,0,0,0,0};
	
//	std::string PObject[3] = {"fatJ_","tagJ1_","tagJ2_"};
//	std::string Data[5] = {"pt","eta","phi","m","d_eta"};
//	TH1I myHist("h1", "Nominal", 10, -4.5, 5.5);
//	TH1I* myHist = new TH1I("h1", "ntuple", 10, -4.5, 5.5);

	TH1::SetDefaultSumw2();	

	Int_t name_p1 = TMath::Nint(Mjj_cut_value);
	Int_t name_p2 = TMath::Nint(tagJ2);

//	TH1D *HFrame[6] = {fh1,fh2,fh3,fh4,fh5,fh6};
	TH1D *HIST2[7][40] = {NULL};

	TChain *tc1 = new TChain("Nominal");
	TChain *tc2 = new TChain("Nominal");
	TChain *tc3 = new TChain("Nominal");
	TChain *tc4 = new TChain("Nominal");
	TChain *tc5 = new TChain("Nominal");
	TChain *tc6 = new TChain("Nominal");

//	tc1->Add("/afs/cern.ch/work/j/jzeng/public/VBS_test/data/data15-*.root");
//	tc1->Add("/afs/cern.ch/work/j/jzeng/public/VBS_test/data/data*.root");
//        tc1->Add("/afs/cern.ch/work/j/jzeng/public/VBS_test/data/MGPy8EG_WZjj_lvqq_EW6.root");
//	tc1->Add("/afs/cern.ch/work/j/jzeng/public/VBS_test/data/MGPy8EvtGen_WWjj_lvqq_EW6.root");
//        printf("%d;\n",tc1->GetNtrees());
//        printf("%lld;\n",tc1->GetEntries());
	
	if (chain == "f") {
//        tc1->Add("/eos/user/j/jzeng/BDT_SYS/ttbar-*.root");
	tc1->Add("/eos/user/j/jzeng/BDT_SYS/ttbar_s/PwHerwigpp-*.root");
//	tc2->Add("/eos/user/j/jzeng/BDT_SYS/ttbar_s/radLo-*.root");
//	tc2->Add("/eos/user/j/jzeng/BDT_SYS/ttbar_s/PwHerwigpp-*.root");
	tc2->Add("/eos/user/j/jzeng/BDT_SYS/ttbar_s/aMcAtNloHerwig-*.root");
//	tc2->Add("/eos/user/j/jzeng/BDT_SYS/ttbar_s/radHi-*.root");
//        tc3->Add("/afs/cern.ch/work/j/jzeng/public/VBS_test/data/W*_v221-*.root");
//        tc4->Add("/afs/cern.ch/work/j/jzeng/public/VBS_test/data/*_improved-*.root");
//        tc5->Add("/afs/cern.ch/work/j/jzeng/public/VBS_test/data/Z*_v221-*.root");
//        tc6->Add("/afs/cern.ch/work/j/jzeng/public/VBS_test/data/singletop*.root");
	}
	else {
//	tc1->Add("/eos/user/j/jzeng/BDT_SYS/ttbar-0.root");
//	tc1->Add("/eos/user/j/jzeng/BDT_SYS/ttbar_s/PwHerwigpp-27.root");
	tc1->Add("/eos/user/j/jzeng/BDT_SYS/ttbar_s/radLo-1.root");
//	tc2->Add("/eos/user/j/jzeng/BDT_SYS/ttbar_s/PwHerwigpp-27.root");
//	tc2->Add("/eos/user/j/jzeng/BDT_SYS/ttbar_s/aMcAtNloHerwig-27.root");
	tc2->Add("/eos/user/j/jzeng/BDT_SYS/ttbar_s/radHi-1.root");
//	tc3->Add("/afs/cern.ch/work/j/jzeng/public/VBS_test/data/We*_v221-7.root");
//	tc4->Add("/afs/cern.ch/work/j/jzeng/public/VBS_test/data/*_improved-2.root");
//	tc5->Add("/afs/cern.ch/work/j/jzeng/public/VBS_test/data/Z*_v221-4.root");
//	tc6->Add("/afs/cern.ch/work/j/jzeng/public/VBS_test/data/singletop_s-0.root");
	}

	printf("%d;\n",tc1->GetNtrees());
        printf("%lld;\n",tc1->GetEntries());
	printf("%d;\n",tc2->GetNtrees());
        printf("%lld;\n",tc2->GetEntries());
//	printf("%d;\n",tc3->GetNtrees());
//        printf("%lld;\n",tc3->GetEntries());
//        printf("%d;\n",tc4->GetNtrees());
//        printf("%lld;\n",tc4->GetEntries());
//        printf("%d;\n",tc5->GetNtrees());
//        printf("%lld;\n",tc5->GetEntries());
//        printf("%d;\n",tc6->GetNtrees());
//        printf("%lld;\n",tc6->GetEntries());
//        printf("%d;\n",tc2->GetNtrees());
//        printf("%lld;\n",tc2->GetEntries());

//	TChain *tc7 = new TChain("Nominal");
//        tc7->Add("/afs/cern.ch/work/j/jzeng/public/VBS_test/data/Py8EG_A14NNPDF23LO_jetjet_JZ*.root");
//	printf("%d;\n",tc7->GetNtrees());
//        printf("%lld;\n",tc7->GetEntries());
	

	if (fatJ_Max < fatJ_Min) {
        printf("Error: fatJ_Max < fatJ_Min! \n");
        return;
        }	

//	gStyle->SetOptFit(10);

	set_hist(HIST2[0]);
        set_hist(HIST2[1]);
        set_hist(HIST2[2]);
        set_hist(HIST2[3]);
        set_hist(HIST2[4]);
        set_hist(HIST2[5]);
	set_hist(HIST2[6]);

	m_fill_hist(tc1, HIST2[0],tagJ1,tagJ2,Mjj_cut_value,region_index,hp_in,lp_in,sel_in);
	m_fill_hist(tc2, HIST2[1],tagJ1,tagJ2,Mjj_cut_value,region_index,hp_in,lp_in,sel_in);
//	fill_hist(tc3, HIST2[2],tagJ1,tagJ2,Mjj_cut_value,region_index,hp_in,lp_in,sel_in);
//	fill_hist(tc4, HIST2[3],tagJ1,tagJ2,Mjj_cut_value,region_index,hp_in,lp_in,sel_in);
//	fill_hist(tc5, HIST2[4],tagJ1,tagJ2,Mjj_cut_value,region_index,hp_in,lp_in,sel_in);
//	fill_hist(tc6, HIST2[5],tagJ1,tagJ2,Mjj_cut_value,region_index,hp_in,lp_in,sel_in);
//	fill_hist(tc7, HIST2[6],tagJ1,tagJ2,Mjj_cut_value);


//	Float_t weight_r= 0.696-(0.000166)*x;
	TH1D* Wjets_mjj = (TH1D*)HIST2[2][5]->Clone("Wjets_mjj");

	for (Int_t k = 1; k <= HIST2[2][5]->GetNbinsX(); k++) {
	Double_t weight_r = 0.696-(0.000166)*(Double_t)HIST2[2][5]->GetBinCenter(k);
        Wjets_mjj->SetBinContent(k,(Double_t)HIST2[2][5]->GetBinContent(k)*weight_r);
        Wjets_mjj->SetBinError(k,(Double_t)HIST2[2][5]->GetBinError(k)*weight_r);

	}
	
	printf("Test: %d \n", HIST2[0][5]->GetNbinsX());	
	printf("Test: %d \n", HIST2[0][5]->GetXaxis()->FindBin(basef));
	printf("Test_c: %d;%d;%d \n", counter_tt[0], counter_tt[1], counter_tt[2]);

	for (int jj=0; jj<nsteps; jj++){
	event_d_eta[0][jj] = HIST2[0][4]->Integral(HIST2[0][4]->GetXaxis()->FindBin(base_deta + width_deta*jj), HIST2[0][4]->GetNbinsX());
        event_d_eta[1][jj] = HIST2[1][4]->Integral(HIST2[1][4]->GetXaxis()->FindBin(base_deta + width_deta*jj), HIST2[1][4]->GetNbinsX());
        event_d_eta[2][jj] = HIST2[2][4]->Integral(HIST2[2][4]->GetXaxis()->FindBin(base_deta + width_deta*jj), HIST2[2][4]->GetNbinsX());
	event_d_eta[3][jj] = HIST2[3][4]->Integral(HIST2[3][4]->GetXaxis()->FindBin(base_deta + width_deta*jj), HIST2[3][4]->GetNbinsX());
	Bg_d_eta[jj] = event_d_eta[1][jj] + event_d_eta[2][jj] + event_d_eta[3][jj];
	sig_d_eta[jj] = event_d_eta[0][jj]/sqrt(Bg_d_eta[jj]);
	cut_d_eta[jj] = base_deta + width_deta*jj;
	
	totalSize[0][jj] = HIST2[0][5]->Integral(HIST2[0][5]->GetXaxis()->FindBin(basef + widthM*jj), HIST2[0][5]->GetNbinsX());
	totalSize[1][jj] = HIST2[1][5]->Integral(HIST2[1][5]->GetXaxis()->FindBin(basef + widthM*jj), HIST2[1][5]->GetNbinsX());
//	totalSize[2][jj] = HIST2[2][5]->Integral(HIST2[2][5]->GetXaxis()->FindBin(basef + widthM*jj), HIST2[2][5]->GetNbinsX());
	totalSize[2][jj] = Wjets_mjj->Integral(Wjets_mjj->GetXaxis()->FindBin(basef + widthM*jj), Wjets_mjj->GetNbinsX());
	totalSize[3][jj] = HIST2[3][5]->Integral(HIST2[3][5]->GetXaxis()->FindBin(basef + widthM*jj), HIST2[3][5]->GetNbinsX());
	totalSize[4][jj] = HIST2[4][5]->Integral(HIST2[4][5]->GetXaxis()->FindBin(basef + widthM*jj), HIST2[4][5]->GetNbinsX());
	totalSize[5][jj] = HIST2[5][5]->Integral(HIST2[5][5]->GetXaxis()->FindBin(basef + widthM*jj), HIST2[5][5]->GetNbinsX());
	Bg_M[jj] = totalSize[1][jj] + totalSize[2][jj] + totalSize[3][jj] + totalSize[4][jj] + totalSize[5][jj];
//	p1_M[jj] = (s+b)*log((b+sb2_M)/((b*b/(s+b)+sb2_M)));
//	p2_M[jj] = log(1.0+sb2_M*s/(b*(b+sb2_M)))*b*b/sb2_M;
	p1_M[jj] = (totalSize[0][jj]+Bg_M[jj])*log(((1.0/Bg_M[jj])+sb2_M)/(1.0/(totalSize[0][jj]+Bg_M[jj])+sb2_M));
        p2_M[jj] = log(1.0+sb2_M*totalSize[0][jj]/(1.0+sb2_M*Bg_M[jj]))/sb2_M;
//	p1_M2[jj] = (totalSize[0][jj]+Bg_M[jj])*log(((1.0/Bg_M[jj])+sb2_M2)/(1.0/(totalSize[0][jj]+Bg_M[jj])+sb2_M2));
//        p2_M2[jj] = log(1.0+sb2_M2*totalSize[0][jj]/(1.0+sb2_M2*Bg_M[jj]))/sb2_M2;
	p1_M2[jj] = (totalSize[0][jj]+Bg_M[jj])*log(((totalSize[0][jj]+Bg_M[jj])*(Bg_M[jj]+sb2_M2*(Bg_M[jj]*Bg_M[jj])))/ (Bg_M[jj]*Bg_M[jj]+(totalSize[0][jj]+Bg_M[jj])*sb2_M2*(Bg_M[jj]*Bg_M[jj])));
	p2_M2[jj] = ((Bg_M[jj]*Bg_M[jj])/(sb2_M2*(Bg_M[jj]*Bg_M[jj])))* log (1+ ((sb2_M2*((Bg_M[jj]*Bg_M[jj]))*totalSize[0][jj])/(Bg_M[jj]*(Bg_M[jj]+sb2_M2*(Bg_M[jj]*Bg_M[jj])))));
	p1_M3[jj] = (totalSize[0][jj]+Bg_M[jj])*log(((1.0/Bg_M[jj])+sb2_M3)/(1.0/(totalSize[0][jj]+Bg_M[jj])+sb2_M3));
        p2_M3[jj] = log(1.0+sb2_M3*totalSize[0][jj]/(1.0+sb2_M3*Bg_M[jj]))/sb2_M3;
//	sig_M[jj] = totalSize[0][jj]/sqrt(Bg_M[jj]);
	sig_M[jj] = sqrt(2*((totalSize[0][jj]+Bg_M[jj])*log(1+totalSize[0][jj]/Bg_M[jj])-totalSize[0][jj]));
//	sig_M[jj] = sqrt(2*(p1_M[jj]-p2_M[jj]));
	sig_M2[jj] = sqrt(2*(p1_M2[jj]-p2_M2[jj]));
	sig_M3[jj] = sqrt(2*(p1_M3[jj]-p2_M3[jj]));
//	sig_M[jj] = sqrt(2*((totalSize[0][jj]+Bg_M[jj])*log((totalSize[0][jj]+Bg_M[jj])*(Bg_M[jj]+sb2_M)/(Bg_M[jj]*Bg_M[jj]+(totalSize[0][jj]+Bg_M[jj])*sb2_M))-((Bg_M[jj]*Bg_M[jj]/sb2_M)*log(1+(sb2_M*totalSize[0][jj])/(Bg_M[jj]*(Bg_M[jj]+sb2_M))))));
	cutM[jj] = basef + widthM*jj;	 	
	}

	tplot_and_save(HIST2[0][0],HIST2[1][0],2/sel_in,r_name,"fatJ_pt",  "fatJ_pt[Gev]");
	tplot_and_save(HIST2[0][1],HIST2[1][1],2/sel_in,r_name,"fatJ_eta", "fatJ_eta");
	tplot_and_save(HIST2[0][2],HIST2[1][2],2/sel_in,r_name,"fatJ_phi", "fatJ_phi");	
	tplot_and_save(HIST2[0][3],HIST2[1][3],2/sel_in,r_name,"fatJ_m",   "fatJ_m[GeV]");
	tplot_and_save(HIST2[0][4],HIST2[1][4],2/sel_in,r_name,"d_eta",    "d_eta");
	tplot_and_save(HIST2[0][6],HIST2[1][6],2/sel_in,r_name,"lvJ_m",    "lvJ_m[Gev]");
	tplot_and_save(HIST2[0][7],HIST2[1][7],1,       r_name,"Mjj_tag",  "Mjj_tag[Gev]");
	tplot_and_save(HIST2[0][8],HIST2[1][8],4/sel_in,r_name,"tagj1_pt", "tagj1_pt[Gev]");
	tplot_and_save(HIST2[0][9],HIST2[1][9],4/sel_in,r_name,"tagj2_pt", "tagj2_pt[Gev]");
	tplot_and_save(HIST2[0][10],HIST2[1][10],4/sel_in,r_name,"dphi_MTM",     "dphi_MTM");
	tplot_and_save(HIST2[0][11],HIST2[1][11],4/sel_in,r_name,"dphi_MET_Jets","dphi_MET_Jets");
	tplot_and_save(HIST2[0][12],HIST2[1][12],4/sel_in,r_name,"dphi_MET_Vhad","dphi_MET_Vhad");
	tplot_and_save(HIST2[0][13],HIST2[1][13],4/sel_in,r_name,"dphi_lep_MET", "dphi_lep_MET");
	tplot_and_save(HIST2[0][14],HIST2[1][14],4/sel_in,r_name,"W_MT",         "W_MT[Gev]");
	tplot_and_save(HIST2[0][15],HIST2[1][15],4/sel_in,r_name,"met",          "met[Gev]");
	tplot_and_save(HIST2[0][16],HIST2[1][16],2/sel_in,r_name,"sigj1_pt",     "sigj1_pt[Gev]");
	tplot_and_save(HIST2[0][17],HIST2[1][17],2/sel_in,r_name,"sigj2_pt",     "sigj2_pt[Gev]");
	tplot_and_save(HIST2[0][18],HIST2[1][18],2/sel_in,r_name,"lvjj_m",       "lvjj_m[Gev]");
	tplot_and_save(HIST2[0][19],HIST2[1][19],1,       r_name,"BDT_merged",   "BDT_merged");
	tplot_and_save(HIST2[0][20],HIST2[1][20],1,       r_name,"BDT_resolved", "BDT_resolved");
	tplot_and_save(HIST2[0][21],HIST2[1][21],2/sel_in,r_name,"lep_eta",      "lep_eta");
        tplot_and_save(HIST2[0][22],HIST2[1][22],1,r_name,"Njets",        "Njets");
        tplot_and_save(HIST2[0][23],HIST2[1][23],1,r_name,"Ntrkjets",     "Ntrkjets");
        tplot_and_save(HIST2[0][24],HIST2[1][24],2/sel_in,r_name,"xiV",          "xiV");
        tplot_and_save(HIST2[0][25],HIST2[1][25],2/sel_in,r_name,"tagj1_width",  "tagj1_width");
        tplot_and_save(HIST2[0][26],HIST2[1][26],2/sel_in,r_name,"tagj2_width",  "tagj2_width");
        tplot_and_save(HIST2[0][27],HIST2[1][27],2/sel_in,r_name,"sigj1_width",  "sigj1_width");
        tplot_and_save(HIST2[0][28],HIST2[1][28],2/sel_in,r_name,"sigj2_width",  "sigj2_width");
        tplot_and_save(HIST2[0][29],HIST2[1][29],2/sel_in,r_name,"res_deta_jj",  "res_deta_jj");
        tplot_and_save(HIST2[0][30],HIST2[1][30],2/sel_in,r_name,"res_dr_lv",    "res_dr_lv");
	tplot_and_save(HIST2[0][31],HIST2[1][31],1,r_name,       "lvjjjj_m",     "lvjjjj_m[Gev]");
	tplot_and_save(HIST2[0][32],HIST2[1][32],1,r_name,       "lvJjj_m",      "lvJjj_m[Gev]");

	printf("tagJ1_cut:%f; tagJ2_cut:%f; \n",tagJ1, tagJ2);
	printf("%s: %f;%f;%f;%f;%f;%f;%f;%f;%f;%f; \n", "M__",sig_M[20],sig_M[21],sig_M[22],sig_M[23],sig_M[24],sig_M[25],sig_M[26],sig_M[27],sig_M[28],sig_M[29]);
	printf("%s: %f;%f;%f;%f;%f;%f;%f;%f;%f;%f; \n", "5%_",sig_M2[20],sig_M2[21],sig_M2[22],sig_M2[23],sig_M2[24],sig_M2[25],sig_M2[26],sig_M2[27],sig_M2[28],sig_M2[29]);
	printf("%s: %f;%f;%f;%f;%f;%f;%f;%f;%f;%f; \n", "10%",sig_M3[20],sig_M3[21],sig_M3[22],sig_M3[23],sig_M3[24],sig_M3[25],sig_M3[26],sig_M3[27],sig_M3[28],sig_M3[29]);
	printf("Signal:    %f;%f;%f;%f;%f; \n", totalSize[0][0],totalSize[0][1],totalSize[0][2],totalSize[0][3],totalSize[0][4]);
        printf("ttbar:     %f;%f;%f;%f;%f; \n", totalSize[1][0],totalSize[1][1],totalSize[1][2],totalSize[1][3],totalSize[1][4]);
//        printf("W+jets:    %f;%f;%f;%f;%f; \n", totalSize[2][0],totalSize[2][1],totalSize[2][2],totalSize[2][3],totalSize[2][4]);
//	printf("diboson:   %f;%f;%f;%f;%f; \n", totalSize[3][0],totalSize[3][1],totalSize[3][2],totalSize[3][3],totalSize[3][4]);
//        printf("Z+jets:    %f;%f;%f;%f;%f; \n", totalSize[4][0],totalSize[4][1],totalSize[4][2],totalSize[4][3],totalSize[4][4]);
//	printf("singletop: %f;%f;%f;%f;%f; \n", totalSize[5][0],totalSize[5][1],totalSize[5][2],totalSize[5][3],totalSize[5][4]);



}
