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

#include "TRatioPlot.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TEntryList.h"
#include "TInterpreter.h"

//class myplotfun {
//
//public:
//
//void fill_hist(TChain *tch, TH1D *HIST[30], Double_t J1_cut, Double_t J2_cut, Double_t Mjj_cut,int region, int hp, int lp);
//void res_fill_hist(TChain *tch, TH1D *HIST[30], Double_t J1_cut, Double_t J2_cut, Double_t Mjj_cut);
//void top_fill_hist(TChain *tch, TH1D *HIST[30], Double_t J1_cut, Double_t J2_cut, Double_t Mjj_cut);
//void normal(TH1D* hh0,TH1D* hh1,TH1D* hh2,TH1D* hh3,TH1D* hh4,TH1D* hh5,TLegend* leg, float scale, int rebin);
//void stack_plot(TH1D* h0,TH1D* h1,TH1D* h2,TH1D* h3,TH1D* h4,TH1D* h5,TLegend* leg, float scale, int rebin);
//void s_normal(TH1D* hh0,TH1D* hh1, TLegend* leg, float scale, int rebin);
//void AddRatioPlot(TH1D* hh0,TH1D* hh1,TH1D* h2,TH1D* hh3,TH1D* hh4,TH1D* hh5, int rebin, TString y_name);
//void s_AddRatioPlot(TH1D* hh0,TH1D* hh1, int rebin);
//void Multi_j(TH1D* hh0,int rebin);
//void CreateSubPad(TCanvas *canvas, Double_t vfrac=0.25);
//void Create2Pad(TCanvas *canvas, Double_t vfrac=0.25);
//void AddLatex();
//void set_hist(TH1D *HIST[30]);
//void plot_and_save(TH1D* hh0,TH1D* hh1,TH1D* hh2,TH1D* hh3,TH1D* hh4,TH1D* hh5,TH1D* hh6, int rebin, TString name, TString unit);
//void splot_and_save(TH1D* hh0,TH1D* hh1, int rebin, TString name, TString unit);
//void top_and_save(TH1D* hh0,TH1D* hh1,TH1D* hh2,TH1D* hh3,TH1D* hh4,TH1D* hh5,TH1D* hh6, int rebin, TString name, TString unit);
//
//};

void m_fill_hist(TChain *tch, TH1D *HIST[40], Double_t J1_cut, Double_t J2_cut, Double_t Mjj_cut,int region, int hp, int lp, int selection)
{

	TTreeReader myReader(tch);

        TTreeReaderValue<float> fatJ_pt(myReader,  "fatJ_pt");
        TTreeReaderValue<float> fatJ_eta(myReader, "fatJ_eta");
        TTreeReaderValue<float> fatJ_phi(myReader, "fatJ_phi");
        TTreeReaderValue<float> fatJ_m(myReader,   "fatJ_m");
        TTreeReaderValue<float> met(myReader,      "met");
//        TTreeReaderValue<float> met_phi(myReader,  "met_phi");
        TTreeReaderValue<float> lep_pt(myReader,   "lep_pt");
        TTreeReaderValue<float> lep_eta(myReader,  "lep_eta");
//        TTreeReaderValue<float> lep_phi(myReader,  "lep_phi");
        TTreeReaderValue<float> weight(myReader,   "weight");

        TTreeReaderValue<float> lvJjjmass(myReader,  "lvJjjmass"); 
	TTreeReaderValue<float> lvJmass(myReader,    "lvJmass");
	TTreeReaderValue<float> lv4jmass(myReader,   "lvjjjjmass");
	TTreeReaderValue<float> lvjjmass(myReader,   "lvjjmass");
        TTreeReaderValue<int> lepiso(myReader,   "lepiso");
        TTreeReaderValue<int> pass_isWZJet(myReader,   "pass_isWZJet");
//	TTreeReaderValue<int> pass_isWZJetLP(myReader,   "pass_isWZJetLP");
	TTreeReaderValue<int> pass_isWZSub_SB(myReader,"pass_isWZSub_MassSideBand");
//	TTreeReaderValue<int> pass_isWZSub_SBLP(myReader,"pass_isWZSub_MassSideBandLP");
        TTreeReaderValue<int> NBjets(myReader,   "NBjets");
	TTreeReaderValue<int> Njets(myReader,   "Njets");
	TTreeReaderValue<int> Ntrkjets(myReader,   "Ntrkjets");
//	TTreeReaderValue<float> dphi_MTM(myReader, "dphi_MET_TrackMET");
//        TTreeReaderValue<float> dphi_MET_Jets(myReader, "dphi_MET_Jets");
//        TTreeReaderValue<float> dphi_MET_Vhad(myReader, "dphi_MET_Vhad");
//        TTreeReaderValue<float> dphi_lep_MET(myReader, "dphi_lep_MET");
//        TTreeReaderValue<float> W_MT(myReader,   "W_MT");

	TTreeReaderValue<float> boosted_mjj_tag(myReader,  "boosted_mjj_tag");

        TTreeReaderValue<float> tagJ1_pt(myReader,  "boosted_tagJ1_pt");
        TTreeReaderValue<float> tagJ1_eta(myReader, "boosted_tagJ1_eta");
        TTreeReaderValue<float> tagJ1_phi(myReader, "boosted_tagJ1_phi");
        TTreeReaderValue<float> tagJ1_m(myReader,   "boosted_tagJ1_m");

        TTreeReaderValue<float> tagJ2_pt(myReader,  "boosted_tagJ2_pt");
        TTreeReaderValue<float> tagJ2_eta(myReader, "boosted_tagJ2_eta");
        TTreeReaderValue<float> tagJ2_phi(myReader, "boosted_tagJ2_phi");
        TTreeReaderValue<float> tagJ2_m(myReader,   "boosted_tagJ2_m");

	TTreeReaderValue<float> sigJ1_pt(myReader,  "sigJ1_pt");
        TTreeReaderValue<float> sigJ2_pt(myReader,  "sigJ2_pt");
        TTreeReaderValue<float> sigJJ_m(myReader,   "sigJJ_m");
	TTreeReaderValue<float> sigJ1_width(myReader,  "sigJ1_width");
	TTreeReaderValue<float> sigJ2_width(myReader,  "sigJ2_width");

	TTreeReaderValue<float> res_tagJ1_pt(myReader,  "resolved_tagJ1_pt");
        TTreeReaderValue<float> res_tagJ1_eta(myReader, "resolved_tagJ1_eta");
//        TTreeReaderValue<float> res_tagJ1_phi(myReader, "resolved_tagJ1_phi");
        TTreeReaderValue<float> res_tagJ1_m(myReader,   "resolved_tagJ1_m");
	TTreeReaderValue<float> res_tagJ1_width(myReader,   "resolved_tagJ1_width");
        TTreeReaderValue<float> res_tagJ2_pt(myReader,  "resolved_tagJ2_pt");
        TTreeReaderValue<float> res_tagJ2_eta(myReader, "resolved_tagJ2_eta");
//        TTreeReaderValue<float> res_tagJ2_phi(myReader, "resolved_tagJ2_phi");
        TTreeReaderValue<float> res_tagJ2_m(myReader,   "resolved_tagJ2_m");
	TTreeReaderValue<float> res_tagJ2_width(myReader,   "resolved_tagJ2_width");

	TTreeReaderValue<float> res_mjj_tag(myReader,       "resolved_mjj_tag");
//	TTreeReaderValue<float> res_dphi_MET_Jets(myReader, "resolved_dphi_MET_Jets");
//        TTreeReaderValue<float> res_dphi_MET_Vhad(myReader, "resolved_dphi_MET_Vhad");
//	TTreeReaderValue<float> res_deta_jj(myReader, "resolved_deta_jj");
	TTreeReaderValue<float> res_dr_lv(myReader, "resolved_dr_lv");

	TTreeReaderValue<int> ONELARGEJET(myReader,   "ONELARGEJET");
	TTreeReaderValue<int> TWOSIGNALJETS(myReader,   "TWOSIGNALJETS");

	TTreeReaderValue<float> BDT_resolved(myReader, "BDT_resolved");
        TTreeReaderValue<float> BDT_merged(myReader, "BDT_merged");
	TTreeReaderValue<float> boosted_xiV(myReader, "boosted_xiV");
	TTreeReaderValue<float> res_xiV(myReader, "resolved_xiV");

        TLorentzVector v1;
        TLorentzVector v2;

	while (myReader.Next()) {

                v1.SetPtEtaPhiM(*tagJ1_pt,*tagJ1_eta,*tagJ1_phi,*tagJ1_m);
                v2.SetPtEtaPhiM(*tagJ2_pt,*tagJ2_eta,*tagJ2_phi,*tagJ2_m);


                if ( *met > 80.0 ) {
                if ( *lep_pt > 27.0 && *lepiso == 1 ) {
		if ( (*NBjets==0 && region!=1)||(*NBjets>0 && region==1) ) {
		
		if ( selection==1){
			if ( *ONELARGEJET==1) {
			if ( *fatJ_pt > 200 && *fatJ_eta < 2.0 ) {
                	if ( (((*pass_isWZJet==1))&&(region!=2))||(((*pass_isWZSub_SB==1))&&(region==2)) ) {
			if ( *tagJ1_pt > J1_cut && *tagJ2_pt > J2_cut ) {
			if ( *boosted_mjj_tag > Mjj_cut ) {
                	HIST[0]->Fill(*fatJ_pt,*weight);
			HIST[1]->Fill(*fatJ_eta,*weight);
                	HIST[2]->Fill(*fatJ_phi,*weight);
                	HIST[3]->Fill(*fatJ_m,*weight);
                	HIST[4]->Fill(abs(*tagJ1_eta - *tagJ2_eta),*weight);
			HIST[6]->Fill(*lvJmass,*weight);
			HIST[7]->Fill(*boosted_mjj_tag,*weight);
			HIST[8]->Fill(*tagJ1_pt,*weight);
			HIST[9]->Fill(*tagJ2_pt,*weight);
//			HIST[10]->Fill(*dphi_MTM,*weight);
//	        	HIST[11]->Fill(*dphi_MET_Jets,*weight);
//	        	HIST[12]->Fill(*dphi_MET_Vhad,*weight);
//        		HIST[13]->Fill(*dphi_lep_MET,*weight);
//        		HIST[14]->Fill(*W_MT,*weight);
			HIST[15]->Fill(*met,*weight);
			HIST[19]->Fill(*BDT_merged, *weight);
			HIST[21]->Fill(*lep_eta,    *weight);
        		HIST[22]->Fill(*Njets,      *weight);
        		HIST[24]->Fill(*boosted_xiV,*weight);
			HIST[32]->Fill(*lvJjjmass,  *weight);	
                	}
                	HIST[5]->Fill((v1+v2).M(),*weight);
                	}
                	}
                	}
			}
                }
		else {
			if ( *TWOSIGNALJETS==1){
			if ( *sigJ1_pt > J1_cut && *sigJ2_pt > J2_cut ) {
			if ( ((*sigJJ_m>64.0 && *sigJJ_m<106.0)&&(region!=2))|| 
			     ((*sigJJ_m<64.0 && *sigJJ_m>106.0)&&(region=2)) ) {
                	if ( *res_tagJ1_pt > 30.0 && *res_tagJ2_pt > 30.0 ) {
                	if ( *res_mjj_tag > Mjj_cut ) {
                	HIST[4]->Fill(abs(*res_tagJ1_eta - *res_tagJ2_eta),*weight);
                	HIST[7]->Fill(*res_mjj_tag,*weight);
                	HIST[8]->Fill(*res_tagJ1_pt,*weight);
                	HIST[9]->Fill(*res_tagJ2_pt,*weight);
//                	HIST[10]->Fill(*dphi_MTM,*weight);
//                	HIST[11]->Fill(*res_dphi_MET_Jets,*weight);
//                	HIST[12]->Fill(*res_dphi_MET_Vhad,*weight);
//                	HIST[13]->Fill(*dphi_lep_MET,*weight);
//                	HIST[14]->Fill(*W_MT,*weight);
                	HIST[15]->Fill(*met,*weight);
                	HIST[16]->Fill(*sigJ1_pt,*weight);
                	HIST[17]->Fill(*sigJ2_pt,*weight);
                	HIST[18]->Fill(*lvjjmass,*weight);
			HIST[20]->Fill(*BDT_resolved,*weight);
			HIST[21]->Fill(*lep_eta,    *weight);
                        HIST[22]->Fill(*Njets,      *weight);
			HIST[23]->Fill(*Ntrkjets,   *weight);
                        HIST[24]->Fill(*res_xiV,    *weight);
			HIST[25]->Fill(*res_tagJ1_width,*weight);
        		HIST[26]->Fill(*res_tagJ2_width,*weight);
        		HIST[27]->Fill(*sigJ1_width,*weight);
        		HIST[28]->Fill(*sigJ2_width,*weight);
//        		HIST[29]->Fill(*res_deta_jj,*weight);
        		HIST[30]->Fill(*res_dr_lv,  *weight);
			HIST[31]->Fill(*lv4jmass,   *weight);
                	}
                	HIST[5]->Fill((v1+v2).M(),*weight);
                	}
			}
			}
		
		}
		}

                }
                }
                }


        }
}

void normal(TH1D* hh0,TH1D* hh1,TH1D* hh2,TH1D* hh3,TH1D* hh4,TH1D* hh5,TLegend* leg, float scale, int rebin)
{
     	hh0->SetStats(0);
	TH1D *hsi   = (TH1D*) hh0->Clone("hsi");
	TH1D *hbk   = (TH1D*) hh1->Clone("hbk");
        hbk->Add(hh2);
        hbk->Add(hh3);
	hbk->Add(hh4);
	hbk->Add(hh5);
	hsi->SetLineColor(kRed);
	hsi->SetStats(0);
	hsi->Rebin(rebin);
	hbk->Rebin(rebin);
	hsi->SetTitle("");
        leg->SetBorderSize(0);
        leg->SetFillStyle (0);
        leg->SetTextSize  (0.04);
        leg->AddEntry(hsi, "Signal",    "lp");
        leg->AddEntry(hbk, "Background",   "lp");
	hsi->Scale(1.0 / hsi->Integral());
        hbk->Scale(1.0 / hbk->Integral());
	hsi->SetMaximum(0.3);
	hsi->Draw("hist same");
	hbk->Draw("hist same");
	leg->Draw();
}

void stack_plot(TH1D* h0,TH1D* h1,TH1D* h2,TH1D* h3,TH1D* h4,TH1D* h5,TLegend* leg, float scale, int rebin)
{
	h0->SetStats(0);
	TH1D *hh0   = (TH1D*) h0->Clone("hh0");
	TH1D *hh1   = (TH1D*) h1->Clone("hh1");
	TH1D *hh2   = (TH1D*) h2->Clone("hh2");
        TH1D *hh3   = (TH1D*) h3->Clone("hh3");
	TH1D *hh4   = (TH1D*) h4->Clone("hh4");
        TH1D *hh5   = (TH1D*) h5->Clone("hh5");
	
        THStack *hs = new THStack("hs","");
        hh1->SetFillColor(kGreen-10);
        hh2->SetFillColor(kViolet-9);
	hh3->SetFillColor(kBlue);
        hh4->SetFillColor(kOrange);
	hh5->SetFillColor(kPink);

	hh0->Rebin(rebin);
	hh1->Rebin(rebin);
	hh2->Rebin(rebin);
        hh3->Rebin(rebin);
	hh4->Rebin(rebin);
        hh5->Rebin(rebin);

	hh2->Scale(0.613);

	hs->Add(hh1);
        hs->Add(hh2);
	hs->Add(hh3);
	hs->Add(hh4);
        hs->Add(hh5);

	hh0->SetMarkerStyle(kFullCircle);
        hh0->SetLineColor(1);
        hh0->SetTitle("");
	hh0->SetLineColor(1);
        hh0->SetMarkerColor(1);
        hh0->SetMarkerSize(1);
	hh0->GetYaxis()->SetTitle("#Events");
        leg->SetBorderSize(0);
        leg->SetFillStyle (0);
        leg->SetTextSize  (0.04);
	leg->AddEntry(hh0, "Signal_MC",    "lp");
//        leg->AddEntry(hh0, "data",    "lp");
        leg->AddEntry(hh1, "t#bar{t}",  "F");
	leg->AddEntry(hh2, "W+jets",    "F");
        leg->AddEntry(hh3, "di-boson",  "F");
	leg->AddEntry(hh4, "Z+jets",    "F");
        leg->AddEntry(hh5, "singletop", "F");
//        hsi->Scale(1.0 / hsi->Integral());
//        hbk->Scale(1.0 / hbk->Integral());
	hh0->SetMaximum(10000.0);
//	hh0->SetMaximum((hh0->GetMaximum())*scale);
        hh0->Draw("e   same");
	hs->Draw("hist same");
        hh0->Draw("e   same");
        leg->Draw();
	gPad->SetLogy();
	gPad->RedrawAxis();
}


void t_normal(TH1D* hh0,TH1D* hh1, TLegend* leg, float scale, int rebin)
{
        hh0->SetStats(0);
        TH1D *hsi   = (TH1D*) hh0->Clone("hsi");
        TH1D *hbk   = (TH1D*) hh1->Clone("hbk");
        hsi->SetLineColor(kRed);
        hsi->SetStats(0);
        hsi->Rebin(rebin);
        hbk->Rebin(rebin);
        hsi->SetTitle("");
        leg->SetBorderSize(0);
        leg->SetFillStyle (0);
        leg->SetTextSize  (0.04);
//	leg->SetHeader("ttbar","C");
//        leg->AddEntry(hsi, "POWHEG+Pythia_8", "lp");
//	leg->AddEntry(hsi, "POWHEG+Herwig",   "lp");
//	leg->AddEntry(hbk, "radLo",  "lp");
//	leg->AddEntry(hbk, "POWHEG+Herwig",   "lp");
//        leg->AddEntry(hbk, "aMC@NLO+Herwig",  "lp");
//	leg->AddEntry(hbk, "radHi",  "lp");
//	leg->SetHeader("W+jets","C");
	leg->SetHeader("Diboson","C");
        leg->AddEntry(hsi, "Sherpa",     "lp");
	leg->AddEntry(hbk, "POWHEG+Pythia_8", "lp");
//        leg->AddEntry(hbk, "MadGraph",   "lp");
	if ( (hbk->Integral())>0 && (hsi->Integral())>0 ){
        hsi->Scale(1.0 / hsi->Integral());
	hbk->Scale(1.0 / hbk->Integral());
//        hbk->Scale(hsi->Integral() / hbk->Integral());
//	hsi->SetMaximum((hsi->GetMaximum())*10);
//        hsi->SetMaximum(0.3);
	}
	hsi->SetMaximum((hsi->GetMaximum())*2.2);
	hsi->SetMinimum(0.0);
        hsi->Draw("same");
        hbk->Draw("same");
        leg->Draw();
}

void AddRatioPlot(TH1D* hh0,TH1D* hh1,TH1D* h2,TH1D* hh3,TH1D* hh4,TH1D* hh5, int rebin, TString y_name)
{

	TH1D *hh2 = (TH1D*) h2->Clone("hh2");
	hh2->Scale(0.613);

	TH1D *htest = (TH1D*) hh0->Clone("htest");
        TH1D *hmc   = (TH1D*) hh1->Clone("hmc");
        hmc->Add(hh2);
        hmc->Add(hh3);
	hmc->Add(hh4);
	hmc->Add(hh5);	
	htest->Rebin(rebin);
        hmc->Rebin(rebin);
        htest->SetTitle("");
	htest->Sumw2();
//	htest->Scale(1.0 / htest->Integral());
//        hmc->Scale(1.0 / hmc->Integral());
        htest->Divide(hmc);
        htest->SetLineColor(1);
        htest->SetMarkerColor(1);
        htest->SetMarkerSize(1);
        htest->GetYaxis()->SetTitle(y_name);
        htest->GetYaxis()->SetTitleSize(.1);
        htest->GetYaxis()->SetTitleOffset(0.35);
        htest->GetYaxis()->SetLabelSize(0.1);
        htest->GetYaxis()->SetRangeUser(0.0,0.2);
        htest->GetXaxis()->SetLabelSize(0.1);
        htest->GetXaxis()->SetTitleSize(0.1);
        htest->GetXaxis()->SetTitleOffset(1);
        htest->Draw("e0 same");
	Int_t nBins          = htest->GetNbinsX();
  	Double_t xLowerLimit = htest->GetBinLowEdge(1);
  	Double_t xUpperLimit = htest->GetBinLowEdge(nBins) + htest->GetBinWidth(nBins);
	TLine *line = new TLine(xLowerLimit, 1., xUpperLimit, 1.);
   	line->SetLineColor(kBlack);
   	line->Draw("same");
   	htest->Draw("e0 same");
}

void AddLatexFun(Double_t p0, Double_t e0, Double_t p1, Double_t e1)
{       
        TString t_p0 = "y = ";
	t_p0 += std::to_string(p1); 
	t_p0 += "(#pm";
	t_p0 += std::to_string(e1);
	t_p0 += ")x#plus";
	t_p0 += std::to_string(p0);
	t_p0 += "(#pm";
	t_p0 += std::to_string(e0);
	t_p0 += ")";
        TLatex latex1;
        latex1.SetNDC();
	latex1.SetTextSize(0.12);
        latex1.SetTextFont(42);
        latex1.DrawLatex(0.35, 0.8, t_p0);
//        latex1.DrawLatex(0.68, 0.78, "Internal");
//        //        latex1.SetTextFont (72);
//        //        latex1.DrawLatex(0.55, 0.78, "ATLAS");

}

void s_AddRatioPlot(TH1D* hh0,TH1D* hh1, int rebin, TString r_name="Syst/Norminal")
{
        TH1D *htest = (TH1D*) hh1->Clone("htest");
        TH1D *hmc   = (TH1D*) hh0->Clone("hmc");
	htest->Rebin(rebin);
        hmc->Rebin(rebin);
        htest->SetTitle("");
	htest->SetStats(0);
	if ( (htest->Integral())>0 && (hmc->Integral())>0 ){
	htest->Scale(1.0 / htest->Integral());
        hmc->Scale(1.0 / hmc->Integral());
	htest->Divide(hmc);
	}
        htest->SetLineColor(1);
        htest->SetMarkerColor(1);
        htest->SetMarkerSize(1);
        htest->GetYaxis()->SetTitle(r_name);
        htest->GetYaxis()->SetTitleSize(.1);
        htest->GetYaxis()->SetTitleOffset(0.35);
        htest->GetYaxis()->SetLabelSize(0.1);
        htest->GetYaxis()->SetRangeUser(0.02,1.98);
        htest->GetXaxis()->SetLabelSize(0.1);
        htest->GetXaxis()->SetTitleSize(0.1);
        htest->GetXaxis()->SetTitleOffset(1);
        htest->Draw("e0 same");
	if ( (htest->Integral())>0 && (hmc->Integral())>0 ){
	htest->Fit("pol1");
	TF1 *fit = htest->GetFunction("pol1");
        Double_t p1 = fit->GetParameter(0);
        Double_t e1 = fit->GetParError(0);
        Double_t p2 = fit->GetParameter(1);
        Double_t e2 = fit->GetParError(1);
	AddLatexFun(p1,e1,p2,e2);
	}
        Int_t nBins          = htest->GetNbinsX();
        Double_t xLowerLimit = htest->GetBinLowEdge(1);
        Double_t xUpperLimit = htest->GetBinLowEdge(nBins) + htest->GetBinWidth(nBins);
        TLine *line = new TLine(xLowerLimit, 1., xUpperLimit, 1.);
        line->SetLineColor(kBlack);
        line->Draw("same");
        htest->Draw("e0 same");
}

void Multi_j(TH1D* hh0,int rebin)
{
	TLegend* leg = new TLegend(0.55, 0.45, 0.95, 0.70);
	hh0->SetStats(0);
	hh0->Rebin(rebin);
	hh0->GetYaxis()->SetLabelSize(0.1);
	hh0->GetXaxis()->SetLabelSize(0.1);
        hh0->GetXaxis()->SetTitleSize(0.1);
	hh0->DrawNormalized("hist");
	leg->SetBorderSize(0);
        leg->SetFillStyle (0);
        leg->SetTextSize  (0.1);
        leg->AddEntry(hh0, "Multi_jets",    "lp");
	leg->Draw();
}

void CreateSubPad(TCanvas *canvas, Double_t vfrac=0.25)
{       
        canvas->SetCanvasSize(canvas->GetWw(),(1.+vfrac)*canvas->GetWh());
	canvas->SetWindowSize(4+canvas->GetWw(),50+canvas->GetWh());	

        Double_t xlow, ylow, xup, yup;
        canvas->GetPad(0)->GetPadPar(xlow,ylow,xup,yup);
        
        canvas->Divide(1,3);
        
        TVirtualPad *upPad = canvas->GetPad(1);
        upPad->SetPad(xlow,yup-0.51*(yup-ylow),xup,yup);
        
	TVirtualPad *midPad = canvas->GetPad(2);
        midPad->SetPad(xlow,ylow+0.25*(yup-ylow),xup,yup-0.51*(yup-ylow));

        TVirtualPad *dwPad = canvas->GetPad(3);
        dwPad->SetPad(xlow,ylow,xup,ylow+0.25*(yup-ylow));
        
        canvas->Update();
        return;
}

void Create2Pad(TCanvas *canvas, Double_t vfrac=0.25)
{
        canvas->SetCanvasSize(canvas->GetWw(),(1.+vfrac)*canvas->GetWh());
        canvas->SetWindowSize(4+canvas->GetWw(),50+canvas->GetWh());

        Double_t xlow, ylow, xup, yup;
        canvas->GetPad(0)->GetPadPar(xlow,ylow,xup,yup);

        canvas->Divide(1,2);

        TVirtualPad *upPad = canvas->GetPad(1);
        upPad->SetPad(xlow,ylow+vfrac*(yup-ylow),xup,yup);

        TVirtualPad *dwPad = canvas->GetPad(2);
        dwPad->SetPad(xlow,ylow,xup,ylow+vfrac*(yup-ylow));

        canvas->Update();
        return;
}

void AddLatex()
{
        TLatex latex1;
        latex1.SetNDC();
        latex1.SetTextFont (42);
        latex1.DrawLatex(0.55, 0.83, "#sqrt{s} = 13 TeV, 36.1 fb^{-1}");
        latex1.DrawLatex(0.68, 0.78, "Internal");
        latex1.SetTextFont (72);
        latex1.DrawLatex(0.55, 0.78, "ATLAS");

}

void set_hist(TH1D *HIST[40])
{
//      const double * mjj_bins[7] = {NULL};
	double BDT_merged_bins[10] = {-1,-0.30,-0.18,-0.06,0.06,0.18,0.30,0.42,0.46,0.5};
	double BDT_resolved_bins[24] = {-1,-0.70,-0.64,-0.58,-0.52,-0.46,-0.40,-0.34,-0.28,-0.22,-0.16,-0.10,-0.04,0.02,0.08,0.14,0.20,0.26,0.32,0.38,0.44,0.50,0.56,0.64};
	double mjj_bins[8] = {0,400,600,1000,1400,1800,2500,4000};
	double lvJjj_m_bins[9]   = {0,1000,1400,1600,2000,2500,3000,4000,5000};
	double lvjjjj_m_bins[28] = {0,1000,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3200,3400,3600,3800,4000,4500,5000};	
//	double mjj_bin[] = {400.0,600.0,1000.0,1400.0,1800.0,2500.0,4000.0};
//	const double * mjj_bins[8] = {0.0,400.0,600.0,1000.0,1400.0,1800.0,2500.0,4000.0};
	HIST[0] = new TH1D("fatJ_pt",    "", 40, 0,  2000);
        HIST[1] = new TH1D("fatJ_eta",   "", 60, -3, 3);
        HIST[2] = new TH1D("fatJ_phi",   "", 100, -5, 5);
        HIST[3] = new TH1D("fatJ_m",     "", 50,  0,  250);
        HIST[4] = new TH1D("d_eta",      "", 100, 0,  10);
        HIST[5] = new TH1D("MJJ",        "", 500, 0,  5000);
        HIST[6] = new TH1D("lvJ_m",      "", 30,  0,  3000);
        HIST[7] = new TH1D("mjj_tag",    "", 7,   mjj_bins); //var bin
        HIST[8] = new TH1D("tagj1_pt",   "", 60,  0,  600);
        HIST[9] = new TH1D("tagj2_pt",   "", 50,  0,  500);
	HIST[10]= new TH1D("dphi_MTM",      "", 70, 0, 3.5);
        HIST[11]= new TH1D("dphi_MET_Jets", "", 70, 0, 3.5);
        HIST[12]= new TH1D("dphi_MET_Vhad", "", 70, 0, 3.5);
        HIST[13]= new TH1D("dphi_lep_MET",  "", 70, 0, 3.5);
        HIST[14]= new TH1D("W_MT",          "", 100, 0, 500);
	HIST[15]= new TH1D("MET",           "", 120, 0, 600);
	HIST[16] = new TH1D("sigj1_pt",     "", 50,  0, 500);
        HIST[17] = new TH1D("sigj2_pt",     "", 50,  0, 500);
	HIST[18] = new TH1D("lvjj_m",       "", 30,  0, 3000);
	HIST[19] = new TH1D("BDT_merged",        "", 9,   BDT_merged_bins);
	HIST[20] = new TH1D("BDT_resolved",      "", 23,  BDT_resolved_bins);
	HIST[21] = new TH1D("lep_eta",           "", 60,  -3, 3);
	HIST[22] = new TH1D("Njets",             "", 10,   0, 10);
        HIST[23] = new TH1D("Ntrkjets",          "", 10,   0, 10);
        HIST[24] = new TH1D("xiV",               "", 20,  -5, 5);
        HIST[25] = new TH1D("tagj1_width",       "", 30,  0, 0.3);
        HIST[26] = new TH1D("tagj2_width",       "", 30,  0, 0.3);
        HIST[27] = new TH1D("sigj1_width",       "", 30,  0, 0.3);
        HIST[28] = new TH1D("sigj2_width",       "", 30,  0, 0.3);
        HIST[29] = new TH1D("res_deta_jj",       "", 12,  0, 6);
	HIST[30] = new TH1D("res_dr_lv",         "", 12,  0, 6);
	HIST[31] = new TH1D("lvjjjj_m",          "", 27,  lvjjjj_m_bins);
	HIST[32] = new TH1D("lvJjj_m",           "", 8,   lvJjj_m_bins);
}

void plot_and_save(TH1D* hh0,TH1D* hh1,TH1D* hh2,TH1D* hh3,TH1D* hh4,TH1D* hh5,TH1D* hh6, int rebin, TString name, TString unit)
{
	std::string output_p2("/afs/cern.ch/work/j/jzeng/public/VBS_test/output_HP/");
//	output_p2 += "test";
	output_p2 += name;
//        output_p2 += "tagj2_";
//        output_p2 += TString::Itoa(name_p1,10);
//        output_p2 += "_";
//        output_p2 += TString::Itoa(name_p2,10);
        output_p2 += ".png";

        TLegend* leg1 = new TLegend(0.55, 0.45, 0.95, 0.70);
        TCanvas *c1   = new TCanvas("c1","",900,600);
        CreateSubPad(c1);
        c1->cd(1);
//        hh0->GetXaxis()->SetTitle("Mjj_tag[Gev]");
	hh0->GetXaxis()->SetTitle(unit);
	hh0->GetXaxis()->SetLabelSize(0.045);
        hh0->GetXaxis()->SetTitleSize(0.045);
        normal(hh0,hh1,hh2,hh3,hh4,hh5, leg1, 1.5, rebin);
        AddLatex();
        c1->cd(2);
        AddRatioPlot(hh0,hh1,hh2,hh3,hh4,hh5,rebin,"");
	c1->cd(3);
	Multi_j(hh6,rebin);
        c1->SaveAs(output_p2.c_str());
}

//void splot_and_save(TH1D* hh0,TH1D* hh1, int rebin, TString name, TString unit)
//{
////        std::string output_p2("/afs/cern.ch/work/j/jzeng/public/VBS_test/output_test/");
//	std::string output_p2("/afs/cern.ch/work/j/jzeng/public/VBS_test/output_W+jets/");
//	output_p2 += name;
//	output_p2 += ".png";
//
//        TLegend* leg1 = new TLegend(0.55, 0.45, 0.95, 0.70);
//        TCanvas *c1   = new TCanvas("c1","",900,600);
//        Create2Pad(c1);
//        c1->cd(1);
//	hh0->GetXaxis()->SetTitle(unit);
//        hh0->GetXaxis()->SetLabelSize(0.045);
//        hh0->GetXaxis()->SetTitleSize(0.045);
//        s_normal(hh0,hh1, leg1, 1.5, rebin);
//        AddLatex();
//        c1->cd(2);
//        s_AddRatioPlot(hh0,hh1, rebin);
//	c1->SaveAs(output_p2.c_str());
//}

void tplot_and_save(TH1D* hh0,TH1D* hh1, int rebin, TString r_name, TString name, TString unit)
{
//        std::string output_p2("/afs/cern.ch/work/j/jzeng/public/VBS_test/output_test/");
	std::string output_p2("/afs/cern.ch/work/j/jzeng/public/VBS_test/output_ttbar/");
	output_p2 += name;
	output_p2 += ".png";

        TLegend* leg1 = new TLegend(0.55, 0.45, 0.95, 0.70);
        TCanvas *c1   = new TCanvas("c1","",900,600);
        Create2Pad(c1);
        c1->cd(1);
	hh0->GetXaxis()->SetTitle(unit);
        hh0->GetXaxis()->SetLabelSize(0.045);
        hh0->GetXaxis()->SetTitleSize(0.045);
        t_normal(hh0,hh1, leg1, 1.5, rebin);
        AddLatex();
        c1->cd(2);
        s_AddRatioPlot(hh0,hh1, rebin,r_name);
	c1->SaveAs(output_p2.c_str());
}



void top_and_save(TH1D* hh0,TH1D* hh1,TH1D* hh2,TH1D* hh3,TH1D* hh4,TH1D* hh5,TH1D* hh6, int rebin, TString name, TString unit)
{
        std::string output_p2("/afs/cern.ch/work/j/jzeng/public/VBS_test/output_test/");
        output_p2 += name;
        output_p2 += ".png";

        TLegend* leg1 = new TLegend(0.65, 0.48, 0.95, 0.73);
        TCanvas *c1   = new TCanvas("c1","",900,600);
        Create2Pad(c1);
        c1->cd(1);
        hh0->GetXaxis()->SetTitle(unit);
//        hh0->GetXaxis()->SetLabelSize(0.045);
//        hh0->GetXaxis()->SetTitleSize(0.045);
        stack_plot(hh0,hh1,hh2,hh3,hh4,hh5, leg1, 1.5, rebin);
        AddLatex();
        c1->cd(2);
	AddRatioPlot(hh0,hh1,hh2,hh3,hh4,hh5,rebin,"Signal/Bg");
        c1->SaveAs(output_p2.c_str());
}

//};





