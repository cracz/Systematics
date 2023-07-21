#include <iostream>
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "Variation.h"
#include "PlotUtils.h"

using namespace PlotUtils;

// Constructor for the "normal" results.
Variation::Variation(TString prefix, TString order_n_str)
{
  ID = prefix;
  fileName = prefix + ".picoDst.result.combined.root";

  file = TFile::Open(fileName);

  initialize(order_n_str);
  fixAttributes(order_n_str);
}


// Destructor
Variation::~Variation()
{
  delete h_vn_pp;
  delete h_vn_pm;
  delete h_vn_kp;
  delete h_vn_km;
  delete h_vn_pr;
  delete h_vn_yCM_00to10_pr;
  delete h_vn_yCM_10to40_pr;
  delete h_vn_yCM_40to60_pr;
  delete h_vn_yCM_00to10_pr_symm;
  delete h_vn_yCM_10to40_pr_symm;
  delete h_vn_yCM_40to60_pr_symm;
  delete h_vn_pT_00to10_pr;
  delete h_vn_pT_10to40_pr;
  delete h_vn_pT_40to60_pr;
  file->Close();
}

// Initialize all histograms
void Variation::initialize(TString order_n_str)
{
  TProfile *p_vn_pp = (TProfile*)file->Get("p_vn_pp");
  TProfile *p_vn_pm = (TProfile*)file->Get("p_vn_pm");
  TProfile *p_vn_kp = (TProfile*)file->Get("p_vn_kp");
  TProfile *p_vn_km = (TProfile*)file->Get("p_vn_km");
  TProfile *p_vn_pr = (TProfile*)file->Get("p_vn_pr");
  p_vn_kp->Rebin();
  p_vn_km->Rebin();
  
  h_vn_pp = p_vn_pp->ProjectionX((TString)p_vn_pp->GetName() +"_"+ ID);
  h_vn_pm = p_vn_pm->ProjectionX((TString)p_vn_pm->GetName() +"_"+ ID);
  h_vn_kp = p_vn_kp->ProjectionX((TString)p_vn_kp->GetName() +"_"+ ID);
  h_vn_km = p_vn_km->ProjectionX((TString)p_vn_km->GetName() +"_"+ ID);
  h_vn_pr = p_vn_pr->ProjectionX((TString)p_vn_pr->GetName() +"_"+ ID);
  
  h_vn_pp = flipHisto(h_vn_pp);
  h_vn_pm = flipHisto(h_vn_pm);
  h_vn_kp = flipHisto(h_vn_kp);
  h_vn_km = flipHisto(h_vn_km);
  h_vn_pr = flipHisto(h_vn_pr);
  
  h_vn_pp = trimCentralityPlot(h_vn_pp);
  h_vn_pm = trimCentralityPlot(h_vn_pm);
  h_vn_kp = trimCentralityPlot(h_vn_kp);
  h_vn_km = trimCentralityPlot(h_vn_km);
  h_vn_pr = trimCentralityPlot(h_vn_pr);
  
  
  //=== vn vs rapidity stuff
  TProfile2D *p2_vn_yCM_cent_pr = (TProfile2D*)file->Get("p2_vn_yCM_cent_pr");
  TProfile2D *p2_vn_yCM_cent_pr_symmetry = (TProfile2D*)file->Get("p2_vn_yCM_cent_pr_symmetry");

  TProfile *p_vn_yCM_00to10_pr = p2_vn_yCM_cent_pr->ProfileY("p_vn_yCM_00to10_pr", 15, 16);
  TProfile *p_vn_yCM_10to40_pr = p2_vn_yCM_cent_pr->ProfileY("p_vn_yCM_10to40_pr", 9, 14);
  TProfile *p_vn_yCM_40to60_pr = p2_vn_yCM_cent_pr->ProfileY("p_vn_yCM_40to60_pr", 5, 8);

  TProfile *p_vn_yCM_00to10_pr_symm = p2_vn_yCM_cent_pr_symmetry->ProfileY("p_vn_yCM_00to10_pr_symm", 15, 16);
  TProfile *p_vn_yCM_10to40_pr_symm = p2_vn_yCM_cent_pr_symmetry->ProfileY("p_vn_yCM_10to40_pr_symm", 9, 14);
  TProfile *p_vn_yCM_40to60_pr_symm = p2_vn_yCM_cent_pr_symmetry->ProfileY("p_vn_yCM_40to60_pr_symm", 5, 8);

  h_vn_yCM_00to10_pr = new TH1D("h_vn_yCM_00to10_pr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  h_vn_yCM_10to40_pr = new TH1D("h_vn_yCM_10to40_pr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  h_vn_yCM_40to60_pr = new TH1D("h_vn_yCM_40to60_pr", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  h_vn_yCM_00to10_pr_symm = new TH1D("h_vn_yCM_00to10_pr_symm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  h_vn_yCM_10to40_pr_symm = new TH1D("h_vn_yCM_10to40_pr_symm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);
  h_vn_yCM_40to60_pr_symm = new TH1D("h_vn_yCM_40to60_pr_symm", ";y-y_{mid};v_{"+order_n_str+"}", 20, -1, 1);

  
  // Convert profiles to histograms
  h_vn_yCM_00to10_pr = p_vn_yCM_00to10_pr->ProjectionX();
  h_vn_yCM_10to40_pr = p_vn_yCM_10to40_pr->ProjectionX();
  h_vn_yCM_40to60_pr = p_vn_yCM_40to60_pr->ProjectionX();

  h_vn_yCM_00to10_pr_symm = p_vn_yCM_00to10_pr_symm->ProjectionX();
  h_vn_yCM_10to40_pr_symm = p_vn_yCM_10to40_pr_symm->ProjectionX();
  h_vn_yCM_40to60_pr_symm = p_vn_yCM_40to60_pr_symm->ProjectionX();



  //=== vn vs pT stuff
  TProfile2D *p2_vn_pT_cent_pr = (TProfile2D*)file->Get("p2_vn_pT_cent_pr");

  TProfile *p_vn_pT_00to10_pr = p2_vn_pT_cent_pr->ProfileY("p_vn_pT_00to10_pr", 15, 16);
  TProfile *p_vn_pT_10to40_pr = p2_vn_pT_cent_pr->ProfileY("p_vn_pT_10to40_pr", 9, 14);
  TProfile *p_vn_pT_40to60_pr = p2_vn_pT_cent_pr->ProfileY("p_vn_pT_40to60_pr", 5, 8);

  h_vn_pT_00to10_pr = new TH1D("h_vn_pT_00to10_pr", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  h_vn_pT_10to40_pr = new TH1D("h_vn_pT_10to40_pr", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  h_vn_pT_40to60_pr = new TH1D("h_vn_pT_40to60_pr", ";p_{T} (GeV);v_{"+order_n_str+"}", 10, 0, 2);
  
  h_vn_pT_00to10_pr = p_vn_pT_00to10_pr->ProjectionX();
  h_vn_pT_10to40_pr = p_vn_pT_10to40_pr->ProjectionX();
  h_vn_pT_40to60_pr = p_vn_pT_40to60_pr->ProjectionX();

}// End initialize()



// Use this with the "Normal" Variation to clean up the plots
void Variation::fixAttributes(TString order_n_str)
{
  if (order_n_str == "3")
    {
      Double_t centralityUpperBounds = 0.15;
      Double_t centralityLowerBounds = -0.15;
      
      h_vn_pp->SetMarkerStyle(20);
      h_vn_pp->SetMarkerSize(2.5);
      h_vn_pp->SetMarkerColor(2);
      h_vn_pp->SetLineColor(2);
      h_vn_pp->SetLineWidth(3);
      h_vn_pp->GetYaxis()->SetTitleOffset(1.7);
      h_vn_pp->GetXaxis()->SetNdivisions(210);
      //h_vn_pp->SetMaximum(centralityUpperBounds);
      //h_vn_pp->SetMinimum(centralityLowerBounds);

      h_vn_pm->SetMarkerStyle(20);
      h_vn_pm->SetMarkerSize(2.5);
      h_vn_pm->SetMarkerColor(4);
      h_vn_pm->SetLineColor(4);
      h_vn_pm->SetLineWidth(3);
      h_vn_pm->GetYaxis()->SetTitleOffset(1.7);
      h_vn_pm->GetXaxis()->SetNdivisions(210);
      //h_vn_pm->SetMaximum(centralityUpperBounds);
      //h_vn_pm->SetMinimum(centralityLowerBounds);

      h_vn_kp->SetMarkerStyle(20);
      h_vn_kp->SetMarkerSize(2.5);
      h_vn_kp->SetMarkerColor(2);
      h_vn_kp->SetLineColor(2);
      h_vn_kp->SetLineWidth(3);
      h_vn_kp->GetYaxis()->SetTitleOffset(1.7);
      h_vn_kp->GetXaxis()->SetNdivisions(210);
      //h_vn_kp->SetMaximum(centralityUpperBounds);
      //h_vn_kp->SetMinimum(centralityLowerBounds);

      h_vn_km->SetMarkerStyle(20);
      h_vn_km->SetMarkerSize(2.5);
      h_vn_km->SetMarkerColor(4);
      h_vn_km->SetLineColor(4);
      h_vn_km->SetLineWidth(3);
      h_vn_km->GetYaxis()->SetTitleOffset(1.7);
      h_vn_km->GetXaxis()->SetNdivisions(210);
      //h_vn_km->SetMaximum(centralityUpperBounds);
      //h_vn_km->SetMinimum(centralityLowerBounds);

  
      h_vn_pr->SetTitle(";Centrality (%);v_{"+order_n_str+"}");
      h_vn_pr->SetMarkerStyle(20);
      h_vn_pr->SetMarkerSize(2.5);
      h_vn_pr->SetMarkerColor(2);
      h_vn_pr->SetLineColor(2);
      h_vn_pr->SetLineWidth(3);
      h_vn_pr->GetYaxis()->SetTitleOffset(1.7);
      h_vn_pr->GetXaxis()->SetNdivisions(210);
      //h_vn_pr->SetMaximum(centralityUpperBounds);
      //h_vn_pr->SetMinimum(centralityLowerBounds);


      //=== vn vs rapidity
      h_vn_yCM_00to10_pr->SetMarkerStyle(20);
      h_vn_yCM_10to40_pr->SetMarkerStyle(20);
      h_vn_yCM_40to60_pr->SetMarkerStyle(20);
      h_vn_yCM_00to10_pr->SetMarkerColor(2);
      h_vn_yCM_10to40_pr->SetMarkerColor(4);
      h_vn_yCM_40to60_pr->SetMarkerColor(8);
      h_vn_yCM_00to10_pr->SetMarkerSize(2);
      h_vn_yCM_10to40_pr->SetMarkerSize(2);
      h_vn_yCM_40to60_pr->SetMarkerSize(2);
      h_vn_yCM_00to10_pr->SetLineColor(2);
      h_vn_yCM_10to40_pr->SetLineColor(4);
      h_vn_yCM_40to60_pr->SetLineColor(8);
      h_vn_yCM_00to10_pr->SetLineWidth(3);
      h_vn_yCM_10to40_pr->SetLineWidth(3);
      h_vn_yCM_40to60_pr->SetLineWidth(3);

      h_vn_yCM_00to10_pr_symm->SetMarkerStyle(20);
      h_vn_yCM_10to40_pr_symm->SetMarkerStyle(20);
      h_vn_yCM_40to60_pr_symm->SetMarkerStyle(20);
      h_vn_yCM_00to10_pr_symm->SetMarkerColor(2);
      h_vn_yCM_10to40_pr_symm->SetMarkerColor(4);
      h_vn_yCM_40to60_pr_symm->SetMarkerColor(8);
      h_vn_yCM_00to10_pr_symm->SetMarkerSize(2);
      h_vn_yCM_10to40_pr_symm->SetMarkerSize(2);
      h_vn_yCM_40to60_pr_symm->SetMarkerSize(2);
      h_vn_yCM_00to10_pr_symm->SetLineColor(2);
      h_vn_yCM_10to40_pr_symm->SetLineColor(4);
      h_vn_yCM_40to60_pr_symm->SetLineColor(8);
      h_vn_yCM_00to10_pr_symm->SetLineWidth(3);
      h_vn_yCM_10to40_pr_symm->SetLineWidth(3);
      h_vn_yCM_40to60_pr_symm->SetLineWidth(3);


      //=== vn vs pT  
      h_vn_pT_00to10_pr->SetMarkerStyle(20);
      h_vn_pT_10to40_pr->SetMarkerStyle(20);
      h_vn_pT_40to60_pr->SetMarkerStyle(20);
      h_vn_pT_00to10_pr->SetMarkerColor(2);
      h_vn_pT_10to40_pr->SetMarkerColor(4);
      h_vn_pT_40to60_pr->SetMarkerColor(8);
      h_vn_pT_00to10_pr->SetMarkerSize(2);
      h_vn_pT_10to40_pr->SetMarkerSize(2);
      h_vn_pT_40to60_pr->SetMarkerSize(2);
      h_vn_pT_00to10_pr->SetLineColor(2);
      h_vn_pT_10to40_pr->SetLineColor(4);
      h_vn_pT_40to60_pr->SetLineColor(8);
      h_vn_pT_00to10_pr->SetLineWidth(3);
      h_vn_pT_10to40_pr->SetLineWidth(3);
      h_vn_pT_40to60_pr->SetLineWidth(3);
    }
}// End fixAttributes()
