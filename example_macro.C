
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include <fstream>
#include <iostream>
#include <regex>
#include <string>

#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"

Double_t zGlobal = 395;

Int_t getMaximumEnergyIndex(TTreeReaderArray<Float_t> &energyVector) {
  Int_t maxIndex = -1;
  Float_t maxEnergy = 0.;
  for (int i = 0; i < energyVector.GetSize(); ++i) {
    if (energyVector[i] > maxEnergy) {
      maxEnergy = energyVector[i];
      maxIndex = i;
    }
  }
  return maxIndex;
}

struct Cluster {
  Cluster(TString _type)
      : theta(0), phi(0), energy(0), x(0), y(0), phiResolution(0),
        thetaResolution(0) {
    type = _type;
  };
  Cluster()
      : theta(0), phi(0), energy(0), x(0), y(0), phiResolution(0),
        thetaResolution(0), type("None"){};

  Cluster addEcalWithSampleCoefficient(const Float_t &sampleCoefficient,
                                       Cluster &emcal) {
    Cluster result;
    result.energy = (sampleCoefficient * emcal.energy) + energy;
    Double_t HcalContribution = energy / result.energy;
    Double_t EcalContribution = 1 - HcalContribution;
    result.theta =
        (EcalContribution * emcal.theta) + (HcalContribution * theta);
    result.phi = (EcalContribution * emcal.phi) + (HcalContribution * phi);

    result.xVertex =
        0; // (EcalContribution * emcal.xVertex) + (HcalContribution * xVertex);
    result.yVertex =
        0; // (EcalContribution * emcal.yVertex) + (HcalContribution * yVertex);

    result.x = (EcalContribution * emcal.x) + (HcalContribution * x);
    result.y = (EcalContribution * emcal.y) + (HcalContribution * y);

    return result;
  }

  Cluster &operator-(const Cluster &mc_position) {
    thetaResolution = (theta - mc_position.theta);
    phiResolution = (phi - mc_position.phi);

    x = zGlobal * tan(theta) * cos(phi);
    y = zGlobal * tan(theta) * sin(phi);

    Double_t xResolution = (x - mc_position.x);
    Double_t yResolution = (y - mc_position.y);

    rResolution = sqrt(xResolution * xResolution + yResolution * yResolution);

    return *this;
  }

  void normalize() {

    if (phi > TMath::Pi())
      phi -= 2 * TMath::Pi();
    if (phi < -TMath::Pi())
      phi += 2 * TMath::Pi();

    x = zGlobal * tan(theta) * cos(phi);
    y = zGlobal * tan(theta) * sin(phi);
  }

  Double_t theta;
  Double_t phi;
  Double_t energy;
  Double_t x;
  Double_t y;

  Double_t xVertex;
  Double_t yVertex;

  Double_t phiResolution;
  Double_t thetaResolution;
  Double_t rResolution;

  TString type;
};
Double_t getMaximum(TH1D *h1, TH1D *h2, TH1D *h3) {
  Double_t max = h1->GetMaximum();
  if (h2->GetMaximum() > max)
    max = h2->GetMaximum();
  if (h3->GetMaximum() > max)
    max = h3->GetMaximum();
  return max;
}
struct ClusterHists {
  ClusterHists(TString _name) {
    hTheta = new TH1D(
        "hTheta_" + _name,
        _name + " Cluster #theta_{cluster}; #theta_{cluster} [rad]; Entries",
        500, 0, 5);
    hPhi = new TH1D(
        "hPhi_" + _name,
        _name + " Cluster #phi_{cluster}; #phi_{cluster} [rad]; Entries", 500,
        -TMath::Pi(), TMath::Pi());
    hThetaResol =
        new TH1D("hThetaResolution_" + _name,
                 _name + " Cluster #theta Resolution ; #Delta#theta; Entries",
                 500, -TMath::Pi(), TMath::Pi());
    hPhiResol =
        new TH1D("hPhiResolution_" + _name,
                 _name + " Cluster #phi Resolution; #Delta#phi; Entries", 500,
                 -TMath::Pi(), TMath::Pi());
    hRxyResol = new TH1D(
        "hRxyResolution_" + _name,
        _name + " Cluster R_{XY} Resolution; #Delta R_{XY}[cm]; Entries", 1000,
        0, 200);
    hEnergy =
        new TH1D("hEnergy_" + _name,
                 _name + " Cluster energy; E [GeV]; Entries", 10000, 0, 10);
    hPos = new TH2D("hPosition_" + _name,
                    _name + " Cluster position x,y; x [mm]; y [mm]; Entries",
                    500, -3000, 3000, 500, -3000, 3000);

    hPhiEnergy = new TH2D("hPhiEnergy_" + _name,
                          _name + " Cluster #phi vs Energy; #phi; Energy", 100,
                          -TMath::Pi(), TMath::Pi(), 100, 0, 10);
    name = _name;
  }

  void Fill(const Cluster &cluster) {
    hTheta->Fill(cluster.theta);
    hPhi->Fill(cluster.phi);
    hEnergy->Fill(cluster.energy);
    hPos->Fill(cluster.xVertex, cluster.yVertex);
    hThetaResol->Fill(cluster.thetaResolution);
    hPhiResol->Fill(cluster.phiResolution);
    hRxyResol->Fill(cluster.rResolution);
  }

  void DrawTogether(TCanvas *can, TString outPdf, ClusterHists &emcalHists,
                    ClusterHists &sumHists, ClusterHists &scatteredHcal) {
    can->cd();
    sumHists.hPhiResol->SetTitle(
        "EMCal, HCal, Sum #phi Resolution; #Delta#phi; Entries");
    sumHists.hPhiResol->SetMarkerColor(kBlack);
    sumHists.hPhiResol->SetLineColor(kBlack);
    sumHists.hPhiResol->GetYaxis()->SetRangeUser(
        0,
        1.2 * getMaximum(sumHists.hPhiResol, emcalHists.hPhiResol, hPhiResol));
    sumHists.hPhiResol->Draw("hist");
    emcalHists.hPhiResol->SetMarkerColor(kBlue);
    emcalHists.hPhiResol->SetLineColor(kBlue);
    emcalHists.hPhiResol->Draw("same");
    hPhiResol->SetMarkerColor(kRed);
    hPhiResol->SetLineColor(kRed);
    hPhiResol->Draw("same");
    scatteredHcal.hPhiResol->SetMarkerColor(kViolet);
    scatteredHcal.hPhiResol->SetLineColor(kViolet);
    scatteredHcal.hPhiResol->Draw("same");

    TLegend leg(0.65, 0.5, 0.8, 0.7);
    leg.AddEntry(emcalHists.hPhiResol, "EMCal only", "l");
    leg.AddEntry(hPhiResol, "HCal only", "l");
    leg.AddEntry(sumHists.hPhiResol, "HCal+EMCal", "l");
    leg.AddEntry(scatteredHcal.hPhiResol, "Scattered HCal", "l");
    leg.Draw();

    can->SaveAs(outPdf);
    sumHists.hThetaResol->SetTitle(
        "EMCal, HCal, Sum #theta Resolution; #Delta#theta; Entries");
    sumHists.hThetaResol->SetMarkerColor(kBlack);
    sumHists.hThetaResol->SetLineColor(kBlack);
    sumHists.hThetaResol->GetYaxis()->SetRangeUser(
        0, 1.2 * getMaximum(sumHists.hThetaResol, emcalHists.hThetaResol,
                            hThetaResol));
    sumHists.hThetaResol->Draw("hist");
    emcalHists.hThetaResol->SetMarkerColor(kBlue);
    emcalHists.hThetaResol->SetLineColor(kBlue);
    emcalHists.hThetaResol->Draw("same");
    hThetaResol->SetMarkerColor(kRed);
    hThetaResol->SetLineColor(kRed);
    hThetaResol->Draw("same");

    scatteredHcal.hThetaResol->SetMarkerColor(kViolet);
    scatteredHcal.hThetaResol->SetLineColor(kViolet);
    scatteredHcal.hThetaResol->Draw("same");

    leg.Draw();
    // TLatex *tl = new TLatex();
    // tl->SetTextSize(0.05);
    // tl->DrawLatexNDC(0.5, 0.5, "No EMCal in");
    can->SaveAs(outPdf);

    sumHists.hTheta->SetTitle("EMCal, HCal, Sum #theta; #theta; Entries");
    sumHists.hTheta->SetMarkerColor(kBlack);
    sumHists.hTheta->SetLineColor(kBlack);
    sumHists.hTheta->GetYaxis()->SetRangeUser(
        0, 1.2 * getMaximum(sumHists.hTheta, emcalHists.hTheta, hTheta));
    sumHists.hTheta->Draw("hist");
    emcalHists.hTheta->SetMarkerColor(kBlue);
    emcalHists.hTheta->SetLineColor(kBlue);

    emcalHists.hTheta->Draw("same");
    hTheta->SetMarkerColor(kRed);
    hTheta->SetLineColor(kRed);
    hTheta->Draw("same");
    scatteredHcal.hTheta->SetMarkerColor(kViolet);
    scatteredHcal.hTheta->SetLineColor(kViolet);
    scatteredHcal.hTheta->Draw("same");

    leg.Draw();
    can->SaveAs(outPdf);

    sumHists.hPhi->SetTitle("EMCal, HCal, Sum  #phi; #phi; Entries");
    sumHists.hPhi->SetMarkerColor(kBlack);
    sumHists.hPhi->SetLineColor(kBlack);
    sumHists.hPhi->GetYaxis()->SetRangeUser(
        0, 1.2 * getMaximum(sumHists.hPhi, emcalHists.hPhi, hPhi));
    sumHists.hPhi->Draw("hist");
    emcalHists.hPhi->SetMarkerColor(kBlue);
    emcalHists.hPhi->SetLineColor(kBlue);
    emcalHists.hPhi->Draw("same");
    hPhi->SetMarkerColor(kRed);
    hPhi->SetLineColor(kRed);
    hPhi->Draw("same");
    scatteredHcal.hPhi->SetMarkerColor(kViolet);
    scatteredHcal.hPhi->SetLineColor(kViolet);
    scatteredHcal.hPhi->Draw("same");

    leg.Draw();
    can->SaveAs(outPdf);

    gPad->SetLogx(1);

    sumHists.hRxyResol->SetTitle(
        "EMCal, HCal, Sum R_{XY} Resolution; #Delta R_{XY}, [mm]; Entries");
    sumHists.hRxyResol->SetMarkerColor(kBlack);
    sumHists.hRxyResol->SetLineColor(kBlack);
    sumHists.hRxyResol->GetYaxis()->SetRangeUser(
        0,
        1.2 * getMaximum(sumHists.hRxyResol, emcalHists.hRxyResol, hRxyResol));
    sumHists.hRxyResol->GetXaxis()->SetRangeUser(0.8, 200);
    sumHists.hRxyResol->Draw("hist");
    emcalHists.hRxyResol->SetMarkerColor(kBlue);
    emcalHists.hRxyResol->SetLineColor(kBlue);

    emcalHists.hRxyResol->Draw("same");
    hRxyResol->SetMarkerColor(kRed);
    hRxyResol->SetLineColor(kRed);
    hRxyResol->Draw("same");

    scatteredHcal.hRxyResol->SetMarkerColor(kViolet);
    scatteredHcal.hRxyResol->SetLineColor(kViolet);
    scatteredHcal.hRxyResol->Draw("same");

    leg.Draw();
    can->SaveAs(outPdf);
    gPad->SetLogx(0);
  }

  Double_t getSigma(TString varname) {
    TH1D *hist;
    TF1 *fitfunc;

    if (varname.Contains("theta")) {
      hist = hThetaResol;
      hist->GetXaxis()->SetRangeUser(-1, 1);

      fitfunc = new TF1("fitfunc", "gaus", -1, 1);
      // fitfunc->SetRange(-0.02, 0.02);
      // fitfunc->SetParLimits(1, -0.1, 0.1);
      // fitfunc->SetParLimits(2, 0, 0.1);
    } else if (varname.Contains("phi")) {
      hist = hPhiResol;
      hist->GetXaxis()->SetRangeUser(-1, 1);
      fitfunc = new TF1("fitfunc", "gaus", -1, 1);
      // fitfunc->SetRange(-0.5, 0.5);
      // fitfunc->SetParLimits(1, -0.5, 0.5);
      // fitfunc->SetParLimits(2, 0, 1);
    }

    else if (varname.Contains("Rxy")) {
      hist = hRxyResol;
      fitfunc = new TF1("fitfunc", "landau", 0, 100);
    }

    else {
      cout << "Wrong variable name" << endl;
      return -1;
    }
    // fitfunc->SetParNames("Area", "Mean", "Sigma");
    fitfunc->SetLineColor(kGreen);

    fitfunc->SetParameters(hist->GetMaximum(), hist->GetMean(),
                           hist->GetStdDev());

    hist->GetYaxis()->SetTitleOffset(1.25);
    if (varname.Contains("Rxy"))
      hist->GetXaxis()->SetRangeUser(0.8, 200);
    hist->SetMarkerStyle(43);
    hist->SetMarkerSize(1.5);
    hist->SetLineColor(kBlue);
    hist->SetMarkerColor(kBlack);

    hist->Fit("fitfunc", "EQMR", "");
    hist->Draw("p");
    Double_t sigma = fitfunc->GetParameter(2);
    Double_t stdDev = hist->GetStdDev(); // sigma * 2.355;

    // Double_t sigma = sigma * 2.355;
    TLatex *tl = new TLatex();
    tl->SetTextSize(0.05);
    tl->DrawLatexNDC(0.62, 0.5, "#" + varname + " " + name);
    tl->DrawLatexNDC(0.62, 0.4, Form("Mean ~ %.2f", hist->GetMean()));
    tl->DrawLatexNDC(0.62, 0.3, Form("Sigma ~ %.2f", sigma));
    tl->DrawLatexNDC(0.62, 0.2, Form("StdDev ~ %.2f", stdDev));

    return stdDev;
  }

  vector<Double_t> Draw(TCanvas *can, TString outPdf) {
    can->cd();
    hTheta->Draw("hist");
    can->SaveAs(outPdf);
    hPhi->Draw("hist");
    can->SaveAs(outPdf);
    hEnergy->Draw("hist");
    can->SaveAs(outPdf);
    hPos->Draw("colz");
    can->SaveAs(outPdf);

    Double_t sigmaThetaSum = getSigma("theta");
    can->SaveAs(outPdf);
    Double_t sigmaPhiSum = getSigma("phi");
    can->SaveAs(outPdf);
    can->SetLogx(1);
    Double_t sigmaRxySum = getSigma("Rxy");
    can->SaveAs(outPdf);
    can->SetLogx(0);

    return {sigmaThetaSum, sigmaPhiSum, sigmaRxySum};
  }

  void Write(TDirectory *output) {
    cout << "Writing histograms" << name << endl;
    output->cd();
    TDirectory *dir = output->mkdir(name);
    dir->cd();
    hTheta->Write();
    hPhi->Write();
    hThetaResol->Write();
    hPhiResol->Write();
    hRxyResol->Write();
    hEnergy->Write();
    hPos->Write();
  }

  TString name;
  TH1D *hTheta;
  TH1D *hPhi;
  TH1D *hThetaResol;
  TH1D *hPhiResol;
  TH1D *hRxyResol;
  TH1D *hEnergy;
  TH2D *hPos;
  TH2D *hPhiEnergy;
};

void example_macro(
    TString inFileName =
        "../eicrecon_neutron_50000events_p5gev_phi48_theta170.24.edm4eic.root",
    TString outFileName = "ana.root") {

  //==========Style of the plot============
  gStyle->SetPalette(1);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleOffset(.85, "X");
  gStyle->SetTitleOffset(.85, "Y");
  gStyle->SetTitleSize(.04, "X");
  gStyle->SetTitleSize(.04, "Y");
  gStyle->SetLabelSize(.04, "X");
  gStyle->SetLabelSize(.04, "Y");
  gStyle->SetHistLineWidth(2);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);

  //=======Reading the root file DD4HEP===========
  TFile *file =
      new TFile(inFileName); // Tree with tracks and hits
                             // Create the tree reader and its data containers
  TTreeReader myReader("events", file); // name of tree and file

  TTreeReaderArray<Float_t> charge(myReader, "MCParticles.charge");
  TTreeReaderArray<Double_t> vx_mc(myReader, "MCParticles.vertex.x");
  TTreeReaderArray<Double_t> vy_mc(myReader, "MCParticles.vertex.y");
  TTreeReaderArray<Double_t> vz_mc(myReader, "MCParticles.vertex.z");
  TTreeReaderArray<Float_t> px_mc(myReader, "MCParticles.momentum.x");
  TTreeReaderArray<Float_t> py_mc(myReader, "MCParticles.momentum.y");
  TTreeReaderArray<Float_t> pz_mc(myReader, "MCParticles.momentum.z");
  TTreeReaderArray<Int_t> status(myReader, "MCParticles.generatorStatus");
  TTreeReaderArray<Int_t> pdg(myReader, "MCParticles.PDG");

  TTreeReaderArray<Float_t> hcal_truth_E(myReader,
                                         "HcalEndcapNTruthClusters.energy");
  TTreeReaderArray<Float_t> hcal_truth_theta(
      myReader, "HcalEndcapNTruthClusters.intrinsicTheta");
  TTreeReaderArray<Float_t> hcal_truth_phi(
      myReader, "HcalEndcapNTruthClusters.intrinsicPhi");
  TTreeReaderArray<Float_t> hcal_truth_x(myReader,
                                         "HcalEndcapNTruthClusters.position.x");
  TTreeReaderArray<Float_t> hcal_truth_y(myReader,
                                         "HcalEndcapNTruthClusters.position.y");

  TTreeReaderArray<Float_t> emcal_truth_E(myReader,
                                          "EcalEndcapNTruthClusters.energy");
  TTreeReaderArray<Float_t> emcal_truth_theta(
      myReader, "EcalEndcapNTruthClusters.intrinsicTheta");
  TTreeReaderArray<Float_t> emcal_truth_phi(
      myReader, "EcalEndcapNTruthClusters.intrinsicPhi");
  TTreeReaderArray<Float_t> emcal_truth_x(
      myReader, "EcalEndcapNTruthClusters.position.x");
  TTreeReaderArray<Float_t> emcal_truth_y(
      myReader, "EcalEndcapNTruthClusters.position.y");

  TTreeReaderArray<Float_t> hcal_reco_E(myReader, "HcalEndcapNClusters.energy");
  TTreeReaderArray<Float_t> hcal_reco_theta(
      myReader, "HcalEndcapNClusters.intrinsicTheta");
  TTreeReaderArray<Float_t> hcal_reco_phi(myReader,
                                          "HcalEndcapNClusters.intrinsicPhi");
  TTreeReaderArray<Float_t> hcal_reco_x(myReader,
                                        "HcalEndcapNClusters.position.x");
  TTreeReaderArray<Float_t> hcal_reco_y(myReader,
                                        "HcalEndcapNClusters.position.y");

  TTreeReaderArray<Float_t> emcal_reco_E(myReader,
                                         "EcalEndcapNClusters.energy");
  TTreeReaderArray<Float_t> emcal_reco_theta(
      myReader, "EcalEndcapNClusters.intrinsicTheta");
  TTreeReaderArray<Float_t> emcal_reco_phi(myReader,
                                           "EcalEndcapNClusters.intrinsicPhi");
  TTreeReaderArray<Float_t> emcal_reco_x(myReader,
                                         "EcalEndcapNClusters.position.x");
  TTreeReaderArray<Float_t> emcal_reco_y(myReader,
                                         "EcalEndcapNClusters.position.y");

  TTreeReaderArray<Float_t> ebarell_truth_E(myReader,
                                            "EcalBarrelTruthClusters.energy");
  TTreeReaderArray<Float_t> ebarell_truth_theta(
      myReader, "EcalBarrelTruthClusters.intrinsicTheta");
  TTreeReaderArray<Float_t> ebarell_truth_phi(
      myReader, "EcalBarrelTruthClusters.intrinsicPhi");
  TTreeReaderArray<Float_t> ebarell_truth_x(
      myReader, "EcalBarrelTruthClusters.position.x");
  TTreeReaderArray<Float_t> ebarell_truth_y(
      myReader, "EcalBarrelTruthClusters.position.y");

  TTreeReaderArray<Float_t> ebarell_reco_E(myReader,
                                           "EcalBarrelClusters.energy");
  TTreeReaderArray<Float_t> ebarell_reco_theta(
      myReader, "EcalBarrelClusters.intrinsicTheta");
  TTreeReaderArray<Float_t> ebarell_reco_phi(myReader,
                                             "EcalBarrelClusters.intrinsicPhi");
  TTreeReaderArray<Float_t> ebarell_reco_x(myReader,
                                           "EcalBarrelClusters.position.x");
  TTreeReaderArray<Float_t> ebarell_reco_y(myReader,
                                           "EcalBarrelClusters.position.y");

  TCanvas *can = new TCanvas("can", "can", 1200, 1000);
  can->SetMargin(0.09, 0.1, 0.1, 0.06);

  TFile *output = new TFile(outFileName, "recreate");
  TH2D *hPtEtaMc =
      new TH2D("hPtEtaMc", "hPtEtaMc;p_{t,gen} (GeV/c);#eta_{gen} ", 1500, 0.,
               10.0, 200, -5, 5);

  ClusterHists hcalRecoHist("HCal_Reco");
  ClusterHists emcalRecoHist("EMCal_Reco");
  ClusterHists ebarellRecHist("EBarell_Reco");
  ClusterHists hcalAndEmcalSumHist("HcalAndEcalSum_Reco");
  ClusterHists scatteredHcal("scattered_Hcal");

  ClusterHists hcalTruthHist("HCal_Truth");
  ClusterHists emcalTruthHist("EMCal_Truth");
  ClusterHists ebarellTruthHist("EBarell_Truth");
  ClusterHists hcalAndEmcalSumTruthHist("HcalAndEcalSum_Truth");
  ClusterHists scatteredHcalTruth("scattered_Hcal_Truth");
  Double_t maxEcalTheta = 0.;
  Double_t minEcalTheta = 1000.;
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////

  Int_t nEvents = myReader.GetEntries() / 1.;
  cout << "Total Events: " << nEvents << endl;

  for (Int_t iEvent = 0; iEvent < nEvents; ++iEvent) {
    myReader.SetEntry(iEvent);
    if (iEvent % 1000 == 0)
      cout << "Event No: " << iEvent << endl;
    // MC Particle
    Cluster mc("mc");
    for (int iParticle = 0; iParticle < charge.GetSize(); ++iParticle) {
      if (status[iParticle] == 1) {
        Double_t p_mc = sqrt(px_mc[iParticle] * px_mc[iParticle] +
                             py_mc[iParticle] * py_mc[iParticle] +
                             pz_mc[iParticle] * pz_mc[iParticle]);
        Double_t pt_mc = sqrt(px_mc[iParticle] * px_mc[iParticle] +
                              py_mc[iParticle] * py_mc[iParticle]);
        Double_t eta_mc =
            -1.0 *
            TMath::Log(TMath::Tan((TMath::ACos(pz_mc[iParticle] / p_mc)) / 2));
        Double_t theta_mc = TMath::ACos(pz_mc[iParticle] / p_mc);
        Double_t phi_mc = TMath::ATan(py_mc[iParticle] / px_mc[iParticle]);
        mc.theta = theta_mc;
        mc.phi = phi_mc;
        mc.energy = p_mc;
        mc.xVertex = vx_mc[iParticle];
        mc.yVertex = vy_mc[iParticle];

        mc.normalize();
        hPtEtaMc->Fill(pt_mc, eta_mc);
      }
    }

    Int_t indexHClusterTruth = getMaximumEnergyIndex(hcal_truth_E);
    Int_t indexHClusterReco = getMaximumEnergyIndex(hcal_reco_E);
    Int_t indexEClusterReco = getMaximumEnergyIndex(emcal_reco_E);
    Int_t indexEClusterTruth = getMaximumEnergyIndex(emcal_truth_E);
    Int_t indexEBarellClusterReco = getMaximumEnergyIndex(ebarell_reco_E);
    Int_t indexEBarellClusterTruth = getMaximumEnergyIndex(ebarell_truth_E);

    // if (indexHClusterReco < 0 && indexEClusterReco < 0)
    //   continue;

    Cluster hcalReco("hcalReco");
    Cluster emcalReco("emcalReco");
    Cluster ebarellReco("ebarellReco");

    Cluster hcalTruth("hcalTruth");
    Cluster emcalTruth("emcalTruth");
    Cluster ebarellTruth("ebarellTruth");

    if (indexHClusterReco >= 0) // only hcal
    {
      hcalReco.energy = hcal_reco_E[indexHClusterReco];
      hcalReco.theta = hcal_reco_theta[indexHClusterReco];
      hcalReco.phi = hcal_reco_phi[indexHClusterReco];

      hcalReco.xVertex = hcal_reco_x[indexHClusterReco];
      hcalReco.yVertex = hcal_reco_y[indexHClusterReco];

      hcalReco = hcalReco - mc;
      hcalReco.normalize();

      if (indexEClusterReco < 0 && indexEBarellClusterReco < 0)
        hcalRecoHist.Fill(hcalReco);
      else
        scatteredHcal.Fill(hcalReco);
    }

    if (indexEBarellClusterReco >= 0) {
      ebarellReco.energy = ebarell_reco_E[indexEBarellClusterReco];
      ebarellReco.theta = ebarell_reco_theta[indexEBarellClusterReco];
      ebarellReco.phi = ebarell_reco_phi[indexEBarellClusterReco];
      ebarellReco.xVertex = ebarell_reco_x[indexEBarellClusterReco];
      ebarellReco.yVertex = ebarell_reco_y[indexEBarellClusterReco];
      ebarellReco = ebarellReco - mc;
      ebarellReco.normalize();

      if (indexHClusterReco < 0)
        ebarellRecHist.Fill(ebarellReco);
    }

    if (indexEClusterReco >= 0) // only hcal
    {
      emcalReco.energy = emcal_reco_E[indexEClusterReco];
      emcalReco.theta = emcal_reco_theta[indexEClusterReco];
      emcalReco.phi = emcal_reco_phi[indexEClusterReco];
      emcalReco.xVertex = emcal_reco_x[indexEClusterReco];
      emcalReco.yVertex = emcal_reco_y[indexEClusterReco];
      emcalReco = emcalReco - mc;
      emcalReco.normalize();

      if (indexHClusterReco < 0)
        emcalRecoHist.Fill(emcalReco);

      if (emcalReco.theta > maxEcalTheta)
        maxEcalTheta = emcalReco.theta;

      if (emcalReco.theta < minEcalTheta)
        minEcalTheta = emcalReco.theta;
    }

    if (indexHClusterReco >= 0 && indexEClusterReco >= 0) {
      Cluster hcalAndEmcalSum = hcalReco;
      hcalAndEmcalSum.type = "hcalAndEmcalSum";
      hcalAndEmcalSum = hcalReco.addEcalWithSampleCoefficient(1.58, emcalReco);
      hcalAndEmcalSum = hcalAndEmcalSum - mc;
      hcalAndEmcalSumHist.Fill(hcalAndEmcalSum);
    }

    if (indexHClusterTruth >= 0) // only hcal
    {
      hcalTruth.energy = hcal_truth_E[indexHClusterTruth];
      hcalTruth.theta = hcal_truth_theta[indexHClusterTruth];
      hcalTruth.phi = hcal_truth_phi[indexHClusterTruth];
      hcalTruth.xVertex = hcal_truth_x[indexHClusterTruth];
      hcalTruth.yVertex = hcal_truth_y[indexHClusterTruth];
      hcalTruth = hcalTruth - mc;
      hcalTruth.normalize();

      if (indexEClusterTruth < 0 && indexEBarellClusterTruth < 0)
        hcalTruthHist.Fill(hcalTruth);
      else
        scatteredHcalTruth.Fill(hcalTruth);
    }

    if (indexEBarellClusterTruth >= 0) {
      ebarellTruth.energy = ebarell_truth_E[indexEBarellClusterTruth];
      ebarellTruth.theta = ebarell_truth_theta[indexEBarellClusterTruth];
      ebarellTruth.phi = ebarell_truth_phi[indexEBarellClusterTruth];
      ebarellTruth.xVertex = ebarell_truth_x[indexEBarellClusterTruth];
      ebarellTruth.yVertex = ebarell_truth_y[indexEBarellClusterTruth];
      ebarellTruth = ebarellTruth - mc;
      ebarellTruth.normalize();

      if (indexHClusterTruth < 0)
        ebarellTruthHist.Fill(ebarellTruth);
    }

    if (indexEClusterTruth >= 0) // only emcal
    {
      emcalTruth.energy = emcal_truth_E[indexEClusterTruth];
      emcalTruth.theta = emcal_truth_theta[indexEClusterTruth];
      emcalTruth.phi = emcal_truth_phi[indexEClusterTruth];
      emcalTruth.xVertex = emcal_truth_x[indexEClusterTruth];
      emcalTruth.yVertex = emcal_truth_y[indexEClusterTruth];
      emcalTruth = emcalTruth - mc;
      emcalTruth.normalize();

      if (indexHClusterTruth < 0)
        emcalTruthHist.Fill(emcalTruth);
    }

    if (indexHClusterTruth >= 0 && indexEClusterTruth >= 0) {
      Cluster hcalAndEmcalSumTruth = hcalTruth;
      hcalAndEmcalSumTruth.type = "hcalAndEmcalSumTruth";
      hcalAndEmcalSumTruth =
          hcalTruth.addEcalWithSampleCoefficient(1.58, emcalTruth);
      hcalAndEmcalSumTruth = hcalAndEmcalSumTruth - mc;
      hcalAndEmcalSumTruthHist.Fill(hcalAndEmcalSumTruth);
    }

  } // Event For loop

  outFileName.ReplaceAll(".root", "");
  TString outPdf = outFileName + "QaCheck.pdf";

  can->cd();
  can->SaveAs(outPdf + "[");
  hPtEtaMc->Draw("colz");
  can->SaveAs(outPdf);

  vector<Double_t> hcalRecoResolutions = hcalRecoHist.Draw(can, outPdf);
  vector<Double_t> ebarellRecoResolutions = ebarellRecHist.Draw(can, outPdf);
  vector<Double_t> emcalRecoResolutions = emcalRecoHist.Draw(can, outPdf);
  vector<Double_t> hcalAndEmcalSumRecoResolutions =
      hcalAndEmcalSumHist.Draw(can, outPdf);
  hcalRecoHist.DrawTogether(can, outPdf, emcalRecoHist, hcalAndEmcalSumHist,
                            scatteredHcal);

  vector<Double_t> hcalTruthResolutions = hcalTruthHist.Draw(can, outPdf);
  vector<Double_t> ebarellTruthResolutions = ebarellTruthHist.Draw(can, outPdf);
  vector<Double_t> emcalTruthResolutions = emcalTruthHist.Draw(can, outPdf);
  vector<Double_t> hcalAndEmcalSumTruthResolutions =
      hcalAndEmcalSumTruthHist.Draw(can, outPdf);
  hcalTruthHist.DrawTogether(can, outPdf, emcalTruthHist,
                             hcalAndEmcalSumTruthHist, scatteredHcalTruth);

  output->cd();

  string phiString = "";
  string thetaString = "";

  // string filename =
  // "eicrecon_neutron_50000events_p5gev_phi78_theta155.872.edm4eic.root"
  std::string filename = inFileName.Data();
  std::regex phiRegex("phi(\\d+\\.?\\d*)");
  std::regex thetaRegex("theta(\\d+\\.?\\d*)");
  std::smatch phiMatch, thetaMatch;

  if (std::regex_search(filename, phiMatch, phiRegex) &&
      std::regex_search(filename, thetaMatch, thetaRegex)) {
    phiString = phiMatch.str(1);
    thetaString = thetaMatch.str(1);
    cout << "Phi =" << phiString << endl;
    cout << "Theta =" << thetaString << endl;
  }

  vector<TString> labels = {"hcal+emcal #theta",      "hcal+emcal #phi",
                            "hcal+emcal #Rxy",        "hcal only #theta",
                            "hcal only #phi",         "hcal only #Rxy",
                            "emcal only #theta",      "emcal only #phi",
                            "emcal only #Rxy",        "hcal+emcal #theta Truth",
                            "hcal+emcal #phi Truth",  "hcal+emcal #Rxy Truth",
                            "hcal only #theta Truth", "hcal only #phi Truth",
                            "hcal only #Rxy Truth",   "emcal only #theta Truth",
                            "emcal only #phi Truth",  "emcal only #Rxy Truth"};

  vector<Double_t> labelsValues = {hcalAndEmcalSumRecoResolutions[0],
                                   hcalAndEmcalSumRecoResolutions[1],
                                   hcalAndEmcalSumRecoResolutions[2],
                                   hcalRecoResolutions[0],
                                   hcalRecoResolutions[1],
                                   hcalRecoResolutions[2],
                                   emcalRecoResolutions[0],
                                   emcalRecoResolutions[1],
                                   emcalRecoResolutions[2],
                                   hcalAndEmcalSumTruthResolutions[0],
                                   hcalAndEmcalSumTruthResolutions[1],
                                   hcalAndEmcalSumTruthResolutions[2],
                                   hcalTruthResolutions[0],
                                   hcalTruthResolutions[1],
                                   hcalTruthResolutions[2],
                                   emcalTruthResolutions[0],
                                   emcalTruthResolutions[1],
                                   emcalTruthResolutions[2]};

  TH1D *hResolution = new TH1D(
      Form("hResolutionPhi%sTheta%s", phiString.data(), thetaString.data()),
      "Sigma values", labels.size(), 0, labels.size());

  for (int i = 0; i < labels.size(); ++i) {
    hResolution->GetXaxis()->SetBinLabel(i + 1, labels[i]);
    hResolution->SetBinContent(i + 1, labelsValues[i]);
  }

  can->cd();
  hResolution->Draw("hist");
  can->SaveAs(outPdf);
  can->SaveAs(outPdf + "]");

  output->cd();
  hResolution->Write();
  TDirectory *dir =
      output->mkdir(Form("Phi%sTheta%s", phiString.data(), thetaString.data()));
  dir->cd();

  hcalRecoHist.Write(dir);
  emcalRecoHist.Write(dir);
  hcalAndEmcalSumHist.Write(dir);
  ebarellRecHist.Write(dir);
  scatteredHcal.Write(dir);

  hcalTruthHist.Write(dir);
  emcalTruthHist.Write(dir);
  hcalAndEmcalSumTruthHist.Write(dir);
  ebarellTruthHist.Write(dir);
  scatteredHcalTruth.Write(dir);

  output->Save();
  output->Close();
}
