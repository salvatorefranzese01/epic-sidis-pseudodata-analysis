//-------------- ePIC extraction table for uTMDs
//--- Based on the previous version of: Lorenzo Polizzi (lorenzo.polizzi@unife.it), Sara Pucillo (sara.pucillo@cern.ch)




#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "TText.h"
#include <TLatex.h>
#include <TPaveStats.h>
#include <TPaveStatsEditor.h>
#include "TPaletteAxis.h"
#include "TPolyLine.h"
#include "TStyle.h"
#include "TColor.h"
#include <cstring> 
#include <fstream>
#include <filesystem> 
#include <set>
#include <tuple>


using namespace std;
namespace fs = std::filesystem;



// Effective MC luminosity from the generated event count and sample cross section.
static double get_lumi_mc(double q_low)
{
    // total_cc in pb
    double total_cc, n_gen;
    if      (q_low < 3.0)    { total_cc = 0.70782586388685187e6; n_gen = 5e6; }  // Q2 1-10  
    else if (q_low < 10.0)   { total_cc = 5.9577669518852711e4;  n_gen = 5e6; }  // Q2 10-100
    else if (q_low < 31)  		{ total_cc = 2.1920402002225015e3;  n_gen = 5e6; }  // Q2 100-1000
    else                      { total_cc = 2.5939732135782342e1;  n_gen = 2e6; }  // Q2 1000-10000
    return n_gen / total_cc;
}

// ------- BINNING -------

// Binning:
// xB  : 1e-4 -> 1
// Q  : 1 -> 100 GeV
// z   : 0.1     -> 1
// Pt  : 0     -> 4 GeV/c
//
// For a given value, return the corresponding bin.
// Bins are numbered starting from 1.
// Returns -1 if the value is outside the range.

int getBin_xB(double xB)
{
    std::vector<double> xB_edges = {1.0e-4, 1.59e-4, 2.51e-4, 3.98e-4, 6.31e-4,
        1.0e-3, 1.59e-3, 2.51e-3, 3.98e-3, 6.31e-3,
        1.0e-2, 1.59e-2, 2.51e-2, 3.98e-2, 6.31e-2,
        1.0e-1, 1.59e-1, 2.51e-1, 3.98e-1, 6.31e-1, 1.0};

    for(size_t i = 0; i < xB_edges.size() - 1; i++)
    {
        if(xB >= xB_edges[i] && xB < xB_edges[i+1])
        {
            return i + 1;
        }
    }
    return -1;  // Outside the bin range.
}


int getBin_Q(double Q)
{
    vector<double> Q_edges = {
    1.0,
    1.33417,
    1.77764,
    2.37065,
    3.16228,
    4.21900,
    5.62139,
    7.49667,
    10.0,
    13.3417,
    17.7764,
    23.7065,
    31.6228,
    100.0
};
    
    for(size_t i = 0; i < Q_edges.size() - 1; i++)
    {
    	if(Q >= Q_edges[i] && Q < Q_edges[i+1])
    	{
      	return i + 1;  
      }
    }
    return -1;  // Outside the bin range.
}


int getBin_z(double z)
{
    std::vector<double> z_edges = {0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    for(size_t i = 0; i < z_edges.size() - 1; i++)
    {
        if(z >= z_edges[i] && z < z_edges[i+1])
        {
            return i + 1;
        }
    }
    return -1;  // Outside the bin range.
}


int getBin_Pt(double Pt)
{
    std::vector<double> Pt_edges = {0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0};

    for(size_t i = 0; i < Pt_edges.size() - 1; i++)
    {
        if(Pt >= Pt_edges[i] && Pt < Pt_edges[i+1])
        {
            return i + 1;
        }
    }
    return -1;  // Outside the bin range.
}



// ------- BIN RANGES -------

std::pair<double,double> getBinRange_xB(int bin)
{
    std::vector<double> xB_edges = {1.0e-4, 1.59e-4, 2.51e-4, 3.98e-4, 6.31e-4,
        1.0e-3, 1.59e-3, 2.51e-3, 3.98e-3, 6.31e-3,
        1.0e-2, 1.59e-2, 2.51e-2, 3.98e-2, 6.31e-2,
        1.0e-1, 1.59e-1, 2.51e-1, 3.98e-1, 6.31e-1, 1.0};

    if(bin < 1 || bin >= xB_edges.size())
        return {-1, -1};

    return {xB_edges[bin-1], xB_edges[bin]};
}


// Q binning
// Q edges corresponding to the requested Q2 binning
std::pair<double,double> getBinRange_Q(int bin)
{
    vector<double> Q_edges = {
    1.0,
    1.33417,
    1.77764,
    2.37065,
    3.16228,
    4.21900,
    5.62139,
    7.49667,
    10.0,
    13.3417,
    17.7764,
    23.7065,
    31.6228,
    100.0
};

    if(bin < 1 || bin >= Q_edges.size())
        return {-1, -1};

    return {Q_edges[bin-1], Q_edges[bin]};
}


std::pair<double,double> getBinRange_z(int bin)
{
    std::vector<double> z_edges = {0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    if(bin < 1 || bin >= z_edges.size())
        return {-1, -1};

    return {z_edges[bin-1], z_edges[bin]};
}


std::pair<double,double> getBinRange_Pt(int bin)
{
    std::vector<double> Pt_edges = {0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0};

    if(bin < 1 || bin >= Pt_edges.size())
        return {-1, -1};

    return {Pt_edges[bin-1], Pt_edges[bin]};
}

// ------- Fine Binning ---------------



//------------- Main function
void epic_extraction_cleaned(int target_pdg = 211, const char* inputDir = "/input/10to100", int beam_e = 10, int beam_p = 250) 
{
	
	// --- Hadron selection (PDG code)
  TString tag, label;
  switch (target_pdg) 
  {
  	case  211: tag = "pipos"; label = "#pi^{+}"; break;
    case -211: tag = "pineg"; label = "#pi^{-}"; break;
    case  321: tag = "kaonpos";  label = "K^{+}"; break;
    case -321: tag = "kaonneg";  label = "K^{-}"; break;
    
    default:
            tag   = Form("pdg%d", target_pdg);
            label = Form("PDG %d", target_pdg);
            break;
  }
  
  cout << tag << " has been selected" << endl;
  
		// Variables
    // --- electron
    double electron_px, electron_py, electron_pz, electron_mom, electron_Theta, electron_Phi, electron_ThetaDeg, electron_E, electron_W, electron_Q2, electron_ass;
    double electron_eta, electron_y;
    // -- hadron
    double hadron_mom, hadron_Q2, hadron_xB, hadron_xF, hadron_z, hadron_PhT, hadron_Phi_h, hadron_Phi_s, hadron_Phi_lab, hadron_Theta, hadron_eta, hadron_y, hadron_W2, hadron_Mx;
    double helicity, eps, hadron_px, hadron_py, hadron_pz, el_px, el_py, el_pz, el_theta, el_phi, el_eta, el_mom, pr_mom, pr_px, pr_py, pr_pz, pr_phi, pr_theta, pr_eta;
    double rec_pdg, good_PID, el_rec_pdg, rec_pdg_mc;
    double el_ass_rec_pdg, el_ass_px, el_ass_py, el_ass_pz, el_ass_theta, el_ass_phi, el_ass_eta, el_ass_mom;
    double hadron_E;
    //
    double hadron_mom_mc, hadron_Q2_mc, hadron_xB_mc, hadron_xF_mc, hadron_z_mc, hadron_PhT_mc;
    double hadron_Phi_h_mc, hadron_Phi_s_mc, hadron_Phi_lab_mc, hadron_Theta_mc, hadron_eta_mc, hadron_y_mc, hadron_W2_mc, hadron_Mx_mc;
    double hel_mc, eps_mc, hadron_px_mc, hadron_py_mc, hadron_pz_mc;
    double hadron_mc_index, index_mc;
    //
    double hadron_mom_all, hadron_Q2_all, hadron_xB_all, hadron_xF_all, hadron_z_all, hadron_PhT_all;
    double hadron_Phi_h_all, hadron_Phi_s_all, hadron_Phi_lab_all, hadron_Theta_all, hadron_eta_all, hadron_y_all, hadron_W2_all, hadron_Mx_all;
    double hel_all, eps_all, hadron_px_all, hadron_py_all, hadron_pz_all;
    double good_PID_all, pdg_all;
    double hadron_index, index_all;
    
    // new variables for the efficiencies (compute on the matched subset)
    double hadron_Q2_mc_matched, hadron_xB_mc_matched, hadron_z_mc_matched, hadron_PhT_mc_matched;
    double hadron_y_mc_matched; 

    //--- input file
    string inputDirStr = inputDir;
    TTree treeHadron(Form("%s_RECO", tag.Data()), Form("RECO %s", label.Data()));
    TTree treeHadron_MC(Form("%s_MC", tag.Data()), Form("MC %s", label.Data()));
    TChain chainElectron("Electron");
    TChain chainHadron_Reco("Hadron Reco");
    TChain chainHadron_MC("Hadron MC");

    int fileCount = 0;

    fs::path inPath(inputDirStr);
    if (fs::exists(inPath) && fs::is_regular_file(inPath) && inPath.extension() == ".root") {
        string filePath = inPath.string();
        chainElectron.Add(Form("%s/ElectronTree_MC", filePath.c_str()));
        chainHadron_MC.Add(Form("%s/HadronTree_MC", filePath.c_str()));
        chainHadron_Reco.Add(Form("%s/HadronTree_RECO", filePath.c_str()));
        fileCount = 1;
    } else {
        for (const auto &entry : fs::directory_iterator(inputDirStr)) {
            if (entry.path().extension() == ".root") {
                string filePath = entry.path().string();
                chainElectron.Add(Form("%s/ElectronTree_MC", filePath.c_str()));
                chainHadron_MC.Add(Form("%s/HadronTree_MC", filePath.c_str()));
                chainHadron_Reco.Add(Form("%s/HadronTree_RECO", filePath.c_str()));
                fileCount++;
            }
        }
    }

    if (fileCount == 0) {
        cerr << "No .root file found in " << inputDirStr << endl;
        return;
    }

    //--- output file
    TString outputFile = Form("epic_2604_%s_%dx%d.root", tag.Data(), beam_e, beam_p);
    TFile outFile(outputFile, "RECREATE");

    // yaml file

    TString yaml_name = Form("%s_%dx%d_26_07.yaml", tag.Data(), beam_e, beam_p);
    std::ofstream yamlFile(yaml_name);

    // csv file

    TString csv_filename = Form("%s_%dx%d_26_07.csv", tag.Data(), beam_e, beam_p);
    std::ofstream csvFile(csv_filename);
    //csvFile << "id, s[GeV^{2}], <Q>[GeV], Q_min[GeV], Q_max[GeV], <xB>, xB_min, xB_max, <z>, z_min, z_max, <pT>, pT_min, pT_max, Delta_bin[GeV^3], xSec_diff[pb/GeV^3], xSec_uncorr_error, xSec_corr_error, y_min, y_max, W2_min[GeV^2], W2_max[GeV^2],target_mass[GeV], product_mass[GeV], N_gen, N_reco, abs_stat_error, abs_sys_error, efficiency\n";  
    
    // with the n_dis counter
    csvFile << "id, s[GeV^{2}], <Q>[GeV], Q_min[GeV], Q_max[GeV], <xB>, xB_min, xB_max, <z>, z_min, z_max, <pT>, pT_min, pT_max, Delta_bin_[GeV^3], xSec_diff[pb/GeV^3], xSec_uncorr_error, xSec_corr_error, y_min, y_max, W2_min[GeV^2], W2_max[GeV^2],target_mass[GeV], product_mass[GeV], N_gen, N_reco, N_reco_matched, N_reco_rel_stat_error, N_reco_abs_sys_error, efficiency, N_DIS, N_DIS_rel_stat_error, N_DIS_abs_sys_error\n";

		//--- Here we collect all the variables from the ttree
    if (chainElectron.GetNtrees() > 0) {
        chainElectron.SetBranchAddress("el_px_mc", &electron_px);
        chainElectron.SetBranchAddress("el_py_mc", &electron_py);
        chainElectron.SetBranchAddress("el_pz_mc", &electron_pz);
        chainElectron.SetBranchAddress("el_mom_mc", &electron_mom);
        chainElectron.SetBranchAddress("el_theta_mc", &electron_Theta);
        chainElectron.SetBranchAddress("el_phi_mc", &electron_Phi);
        chainElectron.SetBranchAddress("el_eta_mc", &electron_eta);
        chainElectron.SetBranchAddress("el_y_mc", &electron_y);
    }

    //Hadron reco
    chainHadron_Reco.SetBranchAddress("hadron_index", &hadron_index);
    chainHadron_Reco.SetBranchAddress("hadron_pdg", &rec_pdg);
    chainHadron_Reco.SetBranchAddress("hadron_pdg_mc", &rec_pdg_mc);
    chainHadron_Reco.SetBranchAddress("hadron_good_PID", &good_PID);
    chainHadron_Reco.SetBranchAddress("hadron_px", &hadron_px);
    chainHadron_Reco.SetBranchAddress("hadron_py", &hadron_py);
    chainHadron_Reco.SetBranchAddress("hadron_pz", &hadron_pz);
    chainHadron_Reco.SetBranchAddress("hadron_mom", &hadron_mom);
    //chainHadron_Reco.SetBranchAddress("hadron_W2", &hadron_W2);
    chainHadron_Reco.SetBranchAddress("hadron_Q2", &hadron_Q2);
    //chainHadron_Reco.SetBranchAddress("hadron_xF", &hadron_xF);
    chainHadron_Reco.SetBranchAddress("hadron_xB", &hadron_xB);
    chainHadron_Reco.SetBranchAddress("hadron_y", &hadron_y);
    chainHadron_Reco.SetBranchAddress("hadron_z", &hadron_z);
    chainHadron_Reco.SetBranchAddress("hadron_PhT", &hadron_PhT);
    chainHadron_Reco.SetBranchAddress("hadron_Phi_lab", &hadron_Phi_lab);
    chainHadron_Reco.SetBranchAddress("hadron_Theta", &hadron_Theta);
    chainHadron_Reco.SetBranchAddress("hadron_eta", &hadron_eta);
    chainHadron_Reco.SetBranchAddress("hadron_Phi_h", &hadron_Phi_h);
    //chainHadron_Reco.SetBranchAddress("Phi_s", &hadron_Phi_s);
    //chainHadron_Reco.SetBranchAddress("helicity", &helicity);
    //chainHadron_Reco.SetBranchAddress("hadron_Mx", &hadron_Mx);
    
    // for the efficiencies
    chainHadron_Reco.SetBranchAddress("hadron_Q2_mc",  &hadron_Q2_mc_matched);
	chainHadron_Reco.SetBranchAddress("hadron_xB_mc",  &hadron_xB_mc_matched);
	chainHadron_Reco.SetBranchAddress("hadron_z_mc",   &hadron_z_mc_matched);
	chainHadron_Reco.SetBranchAddress("hadron_PhT_mc", &hadron_PhT_mc_matched);
	chainHadron_Reco.SetBranchAddress("hadron_y_mc",   &hadron_y_mc_matched);

    //Hadron MC
    chainHadron_MC.SetBranchAddress("hadron_index_mc", &index_mc);
    // Optional: if present, use MC PDG to select the hadron species
    double mc_pdg = 0;
    //bool has_mc_pdg = (chainHadron_MC.SetBranchAddress("hadron_pdg_mc", &mc_pdg) == 0);
    //if (!has_mc_pdg) has_mc_pdg = (chainHadron_MC.SetBranchAddress("hadron_pdg", &mc_pdg) == 0);
    chainHadron_MC.SetBranchAddress("hadron_pdg_mc", &mc_pdg);
    chainHadron_MC.SetBranchAddress("hadron_mom_mc", &hadron_mom_mc);
    chainHadron_MC.SetBranchAddress("hadron_Q2_mc", &hadron_Q2_mc);
    chainHadron_MC.SetBranchAddress("hadron_xB_mc", &hadron_xB_mc);
    //chainHadron_MC.SetBranchAddress("hadron_xF_mc", &hadron_xF_mc);
    chainHadron_MC.SetBranchAddress("hadron_z_mc", &hadron_z_mc);
    chainHadron_MC.SetBranchAddress("hadron_PhT_mc", &hadron_PhT_mc);
    chainHadron_MC.SetBranchAddress("hadron_Phi_lab_mc", &hadron_Phi_lab_mc);
    chainHadron_MC.SetBranchAddress("hadron_Phi_h_mc", &hadron_Phi_h_mc);
    //chainHadron_MC.SetBranchAddress("Phi_s_mc", &hadron_Phi_s_mc);
    chainHadron_MC.SetBranchAddress("hadron_Theta_mc", &hadron_Theta_mc);
    chainHadron_MC.SetBranchAddress("hadron_eta_mc", &hadron_eta_mc);
    chainHadron_MC.SetBranchAddress("hadron_y_mc", &hadron_y_mc);
    //chainHadron_MC.SetBranchAddress("hadron_W2_mc", &hadron_W2_mc);
    //chainHadron_MC.SetBranchAddress("hadron_Mx_mc", &hadron_Mx_mc);
    //chainHadron_MC.SetBranchAddress("helicity_mc", &hel_mc);
    //chainHadron_MC.SetBranchAddress("hadron_epsilon_mc", &eps_mc);
    chainHadron_MC.SetBranchAddress("hadron_px_mc", &hadron_px_mc);
    chainHadron_MC.SetBranchAddress("hadron_py_mc", &hadron_py_mc);
    chainHadron_MC.SetBranchAddress("hadron_pz_mc", &hadron_pz_mc);

    // Interesting hadron --- reco info
    treeHadron.Branch("mc_index", &hadron_index, "mc_index/I");
    treeHadron.Branch("rec_pdg", &rec_pdg, "rec_pdg/D");
    treeHadron.Branch("good_PID", &good_PID, "good_PID/D");
    treeHadron.Branch("px", &hadron_px, "hadron_px/D");
    treeHadron.Branch("py", &hadron_py, "hadron_py/D");
    treeHadron.Branch("pz", &hadron_pz, "hadron_pz/D");
    treeHadron.Branch("E", &hadron_E, "E/D");
    treeHadron.Branch("Mom", &hadron_mom, "Mom/D");
    treeHadron.Branch("Q2", &hadron_Q2, "Q2/D");
    treeHadron.Branch("xB", &hadron_xB, "xB/D");
    treeHadron.Branch("xF", &hadron_xF, "xF/D");
    treeHadron.Branch("z", &hadron_z, "z/D");
    treeHadron.Branch("PhT", &hadron_PhT, "PhT/D");
    treeHadron.Branch("Phi_lab", &hadron_Phi_lab, "Phi_h/D");
    treeHadron.Branch("Phi_h", &hadron_Phi_h, "Phi_h/D");
    treeHadron.Branch("Phi_s", &hadron_Phi_s, "Phi_s/D");
    treeHadron.Branch("theta", &hadron_Theta, "theta/D");
    treeHadron.Branch("eta", &hadron_eta, "eta/D");
    treeHadron.Branch("y", &hadron_y, "y/D");
    treeHadron.Branch("W2", &hadron_W2, "W/D");
    treeHadron.Branch("Mx", &hadron_Mx, "Mx/D");
    treeHadron.Branch("helicity", &helicity, "hel/D");
    treeHadron.Branch("epsilon", &eps, "eps/D");

    // Interesting hadron --- MC info
    treeHadron_MC.Branch("index", &index_mc, "index/I");
    treeHadron_MC.Branch("Mom_mc", &hadron_mom_mc, "Mom_mc/D");
    treeHadron_MC.Branch("Q2_mc", &hadron_Q2_mc, "Q2_mc/D");
    treeHadron_MC.Branch("xB_mc", &hadron_xB_mc, "xB_mc/D");
    treeHadron_MC.Branch("xF_mc", &hadron_xF_mc, "xF_mc/D");
    treeHadron_MC.Branch("z_mc", &hadron_z_mc, "z_mc/D");
    treeHadron_MC.Branch("PhT_mc", &hadron_PhT_mc, "PhT_mc/D");
    treeHadron_MC.Branch("Phi_lab_mc", &hadron_Phi_lab_mc, "Phi_lab_mc/D");
    treeHadron_MC.Branch("Phi_h_mc", &hadron_Phi_h_mc, "Phi_h_mc/D");
    treeHadron_MC.Branch("Phi_s_mc", &hadron_Phi_s_mc, "Phi_s_mc/D");
    treeHadron_MC.Branch("theta_mc", &hadron_Theta_mc, "theta_mc/D");
    treeHadron_MC.Branch("eta_mc", &hadron_eta_mc, "eta_mc/D");
    treeHadron_MC.Branch("y_mc", &hadron_y_mc, "y_mc/D");
    treeHadron_MC.Branch("W_mc", &hadron_W2_mc, "W_mc/D");
    treeHadron_MC.Branch("Mx_mc", &hadron_Mx_mc, "Mx_mc/D");
    treeHadron_MC.Branch("helicity_mc", &hel_mc, "helicity_mc/D");
    treeHadron_MC.Branch("epsilon_mc", &eps_mc, "epsilon_mc/D");
    treeHadron_MC.Branch("pion_px_mc", &hadron_px_mc, "pion_px_mc/D");
    treeHadron_MC.Branch("pion_py_mc", &hadron_py_mc, "pion_py_mc/D");
    treeHadron_MC.Branch("pion_pz_mc", &hadron_pz_mc, "pion_pz_mc/D");

    /*
    treeHadron_all.Branch("index", &index_all, "index/I");
    treeHadron_all.Branch("Mom_all", &hadron_mom_all, "Mom_all/D");
    treeHadron_all.Branch("Q2_all", &hadron_Q2_all, "Q2_all/D");
    treeHadron_all.Branch("xB_all", &hadron_xB_all, "xB_all/D");
    treeHadron_all.Branch("xF_all", &hadron_xF_all, "xF_all/D");
    treeHadron_all.Branch("z_all", &hadron_z_all, "z_all/D");
    treeHadron_all.Branch("PhT_all", &hadron_PhT_all, "PhT_all/D");
    treeHadron_all.Branch("Phi_lab_all", &hadron_Phi_lab_all, "Phi_lab_all/D");
    treeHadron_all.Branch("Phi_h_all", &hadron_Phi_h_all, "Phi_h_all/D");
    treeHadron_all.Branch("Phi_s_all", &hadron_Phi_s_all, "Phi_s_all/D");
    treeHadron_all.Branch("theta_all", &hadron_Theta_all, "theta_all/D");
    treeHadron_all.Branch("eta_all", &hadron_eta_all, "eta_all/D");
    treeHadron_all.Branch("y_all", &hadron_y_all, "y_all/D");
    treeHadron_all.Branch("W_all", &hadron_W2_all, "W_all/D");
    treeHadron_all.Branch("Mx_all", &hadron_Mx_all, "Mx_all/D");
    treeHadron_all.Branch("helicity_all", &hel_all, "helicity_all/D");
    treeHadron_all.Branch("epsilon_all", &eps_all, "epsilon_all/D");
    treeHadron_all.Branch("pion_px_all", &hadron_px_all, "pion_px_all/D");
    treeHadron_all.Branch("pion_py_all", &hadron_py_all, "pion_py_all/D");
    treeHadron_all.Branch("pion_pz_all", &hadron_pz_all, "pion_pz_all/D");
    */


		// 4D keys
  	struct EventData 
  	{
        double Q;
        double xB;
        double z;
        double Pt;
    };
    
    std::map<std::tuple<int,int,int,int>, std::vector<EventData>> data_multi;
    std::map<std::tuple<int,int,int,int>, std::vector<EventData>> data_multi_MC;    
    std::map<std::tuple<int,int>, int> data_DIS_reco;	// for the dis events
    std::map<std::tuple<int,int,int,int>, int> data_matched_truth_count;	// map for the efficiencies (matched events)

    vector<double> total_gen_xB;	
// Count generated and reconstructed hadrons by four-dimensional bin.
    // Hadron of interest
    Long64_t nEntries_had = chainHadron_Reco.GetEntries();
    Long64_t nEntries_hadMC = chainHadron_MC.GetEntries();
    
		// Loop over generated hadrons.
		// Apply truth-level kinematic and species selections before filling each bin.
    for (Long64_t i = 0; i < nEntries_hadMC; i++) 
    {
        chainHadron_MC.GetEntry(i);
        //if(i < 5) cout << "mc_pdg: " << mc_pdg << " target: " << target_pdg << endl;
        if (i % 100000 == 0) cout << "MC entry: " << i << "/" << nEntries_hadMC << endl;
        if(hadron_y_mc <= 0.95 && hadron_y_mc >= 0.01 && hadron_Q2_mc >= 1 && hadron_z_mc < 1)
        {
        		//if(i < 5) cout << "passa filtro cinematico, mc_pdg: " << mc_pdg << endl;
            if (mc_pdg != target_pdg) continue;
            treeHadron_MC.Fill();
            int index_Q_mc = getBin_Q(sqrt(hadron_Q2_mc));
            int index_xB_mc = getBin_xB(hadron_xB_mc);
            int index_z_mc = getBin_z(hadron_z_mc);
            int index_Pt_mc = getBin_Pt(hadron_PhT_mc);
            if(index_Q_mc > 0 && index_xB_mc > 0 && index_z_mc > 0 && index_Pt_mc > 0)
            {
            	//vec_multi_Q2[index_Q2-1][index_xB-1][index_z-1][index_Pt-1].push_back(sqrt(hadron_Q2));
              total_gen_xB.push_back(hadron_xB_mc);
              auto key_mc = std::make_tuple(index_Q_mc, index_xB_mc, index_z_mc, index_Pt_mc);
              EventData ev_mc;
              ev_mc.Q  = sqrt(hadron_Q2_mc);
              ev_mc.xB = hadron_xB_mc;
              ev_mc.z  = hadron_z_mc;
              ev_mc.Pt = hadron_PhT_mc;

              data_multi_MC[key_mc].push_back(ev_mc);
            }
        }
    }
    
    // Loop over reconstructed hadrons.
		// Apply reconstructed selections and fill four-dimensional bins.
    for (Long64_t i = 0; i < nEntries_had; i++) 
    {
        chainHadron_Reco.GetEntry(i);
        hadron_W2 = 0.938272*0.938272 + hadron_Q2*(1.0/hadron_xB - 1.0);  // Compute W squared from reconstructed DIS kinematics.
        if (i % 100000 == 0) cout << "RECO entry: " << i << "/" << nEntries_had << endl;
            
            // Require the reconstructed PID and truth species to match the target hadron.
            if(hadron_y <= 0.95 && hadron_y >= 0.01 && good_PID != -1 && rec_pdg == target_pdg && rec_pdg_mc == target_pdg && hadron_z < 1 && hadron_Q2 >= 1)
            { // Apply reconstructed-level kinematic and particle selections.
            
            	// Count matched reconstructed hadrons in the corresponding generated bin.
				// The generated kinematics must pass the same truth-level selection as N_gen.
				if (hadron_y_mc_matched <= 0.95 && hadron_y_mc_matched >= 0.01 && hadron_Q2_mc_matched >= 1 && hadron_z_mc_matched < 1)
				{
    				int index_Q_tm  = getBin_Q(sqrt(hadron_Q2_mc_matched));
    				int index_xB_tm = getBin_xB(hadron_xB_mc_matched);
    				int index_z_tm  = getBin_z(hadron_z_mc_matched);
    				int index_Pt_tm = getBin_Pt(hadron_PhT_mc_matched);
    				if (index_Q_tm > 0 && index_xB_tm > 0 && index_z_tm > 0 && index_Pt_tm > 0)
    				{
        				auto key_tm = std::make_tuple(index_Q_tm, index_xB_tm, index_z_tm, index_Pt_tm);
        				data_matched_truth_count[key_tm]++;
    				}
				}
               
              treeHadron.Fill();
                
              int index_Q = getBin_Q(sqrt(hadron_Q2));
              int index_xB = getBin_xB(hadron_xB);
              int index_z = getBin_z(hadron_z);
              int index_Pt = getBin_Pt(hadron_PhT);
              
              if(index_Q > 0 && index_xB > 0 && index_z > 0 && index_Pt > 0)
              {
              	//vec_multi_Q2[index_Q2-1][index_xB-1][index_z-1][index_Pt-1].push_back(sqrt(hadron_Q2));
                auto key = std::make_tuple(index_Q, index_xB, index_z, index_Pt);
                EventData ev;
                ev.Q  = sqrt(hadron_Q2);
                ev.xB = hadron_xB;
                ev.z  = hadron_z;
                ev.Pt = hadron_PhT;

                data_multi[key].push_back(ev);
               }
            }
            
		}
		
		treeHadron.Write();
    treeHadron_MC.Write();	
    
    // Collect all populated generated and reconstructed bin keys.
    std::set<std::tuple<int,int,int,int>> all_keys;
    
    for (const auto& [key, _] : data_multi)
    		all_keys.insert(key);

		for (const auto& [key, _] : data_multi_MC)
    		all_keys.insert(key);
    		
    std::cout << "all_keys size before electron loop = "
     << all_keys.size() << endl;
     
     // Count selected DIS electrons in each (Q, xB) bin.
    Long64_t nEntries_el = chainElectron.GetEntries();
    
		for (Long64_t i = 0; i < nEntries_el; i++) 
		{
			
    	chainElectron.GetEntry(i);
    	// Calculate Q squared and xB from the electron.
			double el_E = sqrt(electron_px*electron_px + electron_py*electron_py + electron_pz*electron_pz);
			double el_Q2_calc = 2.0 * beam_e * el_E * (1.0 + electron_pz/el_E);	// Electron beam travels along the negative z direction.
			double el_xB_calc = el_Q2_calc / (4.0 * beam_e * beam_p * electron_y);

    	// Apply the DIS kinematic selection.
    	if (electron_y >= 0.01 && electron_y <= 0.95 && el_Q2_calc >= 1.0) 
    	{
    		// Find the corresponding (Q, xB) bin.
        int iq = getBin_Q(sqrt(el_Q2_calc));
        int ix = getBin_xB(el_xB_calc);
        if (iq > 0 && ix > 0)
				{
    			auto k = std::make_tuple(iq, ix);
    			data_DIS_reco[k]++;  // No z or transverse-momentum binning for DIS counts.
				}
    	}
		}
     
     
     int n_id = 1;
     
 
		//double luminosity_mc   = 5e5 / total_cc;	// 1000 to 10000 has only 500k events instead of 5M
	    
    // -------- Differential cross section --------
    double N_DIS_stat_error = 0.0;
		double N_DIS_sys_error = 0.0;
		double N_rec_stat_error = 0.0;
		double N_rec_sys_error = 0.0;
    
    //loop sulle key
    for (const auto& key : all_keys) 
    {

        auto [q, x, z, pt] = key;
        

        // -------- RECO --------
        std::vector<EventData> values;
        if (data_multi.count(key))
            values = data_multi[key];

        // -------- MC --------
        // N_gen counts selected generated hadrons in this four-dimensional bin.
        int N_gen = data_multi_MC.count(key) ? data_multi_MC[key].size() : 0;

				// N_rec counts selected reconstructed hadrons in this bin.
        int N_rec = values.size();
        if(N_rec == 0) continue;	// Skip bins with no reconstructed hadrons.
        
        //----- matched (eff) -----
        int N_rec_matched = data_matched_truth_count.count(key) ? data_matched_truth_count[key] : 0;

        // -------- medie (RECO) --------
        // Calculate mean reconstructed kinematics within the bin.
        double sum_Q = 0.0, sum_xB = 0.0, sum_z = 0.0, sum_Pt = 0.0;

        for (const auto& ev : values) {
            sum_Q  += ev.Q;
            sum_xB += ev.xB;
            sum_z  += ev.z;
            sum_Pt += ev.Pt;
        }

        double mean_Q  = (N_rec > 0) ? sum_Q  / N_rec : -1;
        double mean_xB = (N_rec > 0) ? sum_xB / N_rec : -1;
        double mean_z  = (N_rec > 0) ? sum_z  / N_rec : -1;
        double mean_Pt = (N_rec > 0) ? sum_Pt / N_rec : -1;
        

        // -------- bin edges --------
        auto [Qmin, Qmax]   = getBinRange_Q(q);
        auto [xBmin, xBmax] = getBinRange_xB(x);
        auto [zmin, zmax]   = getBinRange_z(z);
        auto [ptmin, ptmax] = getBinRange_Pt(pt);

        //double delta_wrong = (Qmax-Qmin)*(xBmax-xBmin)*(zmax-zmin)*(ptmax-ptmin);
        double delta_Q2 = (Qmax * Qmax) - (Qmin * Qmin);;
        double delta = delta_Q2 * (xBmax-xBmin) * (zmax-zmin) * (ptmax-ptmin);
        
		    double luminosity_mc = get_lumi_mc(Qmin);
		    
		    
		    
		    double luminosity_proj = 2500.0;	// 2.5 fb^-1 in pb^-1
				double scale = sqrt(luminosity_mc / luminosity_proj);
				
				if(n_id % 400 == 0)
					cout << " Luminosity MC " << luminosity_mc << " pb^-1" << endl;


        // -------- cross section --------
        // Calculate the differential cross section and its uncertainties.
        double sigma = 0.0;
        double sigma_stat = 0.0;
        double sigma_sys = 0.0;
        double sigma_stat_proj = 0.0;
        double eff = 0.0;
        
        double epsilon = sqrt(0.035*0.035 + 0.015*0.015); // Combined 3.5% and 1.5% systematic contributions.
        
        // DIS counts for multiplicity plots.
        auto key_dis = std::make_tuple(q, x);
				int N_DIS = data_DIS_reco.count(key_dis) ? data_DIS_reco[key_dis] : 0;
				N_DIS_stat_error = (N_DIS > 0) ? (1.0 / sqrt((double)N_DIS)) * scale : 0.0;	// Scale the relative statistical uncertainty to the target luminosity.
				
				N_DIS_sys_error  = N_DIS * epsilon;  // Apply the same systematic fraction to the DIS count.
				
				// Reconstructed-hadron count uncertainties.
				N_rec_stat_error = (N_rec > 0) ? (1.0/sqrt(N_rec)) * scale : 0.0;
				N_rec_sys_error = N_rec * epsilon; 
				    
				
        if (delta > 0 && N_rec > 0 && N_gen > 0) 
        {
            
            // Bin efficiency from truth-matched reconstructed and selected generated hadrons.
            eff = (N_gen > 0) ? (double)N_rec_matched / (double)N_gen : 0.0;
            
            // Efficiency-corrected differential cross section.
            sigma = (double)N_rec / (luminosity_mc * eff * delta);
            

            // Statistical uncertainty from the reconstructed yield and bin efficiency.
            double sigma2 = sigma*sigma;
            double delta_eff = sqrt((eff*(1-eff))/(double)N_gen);
            sigma_stat = sigma*sqrt((1/(double)N_rec) + (delta_eff*delta_eff)/(eff*eff)); 
            
            // Project the statistical uncertainty to the target luminosity.
            double luminosity_proj = 2500.0; // 2.5 fb^-1 = 2500 pb^-1
            //double scale = sqrt(luminosity_mc / luminosity_proj);
            sigma_stat_proj = sigma_stat * scale;

            // Apply the assigned systematic uncertainty.
            
            sigma_sys = sigma * epsilon;
        } 

        double prod_mass = 0;
        if (std::abs(target_pdg) == 321) prod_mass = 0.497;
        else if (std::abs(target_pdg) == 211) prod_mass = 0.139;
       
       // Write YAML and CSV tables.


		// YAML
		//yamlFile << "  - bin: [" << q2 << ", " << x << ", " << z << ", " << pt << "]\n";
		yamlFile << "  - bin:\n";
		yamlFile << "      id: " << n_id << "\n";
		yamlFile << "      s_GeV2: " << 4.0 * beam_e * beam_p << "\n";
		yamlFile << "      Q:\n";
		yamlFile << "        Q_mean: " << mean_Q << "\n";
		yamlFile << "        Q_min: " << Qmin << "\n";
		yamlFile << "        Q_max: " << Qmax << "\n";
		yamlFile << "      xB:\n";
		yamlFile << "        xB_mean: "<< mean_xB <<"\n";
		yamlFile << "        xB_min: " << xBmin << "\n";
		yamlFile << "        xB_max: " << xBmax << "\n";
		yamlFile << "      z:\n";
		yamlFile << "        z_mean: " << mean_z << "\n";
		yamlFile << "        z_min: " << zmin << "\n";
		yamlFile << "        z_max: " << zmax << "\n";
		yamlFile << "      Pt:\n";
		yamlFile << "        Pt_mean: " << mean_Pt << "\n";
		yamlFile << "        Pt_min: " << ptmin << "\n";
		yamlFile << "        Pt_max: " << ptmax << "\n";
		yamlFile << "      Delta_bin: "<< delta <<"\n";
		yamlFile << "      xSec_diff: "<< sigma <<"\n";
		yamlFile << "      xSec_uncorr_error: " << sigma_stat_proj << "\n";
		yamlFile << "      xSec_corr_error: " << sigma_sys << "\n";
		yamlFile << "      y_min: " << 0.01 << "\n";
		yamlFile << "      y_max: " << 0.95 << "\n";
		yamlFile << "      W2_min: " << 10.0 << "\n";
		yamlFile << "      W2_max: " << 10000.0 << "\n";
		yamlFile << "      Target_mass_GeV: " << 0.938 << "\n";
		yamlFile << "      Product_mass_GeV: " << prod_mass << "\n";
		yamlFile << "      Num_gen: " << N_gen << "\n";
		yamlFile << "      Num_reco: " << N_rec << "\n";
		yamlFile << "      Num_reco_matched: " << N_rec_matched << "\n";
		yamlFile << "      N_reco_rel_stat_error: " << N_rec_stat_error << "\n";
		yamlFile << "      N_reco_abs_sys_error: " << N_rec_sys_error << "\n";
		yamlFile << "      efficiency: " << eff << "\n";
		yamlFile << "      N_DIS: " << N_DIS << "\n";
		yamlFile << "      N_DIS_rel_stat_error: " << N_DIS_stat_error << "\n";
		yamlFile << "      N_DIS_abs_sys_error: " << N_DIS_sys_error << "\n";
		        
        // CSV
        csvFile << n_id++ << "," << 4.0 * beam_e * beam_p << "," << mean_Q << "," << Qmin << "," << Qmax << "," << mean_xB << "," << xBmin << "," << xBmax << "," << 
        mean_z << "," << zmin << "," << zmax << "," << mean_Pt << "," << ptmin << "," << ptmax << "," << delta << "," << sigma << "," << sigma_stat_proj << "," << sigma_sys << "," << 
        0.01 << "," << 0.95 << "," << 10.0 << "," << 10000.0 << "," << 0.938 << "," << prod_mass << "," << N_gen << "," << N_rec << "," << N_rec_matched << "," << N_rec_stat_error << "," << N_rec_sys_error << "," << eff << "," << N_DIS << ","<< N_DIS_stat_error << "," << N_DIS_sys_error << "\n";
   

    }

    //outFile.Write();
    outFile.Close();
    //chain.Close();
    TString o_name = Form("%s_%dx%d_26_06", tag.Data(), beam_e, beam_p);
    cout << "-------------------------------------------" << endl;
    cout << "ROOT output file: " << outputFile << endl;
}
