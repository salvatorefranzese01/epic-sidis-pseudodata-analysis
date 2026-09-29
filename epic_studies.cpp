//-------------- Creation of Trees for electron and hadron
//--- Authors: Lorenzo Polizzi (lorenzo.polizzi@unife.it), Sara Pucillo (sara.pucillo@cern.ch), Nicolò Valle (nicolo.valle@cern.ch)
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TAxis.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <set>
#include <vector>
#include <TMath.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TPaveStatsEditor.h>
#include "TLorentzVector.h"
#include <cmath>
#include <random>
#include <fstream>

//------------------------ Functions ------------------------

//------------- Generation of logarithmically spaced bin edges between xmin and xmax
std::vector<double> CreateLogBinning(int nbins, double xmin, double xmax) {
    std::vector<double> bin_edges(nbins + 1);
    double logxmin = std::log10(xmin);
    double logxmax = std::log10(xmax);
    double bin_width = (logxmax - logxmin) / nbins;
    for (int i = 0; i <= nbins; ++i) {
        bin_edges[i] = std::pow(10, logxmin + i * bin_width);
    }
    return bin_edges;
}

//------------- Main function
void epic_studies(const char* fileList){

    //--- Open and read the input file
    std::ifstream inputFile(fileList);
    if (!inputFile.is_open()) {
        std::cerr << "Error: unable to open the list file!" << fileList << std::endl;
        return;
    }
    std::string outputFile;
    std::vector<std::string> fileNames;
    std::string line;
    bool firstLine = true;
    while (std::getline(inputFile, line)) {
        if (firstLine) {
            outputFile = line;          // first line: name of the output file
            firstLine = false;
        } else {
            fileNames.push_back(line);  // other lines: input files
        }
    }
    inputFile.close();
    
    std::cout << "Output: " << outputFile << std::endl;
std::cout << "Number of input files: " << fileNames.size() << std::endl;

for (auto &f : fileNames)
    std::cout << "INPUT=[" << f << "]" << std::endl;

    //--- Set up input file chain
    TChain *mychain = new TChain("events");
    for (const auto& file : fileNames) {
        mychain->Add(file.c_str());
    }
    // Check that the input chain contains events.
    if (mychain->GetEntries() == 0) {
        std::cerr << "Error: No events found in the chain!" << std::endl;
        return;
    }
    //--- Recreate the output file
    TFile *ofile = TFile::Open(outputFile.c_str(), "RECREATE");

    // Initialize the tree reader.
    TTreeReader tree_reader(mychain);
    if (mychain->GetEntries() == 0) {
        std::cerr << "Error: No events found in the chain!" << std::endl;
        return;
    }

    //--- Get Particle Information
    TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
    TTreeReaderArray<double> partMomX(tree_reader, "MCParticles.momentum.x");
    TTreeReaderArray<double> partMomY(tree_reader, "MCParticles.momentum.y");
    TTreeReaderArray<double> partMomZ(tree_reader, "MCParticles.momentum.z");
    TTreeReaderArray<double> partMass(tree_reader, "MCParticles.mass");
    TTreeReaderArray<int> partPdg(tree_reader, "MCParticles.PDG");
    TTreeReaderArray<int> parentsIndex(tree_reader, "_MCParticles_parents.index");
    TTreeReaderArray<int> daughterIndex(tree_reader, "_MCParticles_daughters.index");
    TTreeReaderArray<unsigned int> par(tree_reader, "MCParticles.parents_end");

    //--- Get Detector Information
    //TTreeReaderArray<float> dRICHx(tree_reader, "_DRICHAerogelTracks_points.position.x");
    //TTreeReaderArray<float> dRICHy(tree_reader, "_DRICHAerogelTracks_points.position.y");
    //TTreeReaderArray<float> dRICHz(tree_reader, "_DRICHAerogelTracks_points.position.z");
    //TTreeReaderArray<float> dRICH_momX(tree_reader, "_DRICHAerogelTracks_points.momentum.x");
    //TTreeReaderArray<float> dRICH_momY(tree_reader, "_DRICHAerogelTracks_points.momentum.y");
    //TTreeReaderArray<float> dRICH_momZ(tree_reader, "_DRICHAerogelTracks_points.momentum.z");
    //TTreeReaderArray<float> dRICH_theta(tree_reader, "_DRICHAerogelTracks_points.theta");
    //TTreeReaderArray<float> dRICH_phi(tree_reader, "_DRICHAerogelTracks_points.phi");

    //--- Get Reconstructed Track Information
    TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedRealPIDParticles.momentum.x");
    TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedRealPIDParticles.momentum.y");
    TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedRealPIDParticles.momentum.z");
    TTreeReaderArray<float> trackMass(tree_reader, "ReconstructedChargedParticles.mass");
    TTreeReaderArray<int> recPdg(tree_reader, "ReconstructedChargedParticles.PDG");
    TTreeReaderArray<float> goodnessOfPID(tree_reader, "ReconstructedChargedRealPIDParticles.goodnessOfPID");
    TTreeReaderArray<float> chi2(tree_reader, "CentralCKFTracks.chi2");

    //--- Get Associations Between MCParticles and ReconstructedChargedParticles
    TTreeReaderArray<unsigned int> scatElAssoc(tree_reader, "MCScatteredElectronAssociations_objIdx.collectionID"); 
    //TTreeReaderArray<unsigned int> simuAssoc(tree_reader, "ReconstructedChargedParticleAssociations_sim.index");
    TTreeReaderArray<unsigned int> recoAssoc(tree_reader, "_ReconstructedChargedParticleAssociations_rec.index");
    TTreeReaderArray<unsigned int> simuAssoc(tree_reader, "_ReconstructedChargedParticleAssociations_sim.index");

    //--- Graphs information
    const int nbins = 120;
    const double xmin_xB = 1e-4, xmax_xB = 1;
    const double xmin_Q2 = 1, xmax_Q2 = 100.;
    const double xmin_a = 1, xmax_a = 10000.;
    auto make_bins = [](int bins, double min, double max) {
        return CreateLogBinning(bins, min, max);
      };
    const auto log_bins_Q2 = make_bins(nbins, xmin_Q2, xmax_Q2);
    const auto log_bins_xB = make_bins(nbins, xmin_xB, xmax_xB);
    const auto log_bins_pr = make_bins(nbins, xmin_a, xmax_a);

    //--- Variables
    //--- electron
    double el_px_mc, el_py_mc, el_pz_mc, el_mom_mc, el_y_mc, el_theta_mc, el_phi_mc, el_eta_mc;
    double el_px, el_py, el_pz, el_mom, el_pdg, el_theta, el_phi, el_eta, el_y;
    //--- proton
    double pr_mom_mc, pr_px_mc, pr_py_mc, pr_pz_mc, pr_phi_mc, pr_theta_mc, pr_eta_mc;
    double pr_mom, pr_px, pr_py, pr_pz, pr_phi, pr_theta, pr_eta;
    //--- hadron
    double hadron_px_mc, hadron_py_mc, hadron_pz_mc, hadron_mom_mc, hadron_Q2_mc, hadron_xB_mc, hadron_xF_mc, hadron_z_mc, hadron_PhT_mc;
    double hadron_Phi_lab_mc, hadron_Phi_h_mc, hadron_Phi_s_mc, hadron_Theta_mc, hadron_eta_mc, hadron_y_mc, hadron_W_mc, hadron_Mx_mc;
    double hadron_pdg_mc;
    double hel_mc, eps_mc;
    double hadron_px, hadron_py, hadron_pz, hadron_mom, hadron_Q2, hadron_xB, hadron_xF, hadron_z, hadron_PhT;
    double hadron_Phi_lab, hadron_Phi_h, hadron_Phi_s, hadron_Theta, hadron_eta, hadron_y, hadron_W, hadron_Mx;
    double hel, eps, hadron_pdg;
    double hadron_goodPID, hadron_index, hadron_index_mc;
    //int hadron_count;
    

    //--- Trees
    //--- MC
    TTree ElectronTreeMC("ElectronTree_MC", "");
    ElectronTreeMC.Branch("el_px_mc", &el_px_mc, "el_px_mc/D");
    ElectronTreeMC.Branch("el_py_mc", &el_py_mc, "el_py_mc/D");
    ElectronTreeMC.Branch("el_pz_mc", &el_pz_mc, "el_pz_mc/D");
    ElectronTreeMC.Branch("el_mom_mc",&el_mom_mc, "el_mom_mc/D");
    ElectronTreeMC.Branch("el_y_mc", &el_y_mc, "el_y_mc/D");
    ElectronTreeMC.Branch("el_theta_mc", &el_theta_mc, "el_theta_mc/D");
    ElectronTreeMC.Branch("el_phi_mc", &el_phi_mc, "el_phi_mc/D");
    ElectronTreeMC.Branch("el_eta_mc", &el_eta_mc, "el_eta_mc/D");

    TTree HadronTreeMC("HadronTree_MC", "");
    HadronTreeMC.Branch("hadron_index_mc", &hadron_index_mc, "hadron_index_mc/D");
    HadronTreeMC.Branch("hadron_pdg_mc", &hadron_pdg_mc, "hadron_pdg_mc/D");
    HadronTreeMC.Branch("hadron_px_mc", &hadron_px_mc, "hadron_px_mc/D");
    HadronTreeMC.Branch("hadron_py_mc", &hadron_py_mc, "hadron_py_mc/D");
    HadronTreeMC.Branch("hadron_pz_mc", &hadron_pz_mc, "hadron_pz_mc/D");
    HadronTreeMC.Branch("hadron_mom_mc", &hadron_mom_mc, "hadron_mom_mc/D");
    HadronTreeMC.Branch("hadron_Q2_mc", &hadron_Q2_mc, "hadron_Q2_mc/D");
    HadronTreeMC.Branch("hadron_xB_mc", &hadron_xB_mc, "hadron_xB_mc/D");
    HadronTreeMC.Branch("hadron_xF_mc", &hadron_xF_mc, "hadron_xF_mc/D");
    HadronTreeMC.Branch("hadron_z_mc", &hadron_z_mc, "hadron_z_mc/D");
    HadronTreeMC.Branch("hadron_PhT_mc", &hadron_PhT_mc, "hadron_PhT_mc/D");
    HadronTreeMC.Branch("hadron_Phi_lab_mc", &hadron_Phi_lab_mc, "hadron_Phi_lab_mc/D");
    HadronTreeMC.Branch("hadron_Phi_h_mc", &hadron_Phi_h_mc, "hadron_Phi_h_mc/D");
    //HadronTreeMC.Branch("hadron_Phi_s_mc", &hadron_Phi_s_mc, "hadron_Phi_s_mc/D");
    HadronTreeMC.Branch("hadron_Theta_mc", &hadron_Theta_mc, "hadron_Theta_mc/D");
    HadronTreeMC.Branch("hadron_eta_mc", &hadron_eta_mc, "hadron_eta_mc/D");
    HadronTreeMC.Branch("hadron_y_mc", &hadron_y_mc, "hadron_y_mc/D");
    //HadronTreeMC.Branch("hadron_W_mc", &hadron_W_mc, "hadron_W_mc/D");
    //HadronTreeMC.Branch("hadron_Mx_mc", &hadron_Mx_mc, "hadron_Mx_mc/D");
    //HadronTreeMC.Branch("helicity_mc", &hel_mc, "helicity_mc/D");
    //HadronTreeMC.Branch("epsilon_mc", &eps_mc, "epsilon_mc/D");

    //--- Reco
    TTree ElectronTreeRECO("ElectronTree_RECO", "");
    ElectronTreeRECO.Branch("el_px", &el_px, "el_px/D");
    ElectronTreeRECO.Branch("el_py", &el_py, "el_py/D");
    ElectronTreeRECO.Branch("el_pz", &el_pz, "el_pz/D");
    ElectronTreeRECO.Branch("el_mom", &el_mom, "el_mom/D");
    ElectronTreeRECO.Branch("el_pdg", &el_pdg, "el_pdg/D");
    ElectronTreeRECO.Branch("el_y", &el_y, "el_y/D");
    ElectronTreeRECO.Branch("el_theta", &el_theta, "el_theta/D");
    ElectronTreeRECO.Branch("el_phi", &el_phi, "el_phi/D");
    ElectronTreeRECO.Branch("el_eta", &el_eta, "el_eta/D");

    TTree HadronTreeRECO("HadronTree_RECO", "");
    HadronTreeRECO.Branch("hadron_index", &hadron_index, "hadron_index/D");
    HadronTreeRECO.Branch("hadron_px", &hadron_px, "hadron_px/D");
    HadronTreeRECO.Branch("hadron_py", &hadron_py, "hadron_py/D");
    HadronTreeRECO.Branch("hadron_pz", &hadron_pz, "hadron_pz/D");
    HadronTreeRECO.Branch("hadron_pdg", &hadron_pdg, "hadron_pdg/D");
    HadronTreeRECO.Branch("hadron_pdg_mc", &hadron_pdg_mc, "hadron_pdg_mc/D");
    HadronTreeRECO.Branch("hadron_good_PID", &hadron_goodPID, "hadron_good_PID/D");
    HadronTreeRECO.Branch("hadron_mom", &hadron_mom, "hadron_mom/D");
    HadronTreeRECO.Branch("hadron_Q2", &hadron_Q2, "hadron_Q2/D");
    HadronTreeRECO.Branch("hadron_xB", &hadron_xB, "hadron_xB/D");
    HadronTreeRECO.Branch("hadron_xF", &hadron_xF, "hadron_xF/D");
    HadronTreeRECO.Branch("hadron_z", &hadron_z, "hadron_z/D");
    HadronTreeRECO.Branch("hadron_PhT", &hadron_PhT, "hadron_PhT/D");
    HadronTreeRECO.Branch("hadron_Phi_lab", &hadron_Phi_lab, "hadron_Phi_lab/D");
    HadronTreeRECO.Branch("hadron_Phi_h", &hadron_Phi_h, "hadron_Phi_h/D");
    //HadronTreeRECO.Branch("hadron_Phi_s", &hadron_Phi_s, "hadron_Phi_s/D");
    HadronTreeRECO.Branch("hadron_Theta", &hadron_Theta, "hadron_Theta/D");
    HadronTreeRECO.Branch("hadron_eta", &hadron_eta, "hadron_eta/D");
    HadronTreeRECO.Branch("hadron_y", &hadron_y, "hadron_y/D");
    //HadronTreeRECO.Branch("hadron_W", &hadron_W, "hadron_W/D");
    //HadronTreeRECO.Branch("hadron_Mx", &hadron_Mx, "hadron_Mx/D");
    //HadronTreeRECO.Branch("helicity", &hel, "helicity/D");
    //HadronTreeRECO.Branch("epsilon", &eps, "epsilon/D");
    
    // Store matched generated-hadron kinematics in the reconstructed tree for resolution studies.
    HadronTreeRECO.Branch("hadron_Q2_mc",  &hadron_Q2_mc);
		HadronTreeRECO.Branch("hadron_xB_mc",  &hadron_xB_mc);
		HadronTreeRECO.Branch("hadron_z_mc",   &hadron_z_mc);
		HadronTreeRECO.Branch("hadron_PhT_mc", &hadron_PhT_mc);
		HadronTreeRECO.Branch("hadron_y_mc",   &hadron_y_mc);


    //double count3 = 0;
    double event = -1;
    //--- output for debug
    std::ofstream csv("output.csv");
    csv << "event , status, parents, reco_pdg, true_pdg, px_reco, py_reco, pz_reco, mom_reco, px_true, py_true, pz_true, mom_true\n";

    //------------- pdg list
    //--- 11 electron
    //--- -211 negative pion
    //--- 211 positive pion
    //--- 321 positive kaon
    //--- -321 negative kaon
    //--- 2212 proton

    //------------- Pythia info
    // status = 4 is the beam (ref 1767) in HepMC
    // status = 1 stable particle in the final state

    //------------- Loop over events

    TLorentzVector ElectronScattered;
    TLorentzVector ElectronScattered_mc;
    TLorentzVector ProtonScattered;
    TLorentzVector ProtonScattered_mc;
    TLorentzVector MC_ProtonBeam;
    TLorentzVector MC_ElectronBeam;
    

    while(tree_reader.Next()) {

        //--- Useful vectors or variables
        std::vector<TVector3> recScatElectron;
        /* dead code, i get the energies from the beam directly 
        TVector3 ElBeam(0.,0.,-18.);
        const double m_p = 0.9382720813;  // GeV
        const double m_e = 0.00051099895; // GeV
        double pz_e = -10.0;
        double p_e = std::abs(pz_e);
        double E_e = sqrt(p_e*p_e + m_e*m_e);
        double pz_p = 100.0;
        double p_p = std::abs(pz_p);
        double E_p = sqrt(p_p*p_p + m_p*m_p);
        TLorentzVector ElectronBeam(0., 0., pz_e, E_e);
        TLorentzVector ProtonBeam(0., 0., pz_p, E_p);
        TLorentzVector Lab = ElectronBeam + ProtonBeam;
        */
        double currentPhi = 0;
        double currentMom = 0;
        TVector3 ScatElectron_mc;
        //double count2 = 0;
        //count3++;
        //double count = 0;
        //if (event > 20) continue; // used to look at small data
        //--- Loop over generated particles
        for(unsigned int i=0; i<partGenStat.GetSize(); i++){
            int pdg = (partPdg[i]);
            int pdg2 = (std::abs(partPdg[i])); //in order to consider both particle and antiparticle
            // status = 4 is the beam (ref 1767) in HepMC
            // status = 1 stable particle in the final state
            if(partGenStat[i] != -1){    // stable
                if(parentsIndex[i] > -1){
                    if (partGenStat[i] == 4 && pdg == 2212){ //beam+proton
                        TVector3 Mom(partMomX[i],partMomY[i],partMomZ[i]);
                        double E = sqrt(0.9382720813*0.9382720813 + partMomX[i]*partMomX[i] + partMomY[i]*partMomY[i] + partMomZ[i]*partMomZ[i]);
                        MC_ProtonBeam.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], E);
                        event++;
                        //csv << "----------------------------- event " << event << " ----------------------------" << endl;
                    }
                    if(partGenStat[i] == 4 && pdg == 11){ //beam+electron
                        TVector3 Mom(partMomX[i],partMomY[i],partMomZ[i]);
                        double E = sqrt(0.00051099895*0.00051099895 + partMomX[i]*partMomX[i] + partMomY[i]*partMomY[i] + partMomZ[i]*partMomZ[i]);
                        MC_ElectronBeam.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], E);
                    }
                    if(pdg2 == 11 || pdg2 == 211 || pdg2 == 321 || pdg2 == 2212){ //all hadrons
                        for(unsigned int j=0; j<recPdg.GetSize(); j++){
                            if(simuAssoc[j] == i){
                                TVector3 ElMom(partMomX[i],partMomY[i],partMomZ[i]);
                                TVector3 MomReco(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]);
                                csv << event << ",  " << partGenStat[i] << ",  "  << parentsIndex[i] << ",  " << recPdg[j] << ",  " << pdg << ",  " << MomReco.X() << ",  " << MomReco.Y() << ",  " << MomReco.Z() << ",  " << MomReco.Mag() << ",  " << ElMom.X() << ",  " << ElMom.Y()  << ",  " << ElMom.Z() << ",  " << ElMom.Mag()  << "\n";
                                //count++;
                            }
                        }
                    }
                }
                if(pdg == 11 || pdg2 == 211 || pdg2 == 321 || pdg2 == 2212){ //all hadrons
                    //count2++;
                    if(parentsIndex[i] > -1){
                        //--- loops on the reconstructed particles and associations with the simulated particle
                        for(unsigned int j=0; j<recPdg.GetSize(); j++){
                            if(simuAssoc[j] == i){ // Find association index matching the index of the generated particle we are looking at
                                if (pdg == 11 ){ //electron
                                    TVector3 ElMom(partMomX[i],partMomY[i],partMomZ[i]);
                                    el_eta_mc = ElMom.PseudoRapidity();
                                    double angleR = ElMom.Theta();
                                    double angle = angleR * (180.0 / TMath::Pi());
                                    el_px_mc = ElMom.X();
                                    el_py_mc = ElMom.Y();
                                    el_pz_mc = ElMom.Z();
                                    el_theta_mc = ElMom.Theta();
                                    el_phi_mc = ElMom.Phi();
                                    el_mom_mc = ElMom.Mag();
                                    currentPhi = ElMom.Theta();
                                    currentMom = ElMom.Mag();
                                    el_pdg = recPdg[j];
                                    ScatElectron_mc.SetXYZ(partMomX[i],partMomY[i],partMomZ[i]);
                                    ElectronScattered_mc.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], el_mom_mc);
                                    TLorentzVector el_q_mc = MC_ElectronBeam - ElectronScattered_mc;
                                    el_y_mc = (MC_ProtonBeam.Dot(el_q_mc))/(MC_ProtonBeam.Dot(MC_ElectronBeam));

                                    // Compute reconstructed electron kinematics before applying the electron selection.
                                    int recpdg = (recPdg[j]);
                                    TVector3 recElmom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]);
                                    double eta = recElmom.PseudoRapidity();
                                    el_px = recElmom.X(), el_py = recElmom.Y(), el_pz = recElmom.Z();
                                    el_mom = recElmom.Mag();
                                    double theta = recElmom.Theta();
                                    double phi = recElmom.Phi();
                                    // Store reconstructed electron angles and pseudorapidity.
                                    el_eta = eta;
                                    el_theta = theta;
                                    el_phi = phi;
                                    ElectronScattered.SetPxPyPzE(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]], el_mom);
                                    TLorentzVector el_q = MC_ElectronBeam - ElectronScattered;
                                    el_y = (MC_ProtonBeam.Dot(el_q))/(MC_ProtonBeam.Dot(MC_ElectronBeam));

                                    // Select generated electrons within the kinematic and acceptance ranges.
                                    if(el_mom_mc >= 1 && el_eta_mc >= -3.5 && el_eta_mc <= 3.5 && el_y_mc >= 0.01 && el_y_mc <= 0.99){
                                        //---filling
                                        ElectronTreeMC.Fill();
                                    }
                                    // Select reconstructed electrons using reconstructed kinematics.
                                    if(el_mom >= 1 && el_eta >= -3.5 && el_eta <= 3.5 && el_y >= 0.01 && el_y <= 0.99){
                                        //---filling
                                        ElectronTreeRECO.Fill();
                                    }
                                } if (std::abs(recPdg[j]) == 211 || std::abs(recPdg[j]) == 321 || std::abs(recPdg[j]) == 2212){ //important request
                                    TVector3 particle(partMomX[i],partMomY[i],partMomZ[i]);
                                    TVector3 reco_particle(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]);
                                    double eta_particle = particle.PseudoRapidity();
                                    double eta = reco_particle.PseudoRapidity();
                                    hadron_pdg_mc = partPdg[i];
                                    hadron_index = i;
                                    // Apply the hadron acceptance and reconstructed-electron y selection.
                                    if(eta >= -3.5 && eta <= 3.5 && el_y >= 0.01 && el_y <= 0.99){
                                        hadron_pdg = (recPdg[j]);
                                        double goodPID_cont = goodnessOfPID[j];
                                        hadron_mom = reco_particle.Mag();
                                        // Virtual photon from the reconstructed scattered electron.
                                        TLorentzVector photon_hadron_noBoost = MC_ElectronBeam - ElectronScattered;
                                        double mass = partMass[i];
                                        double E_hadron = sqrt(hadron_mom*hadron_mom + mass*mass);
                                        TLorentzVector hadron_noBoost(reco_particle.X(), reco_particle.Y(), reco_particle.Z(), E_hadron);
                                        hadron_y = el_y;
                                        hadron_px = reco_particle.X(); hadron_py = reco_particle.Y(); hadron_pz = reco_particle.Z();
                                        hadron_index++;
                                        hadron_mom = reco_particle.Mag();
                                        hadron_Q2 = -photon_hadron_noBoost.M2();
                                        hadron_xB = hadron_Q2 / (2 * MC_ProtonBeam.Dot(photon_hadron_noBoost));
                                        hadron_eta = eta;
                                        hadron_Theta = reco_particle.Theta(); hadron_Phi_lab = reco_particle.Phi();
                                        hadron_z = (MC_ProtonBeam * hadron_noBoost) / (MC_ProtonBeam * photon_hadron_noBoost);
                                        hadron_goodPID = goodPID_cont;
                                        
                                        // Store the corresponding generated-hadron kinematics.
                                        TLorentzVector photon_hadron_noBoost_mc = MC_ElectronBeam - ElectronScattered_mc;
                                        double E_hadron_mc = sqrt(particle.Mag()*particle.Mag() + mass*mass);
                                        TLorentzVector hadron_noBoost_mc(particle.X(), particle.Y(), particle.Z(), E_hadron_mc);                                       
                                        
                                        hadron_Q2_mc = -photon_hadron_noBoost_mc.M2();
                                        hadron_xB_mc = hadron_Q2_mc / (2. * MC_ProtonBeam.Dot(photon_hadron_noBoost_mc));
                                        hadron_z_mc = (MC_ProtonBeam * hadron_noBoost_mc) / (MC_ProtonBeam * photon_hadron_noBoost_mc);
                                        
                                        //if (hadron_Q2 < 1) std::cout << "ev: " << event << "   mom: " << hadron_mom << "  el: " << ElectronScattered_mc.Px() << ", "  << ElectronScattered_mc.Py() << ", " <<  ElectronScattered_mc.Pz() << ", " << ElectronScattered_mc.E()
                                        //<< "    " << MC_ElectronBeam.Px() << ", " << MC_ElectronBeam.Py() << ", " << MC_ElectronBeam.Pz() << ", " << MC_ElectronBeam.E() << std::endl;
                                        TLorentzVector hadron_4vec = hadron_noBoost;
                                        TLorentzVector photon_hadron_4vec = photon_hadron_noBoost;
                                        // boost gamma*N - mandatory for the extraction of P_hT, Phi_h and Phi_s
                                        TLorentzVector gammaN = photon_hadron_noBoost + MC_ProtonBeam;
                                        TVector3 boost_gammaN = -gammaN.BoostVector();
                                        hadron_4vec.Boost(boost_gammaN);
                                        photon_hadron_4vec.Boost(boost_gammaN);
                                        TVector3 zAxis = photon_hadron_4vec.Vect().Unit();
                                        hadron_PhT = hadron_4vec.Perp(zAxis);
                                        // calculation of Phi_h
                                        TVector3 yAxis = (MC_ElectronBeam.Vect().Cross(ElectronScattered_mc.Vect()));
                                        TVector3 xAxis = yAxis.Cross(zAxis);
                                        TVector3 PhT_vector = hadron_4vec.Vect() - (hadron_4vec.Vect().Dot(zAxis))*zAxis;
                                        
                                        // Compute generated transverse momentum in the photon-proton frame.
																				TLorentzVector hadron_4vec_mc = hadron_noBoost_mc;
																				TLorentzVector photon_hadron_4vec_mc = photon_hadron_noBoost_mc;
																				TLorentzVector gammaN_mc = photon_hadron_noBoost_mc + MC_ProtonBeam;
																				TVector3 boost_gammaN_mc = -gammaN_mc.BoostVector();
																				hadron_4vec_mc.Boost(boost_gammaN_mc);
																				photon_hadron_4vec_mc.Boost(boost_gammaN_mc);
																				TVector3 zAxis_mc = photon_hadron_4vec_mc.Vect().Unit();
																				hadron_PhT_mc = hadron_4vec_mc.Perp(zAxis_mc);

																				// Calculate the reconstructed hadron azimuthal angle.
                                        double hadron_Phx = PhT_vector.Dot(xAxis);
                                        double hadron_Phy = PhT_vector.Dot(yAxis);
                                        hadron_Phi_h = TMath::ATan2(hadron_Phy, -hadron_Phx);
                                        // let's save also y_mc so we can compute the resolution for y too
                                        hadron_y_mc = (MC_ProtonBeam.Dot(photon_hadron_noBoost_mc)) / (MC_ProtonBeam.Dot(MC_ElectronBeam));
                                        
                                        //---filling
                                        HadronTreeRECO.Fill();
                                    }
                                }
                            }
                        } // no reco cycle
                        if ((pdg2 == 211 || pdg2 == 321 || pdg2 == 2212) && partGenStat[i] == 1){
                            TVector3 mc_particle(partMomX[i],partMomY[i],partMomZ[i]);
                            hadron_index_mc = i;
                            hadron_px_mc = mc_particle.X(), hadron_py_mc = mc_particle.Y(), hadron_pz_mc = mc_particle.Z();
                            hadron_pdg_mc = partPdg[i], hadron_mom_mc = mc_particle.Mag();
                            TLorentzVector photon_hadron_noBoost_mc = MC_ElectronBeam - ElectronScattered_mc;
                            double mass = partMass[i];
                            double E_hadron_mc = sqrt(hadron_mom_mc*hadron_mom_mc + mass*mass);
                            TLorentzVector hadron_noBoost_mc(mc_particle.X(), mc_particle.Y(), mc_particle.Z(), E_hadron_mc);
                            hadron_Q2_mc = -photon_hadron_noBoost_mc.M2();
                            hadron_xB_mc = hadron_Q2_mc / (2 * MC_ProtonBeam.Dot(photon_hadron_noBoost_mc));
                            hadron_eta_mc = mc_particle.PseudoRapidity(), hadron_Theta_mc = mc_particle.Theta(); hadron_Phi_lab_mc = mc_particle.Phi();
                            hadron_z_mc = (MC_ProtonBeam * hadron_noBoost_mc) / (MC_ProtonBeam * photon_hadron_noBoost_mc);
                            //
                            TLorentzVector hadron_4vec_mc = hadron_noBoost_mc;
                            TLorentzVector photon_hadron_4vec_mc = photon_hadron_noBoost_mc;
                            // boost gamma*N - mandatory for the extraction of P_hT, Phi_h and Phi_s
                            TLorentzVector gammaN_mc = photon_hadron_noBoost_mc + MC_ProtonBeam;
                            TVector3 boost_gammaN_mc = -gammaN_mc.BoostVector();
                            hadron_4vec_mc.Boost(boost_gammaN_mc);
                            photon_hadron_4vec_mc.Boost(boost_gammaN_mc);
                            TVector3 zAxis_mc = photon_hadron_4vec_mc.Vect().Unit();
                            hadron_PhT_mc = hadron_4vec_mc.Perp(zAxis_mc);
                            //
                            hadron_y_mc = (MC_ProtonBeam.Dot(photon_hadron_noBoost_mc))/(MC_ProtonBeam.Dot(MC_ElectronBeam));
                            if(hadron_eta_mc >= -3.5 && hadron_eta_mc <= 3.5 && hadron_Q2_mc >= 1 && hadron_z_mc < 1 && hadron_z_mc > 0 && hadron_xB_mc > 0 && hadron_xB_mc < 1) HadronTreeMC.Fill();
                        }
                    }
                }
            } //end of the loop over stable generated particles
        } //end of the loop over generated particles
    } //end of the loop over events
    
    
    
    


    //ElectronTreeMC.Write();
    //HadronTreeMC.Write();
    //ElectronTreeRECO.Write();
    //HadronTreeRECO.Write();

    std::cout << " " << std::endl;
    std::cout << "ROOT output file: " << outputFile << endl;
    std::cout << "______________________________________________________________________________________" << std::endl;
    std::cout << " " << std::endl;


    ofile->Write();
    ofile->Close();

    mychain->Delete();
    ofile->Delete();
  }
