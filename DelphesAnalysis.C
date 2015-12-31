#define DelphesAnalysis_cxx

#include "DelphesAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <vector>
#include <TLorentzVector.h>
#include <string>
#include <iostream>
void DelphesAnalysis::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L DelphesAnalysis.C
//      Root > DelphesAnalysis t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by   b_branchname->GetEntry(ientry); //read only this branch
TH1D *R_b1_j1 = new TH1D("deltaR_b1_j1","Delta R of b1 and j1",50,0,0.8);
TH1D *R_b2_j1 = new TH1D("deltaR_b2_j1","Delta R of b2 and j1",50,0,0.8);
TH1D *R_b1_j2 = new TH1D("deltaR_b1_j2","Delta R of b1 and j2",50,0,0.8);
TH1D *R_b2_j2 = new TH1D("deltaR_b2_j2","Delta R of b2 and j2",50,0,0.8);
TH1D *ZMassMatched = new TH1D("ZmassMatched","Z mass with Truth Match",50,30,150);   
TH1D *h2 = new TH1D("muonPt","LL Pt",50,0,300);
TH1D *Wmass = new TH1D("Wmass","W Transverse Mass ",60,0,180);
TH1D *Zmass= new TH1D("Zmass"," Z mass ",60,0,180);
//TH1D *Zmass= new TH1D("ZmassGenJet","Z mass ",60,0,180);
TH1D *ZmassGen= new TH1D("ZmassGen", " Z mass",60,0,180);
TH1D *WmassCut = new TH1D("WMass","W mass after cuts",60,0,180);
TH1D *integratedZ = new TH1D("IntegratedZ"," fully cut Z mass" ,30,50,105);
TH1D *integratedW = new TH1D("IntegratedW"," fully cut W mass" ,30,55,90);
TH1D *numJetsAfterCuts = new TH1D("NumJetsAfterCuts","Number of jets after cuts",15,0,15);
TH1D *numbJets = new TH1D("numbJets","Number of of b jets",5,0,5);

 int muonCount;
 int bCount;
int bJetCount;
vector <int> index;
TH1D *h1 = new TH1D("Met","Missing Et",50,0,300);
TCanvas *c1 = new TCanvas("c1","MET",700,400);
int integrator=0;

std::string line;

double dphiZ;
double dphiMix;
vector< vector<double> > yMatrix;
vector< vector<double> > phiMatrix;
vector< vector<double> > IDMatrix;
vector< vector<double> > MIDMatrix;
double ZMassMatrix[10000];
double jetMatrix[10000][2][2] ;
//The above will be my selected jets in question. First index will be the event, second will
//indicate the jet (first or second), and the last two will be y and phi in that order.

yMatrix.resize(10000,std::vector<double>(12,0));
phiMatrix.resize(10000,std::vector<double>(12,0));
IDMatrix.resize(10000,std::vector<double>(12,0));
MIDMatrix.resize(10000,std::vector<double>(12,0));

if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();   

  Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

     Long64_t ientry = LoadTree(jentry);
//cout << " FRESH EVENT " << jentry << endl ;
	muonCount=0;
	bJetCount=0;
	bCount=0;

	if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

index.clear();

	for ( int j =0;j<Jet_size;j++)//This is the Jet loop
{	
		if(Jet_BTag[j]==1){
		bJetCount++;
		index.push_back (j);
//		cout << " this is bjet number " << j <<" out of " << Jet_size << " it has Pt= " << Jet_PT[j] << endl;
				} 
}


double dumy=0.0;

	for(int i =0;i<Muon_size;i++)//this is the Muon Loop
	{	
		if(Muon_PT[i]>dumy) 
		 dumy = Muon_PT[i];
		
	}
//cout << "testline 7 " << endl;
if(dumy!=0.0)
 h2->Fill(dumy);
 h1->Fill(MissingET_MET[0]); 



//These are the prelimenary cuts that preceed the jet sorting
//There will be more cuts after the jet sorting that are dependent on this process
if(Muon_size>=1){
if(bJetCount>=2){
if(Jet_size>2){
TLorentzVector j1;
TLorentzVector j2;

int highestNobPT;

if (index[0]!=0 && index[1]!=0)
	 highestNobPT=0;

else if(index[0]!=1 && index[1]!=1)
	 highestNobPT=1;

else highestNobPT=2;

std::vector<int> sortJets ;
sortJets.push_back(index[0]);
sortJets.push_back(index[1]);
sortJets.push_back(highestNobPT);

//The ram variable is an memory location used for the sorting of the sortJets array
int ram;
for(int i=0;i<sortJets.size();i++)
{
for(int j=0;j<sortJets.size();j++)
	{
		if(sortJets[j]>sortJets[i])
		{
			ram=sortJets[i];
			sortJets[i]=sortJets[j];
			sortJets[j]=ram;
		}
	}

}
//At this point the array sortJets has the three Jets I like sorted by PT
//This declaration uses the jet sorting algorithm
j1.SetPtEtaPhiM(Jet_PT[sortJets[1]],Jet_Eta[sortJets[1]],Jet_Phi[sortJets[1]],Jet_Mass[sortJets[1]]);
j2.SetPtEtaPhiM(Jet_PT[sortJets[2]],Jet_Eta[sortJets[2]],Jet_Phi[sortJets[2]],Jet_Mass[sortJets[2]]);

//This jet declarations will simply take the second and third highest Pt jets in the event
//j2.SetPtEtaPhiM(Jet_PT[2],Jet_Eta[2],Jet_Phi[2],Jet_Mass[2]);
//j1.SetPtEtaPhiM(Jet_PT[1],Jet_Eta[1],Jet_Phi[1],Jet_Mass[1]);

//This is where I implement a rapidity cut and this is the best place for all subsequent added cuts
if(fabs(j1.Rapidity()-j2.Rapidity())<=0.9){
jetMatrix[jentry][0][0]=j1.Rapidity();
jetMatrix[jentry][1][0]=j2.Rapidity();
jetMatrix[jentry][0][1]=j1.Phi();
jetMatrix[jentry][1][1]=j2.Phi();

TLorentzVector sum;
sum= j1+j2;

//This ZMassMatrix is is just to catalogue the masses of the Z boson for use in the later steps where I 
//compare to the Truth data
ZMassMatrix[jentry]=sum.M();

//The cut W mass
//I calculate the W mass by using the kinematics of the Leading Muon and the Missing ET
//This could be done more streamlined with TLorentz but I did it for the excercise
double thetaMuon=2*atan(exp(-Muon_Eta[0]));
double thetaMET=2*atan(exp(-MissingET_Eta[0]));

double PxMuon=Muon_PT[0]*cos(Muon_Phi[0]);
double PyMuon=Muon_PT[0]*sin(Muon_Phi[0]);
double PzMuon=Muon_PT[0]/tan(thetaMuon);

double PxMET=MissingET_MET[0]*cos(MissingET_Phi[0]);
double PyMET=MissingET_MET[0]*sin(MissingET_Phi[0]);
double PzMET=MissingET_MET[0]/tan(thetaMET);

double MuonMag=sqrt(pow(PxMuon,2)+pow(PyMuon,2)+pow(PzMuon,2));
double METMag=sqrt(pow(PxMET,2)+pow(PyMET,2)+pow(PzMET,2));


double Wmt2= pow(Muon_PT[0]+MissingET_MET[0],2)-pow(PxMuon+PxMET,2)-pow(PyMuon+PyMET,2);
WmassCut->Fill(sqrt(Wmt2));

//This is the "Integrator" portion of the event where I count the number of events that fall on both the Z and W mass shells simultaneously
if(sum.M()<=100 && sum.M()>=55 && sqrt(Wmt2)<=85 && sqrt(Wmt2)>=60)
{
integrator++;
integratedW->Fill(sqrt(Wmt2));
integratedZ->Fill(sum.M());
}

Zmass->Fill(sum.M());


numJetsAfterCuts->Fill(Jet_size);
numbJets->Fill(bJetCount);

//These last brackets below correspond to the cutting criteria
}
}
}
}

//Begin reconstruction of W

//if(Muon_size>=1){
//cout << "this is where i Fail?" << endl;
//this section is where I plot the W tranverse mass but without the overall cuts
//Same analysis as above but without the cuts and different histogram
double thetaMuon=2*atan(exp(-Muon_Eta[0]));
double thetaMET=2*atan(exp(-MissingET_Eta[0]));

double PxMuon=Muon_PT[0]*cos(Muon_Phi[0]);
double PyMuon=Muon_PT[0]*sin(Muon_Phi[0]);
double PzMuon=Muon_PT[0]/tan(thetaMuon);

double PxMET=MissingET_MET[0]*cos(MissingET_Phi[0]);
double PyMET=MissingET_MET[0]*sin(MissingET_Phi[0]);
double PzMET=MissingET_MET[0]/tan(thetaMET);

double MuonMag=sqrt(pow(PxMuon,2)+pow(PyMuon,2)+pow(PzMuon,2));
double METMag=sqrt(pow(PxMET,2)+pow(PyMET,2)+pow(PzMET,2));


double Wmt2= pow(Muon_PT[0]+MissingET_MET[0],2)-pow(PxMuon+PxMET,2)-pow(PyMuon+PyMET,2);
Wmass->Fill(sqrt(Wmt2));

      // if (Cut(ientry) < 0) continue;
   }//EVENT LOoP CLOSER


//Here begins the Truth analysis
//I begin by reading a file in that contains the data from the truth,
//putting that information into some matricies,
//then comparing to the information I collected above in the event loops of the delphes file

ifstream file("TruthDeltaR.txt");
//  file = new ifstream("TruthDeltaR.txt");
  int eventNum=0;
  int partNum=0;
  double dumy=0;
  string thing;
  if(file.is_open()){
  while(file.good()){
  getline(file,line);
  if(file.eof())break;
  if(line == "E"){
  cout << " Truth Event Comparison " << endl;
  cout <<line << endl;
  getline(file,line);
  cout <<line<<endl;
  istringstream(line) >> eventNum;
  }
  if(line == "P"){
//  cout << "line contains = " << line <<endl;
  getline(file,line);
  istringstream(line) >> partNum;
//  cout << "line contains = " << line <<endl;
  getline(file,line);
  istringstream(line) >> IDMatrix[eventNum][partNum];
//  cout << "line contains = " << line <<endl;
  getline(file,line);
  istringstream(line) >> MIDMatrix[eventNum][partNum];
//  cout << "line contains = " << line <<endl;
  getline(file,line);
  istringstream(line) >> yMatrix[eventNum][partNum];
//  cout << "line contains = " << line <<endl;
  getline(file,line);
  istringstream(line) >> phiMatrix[eventNum][partNum];
//  cout << "line contains = " << line <<endl;
  }

}
}

double dR_b1_j1=0;
double dR_b2_j2=0;
double dR_b1_j2=0;
double dR_b2_j1=0;

double lowestR1=100;
double lowestR2=100;
int indexOfR1=0;
int indexOfR2=0;
int goodMatchCount=0;
for(int i =0;i<10000;i++){
cout << "next event comparison " << endl;


 
if(jetMatrix[i][0][0]!=0 && jetMatrix[i][0][1]!=0 && jetMatrix[i][1][0]!=0 && jetMatrix[i][1][1]!=0 ){
 dR_b1_j1 = sqrt(pow(yMatrix[i][8]-jetMatrix[i][0][0],2)+pow(phiMatrix[i][8]-jetMatrix[i][0][1],2));  
 dR_b1_j2 = sqrt(pow(yMatrix[i][8]-jetMatrix[i][1][0],2)+pow(phiMatrix[i][8]-jetMatrix[i][1][1],2));
 dR_b2_j1 = sqrt(pow(yMatrix[i][9]-jetMatrix[i][0][0],2)+pow(phiMatrix[i][9]-jetMatrix[i][0][1],2));
 dR_b2_j2 = sqrt(pow(yMatrix[i][9]-jetMatrix[i][1][0],2)+pow(phiMatrix[i][9]-jetMatrix[i][1][1],2));

R_b1_j1->Fill(dR_b1_j1);
R_b1_j2->Fill(dR_b1_j2);
R_b2_j1->Fill(dR_b2_j1);
R_b2_j2->Fill(dR_b2_j2);

if(dR_b1_j1 <= 0.3 && dR_b2_j2 <=0.3  || dR_b1_j2<=0.3 && dR_b2_j1 <= 0.3 )
ZMassMatched->Fill(ZMassMatrix[i]);
}
}

file.close();
//delete file;

//Here is just plotting the different hisograms
c1->Divide(1,2);
c1->cd(1);
h1->Draw();
c1->cd(2);
h2->Draw();
TCanvas *c2 = new TCanvas("w mass","W mass",1000,600);
c2->Divide(1,2);
c2->cd(1);
WmassCut->SetLabelSize(0.06);
WmassCut->Draw("e");
c2->cd(2);
Zmass->SetLabelSize(0.06);
Zmass->Draw("e");
//cout << " FINALLY = " << goodMatchCount <<endl;
/*
TCanvas *cr = new TCanvas("delta R","delta R",1000,600);
cr->Divide(2,2);
cr->cd(1);
R_b1_j1->Draw("e");
cr->cd(2);
R_b1_j2->Draw("e");
cr->cd(3);
R_b2_j1->Draw("e");
cr->cd(4);
R_b2_j2->Draw("e");
*/
TCanvas *zCheck = new TCanvas("Z_mass_Matched","Z mass matched to Truth",600,400);
ZMassMatched->Draw("e");

TCanvas *jets = new TCanvas("Jets stuff","Jets of selected events",600,600);
jets->Divide(1,2);
jets->cd(1);
numJetsAfterCuts->Draw();
jets->cd(2);
numbJets->Draw();

/*
TCanvas *w = new TCanvas("W mass", " W mass",600,400);
Wmass->Draw("e");
*/
/*
TCanvas *integrated = new TCanvas("integrated masses"," Fully cut masses",600,800);
integrated->Divide(1,2);
integrated->cd(1);
integratedW->Draw("e");
integrated->cd(2);
integratedZ->Draw("e");
cout << endl << " Number of events between 55 and 100 Gev = " << integrator << endl;
*/
}
