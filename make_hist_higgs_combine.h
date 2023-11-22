#include <TF1.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <TGraphAsymmErrors.h>
#include <TVector3.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TMath.h>
#include <fstream>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <THStack.h>
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/Point3D.h"

#define LUMI2016APV 19.52
#define LUMI2016 16.81
#define LUMI2017 41.48
#define LUMI2018 59.83
//#define LUMI2018 137.65
#define mZ 91

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > LV;

struct LepInfo{
        float pt;
        float eta;
        float phi;
	float tight;
	int pdgid;
	float mt;
};

struct FatJetInfo{
	float sdmass;
        float pt;
        float eta;
        float phi;
        float tau21;
        float tau2;
        float tau1;
        float deepMDW;
        float deepW;
	float deepMDZ;
        float deepZ;
        float deepT;
        float deepMDbb;
	float WPid;
};

struct JetInfo{
        float pt;
        float eta;
        float phi;
        bool passbloose;
        bool passbmedium;
        bool overlapwithfatjet;
};

struct Event{
	float scale;
        vector<LepInfo> lep;
        vector<FatJetInfo> fatjet;
        vector<JetInfo> jet;
	int run;
	unsigned long long event;
        float mll;
        float ht;
        float st;
	float stjesup;
	float stjesdn;
	float stjerup;
	float stjerdn;
	float deltaRll;
        float deltaphill;
	float njet_overlapremoved;
	float nbloose;
	float nbmedium;
	float met;
        bool ist;
        bool isW;
        bool issingleWq;
};


class HistCollection{
        protected:
	TH1F* htotal;
	TH1F* hpuWeightup;
	TH1F* hpuWeightdn;
	TH1F* hprefireWeightup;	
	TH1F* hprefireWeightdn;
	TH1F* htriggerWeightup;
	TH1F* htriggerWeightdn;
	TH1F* hlepSFelup;
	TH1F* hlepSFeldn;
        TH1F* hlepSFmuup;
	TH1F* hlepSFmudn;
	TH1F* hbtagSFup;
	TH1F* hbtagSFdn;
	TH1F* hjesup;
	TH1F* hjesdn;
	TH1F* hjerup;
	TH1F* hjerdn;
	TH1F* hqcdscale[6];
	TH1F* hpdf[100];

        TH1F* htotalEFT[216];
	TH1F* hpuWeightupEFT[216];
        TH1F* hpuWeightdnEFT[216];
        TH1F* hprefireWeightupEFT[216];
        TH1F* hprefireWeightdnEFT[216];
        TH1F* htriggerWeightupEFT[216];
        TH1F* htriggerWeightdnEFT[216];
        TH1F* hlepSFelupEFT[216];
        TH1F* hlepSFeldnEFT[216];
        TH1F* hlepSFmuupEFT[216];
        TH1F* hlepSFmudnEFT[216];
        TH1F* hbtagSFupEFT[216];
        TH1F* hbtagSFdnEFT[216];
        TH1F* hjesupEFT[216];
        TH1F* hjesdnEFT[216];
        TH1F* hjerupEFT[216];
        TH1F* hjerdnEFT[216];


        vector<unsigned long long> checkDuplicates;
/*	
	struct{
		vector<unsigned long long> checkrun;
		vector<int> checkevt;
		int size;
	}checkDuplicates;
*/	
	virtual bool Cut(const Event &evt);
	virtual bool LepID(const LepInfo &lep);
	virtual bool FatJetID(const FatJetInfo &fatjet);
	virtual bool JetID(const JetInfo &jet);
        public:
        HistCollection();
        ~HistCollection();
        bool Write(const char* path);
        bool Fill(const float lumi,const float cross_section,const char* path);
};

HistCollection::HistCollection(){
	checkDuplicates.clear();
//	checkDuplicates.checkrun.clear();
//	checkDuplicates.checkevt.clear();
//	checkDuplicates.size=0;
	float x[4]={0,800,1400,1500};
	htotal=		new TH1F("htotal","",3,x);
	hpuWeightup=	new TH1F("pu_weight_Up","",3,x);
	hpuWeightdn=	new TH1F("pu_weight_Down","",3,x);
	hprefireWeightup=new TH1F("prefire_weight_Up","",3,x);
	hprefireWeightdn=new TH1F("prefire_weight_Down","",3,x);
	htriggerWeightup=new TH1F("trigger_weight_Up","",3,x);
	htriggerWeightdn=new TH1F("trigger_weight_Down","",3,x);
	hlepSFelup=	new TH1F("lep_sf_el_tight_Up","",3,x);
	hlepSFeldn=	new TH1F("lep_sf_el_tight_Down","",3,x);
	hlepSFmuup=	new TH1F("lep_sf_mu_tight_Up","",3,x);
	hlepSFmudn=	new TH1F("lep_sf_mu_tight_Down","",3,x);
	hbtagSFup=	new TH1F("btag_sf_medium_Up","",3,x);
       	hbtagSFdn=	new TH1F("btag_sf_medium_Down","",3,x);
	hjesup=		new TH1F("jes_Up","",3,x);
        hjesdn=		new TH1F("jes_Down","",3,x);
        hjerup=		new TH1F("jer_Up","",3,x);
        hjerdn=		new TH1F("jer_Down","",3,x);
        TString wcpoint_large[18]={"_m100p0","_m50p0","_m10p0","_m7p0","_m5p0","_m3p0","_m1p0","_m0p5","_m0p1","_0p1","_0p5","_1p0","_3p0","_5p0","_7p0","_10p0","_50p0","_100p0"};
        TString wcpoint_small[18]={"_m10p0","_m5p0","_m1p0","_m0p7","_m0p5","_m0p3","_m0p1","_m0p05","_m0p01","_0p01","_0p05","_0p1","_0p3","_0p5","_0p7","_1p0","_5p0","_10p0"};
        TString wc[12]={"cW","cHbox","cHDD","cHW","cHB","cHWB","cHl3","cHq1","cHq3","cHu","cHd","cll1"};
	for(unsigned int i=0;i<6;i++){
		hqcdscale[i]=new TH1F(TString::Format("qcd_%d",i),"",3,x);
	}
	for(unsigned int i=0;i<100;i++){
		hpdf[i]=new TH1F(TString::Format("pdf_%d",i),"",3,x);
	}
	for(unsigned int i=0;i<12;i++){
		TString *wcp=wcpoint_large;
		if(i==0 || i==3){
			wcp=wcpoint_small;
		}
		for(unsigned int j=0;j<18;j++){
			htotalEFT[i*18+j]=		new TH1F("htotal_"+wc[i]+wcp[j],"",3,x);
			hpuWeightupEFT[i*18+j]=		new TH1F(wc[i]+wcp[j]+"_pu_weight_Up","",3,x);
			hpuWeightdnEFT[i*18+j]=		new TH1F(wc[i]+wcp[j]+"_pu_weight_Down","",3,x);
			hprefireWeightupEFT[i*18+j]=	new TH1F(wc[i]+wcp[j]+"_prefire_weight_Up","",3,x);
			hprefireWeightdnEFT[i*18+j]=	new TH1F(wc[i]+wcp[j]+"_prefire_weight_Down","",3,x);
			htriggerWeightupEFT[i*18+j]=	new TH1F(wc[i]+wcp[j]+"_trigger_weight_Up","",3,x);
			htriggerWeightdnEFT[i*18+j]=	new TH1F(wc[i]+wcp[j]+"_trigger_weight_Down","",3,x);
			hlepSFelupEFT[i*18+j]=		new TH1F(wc[i]+wcp[j]+"_lep_sf_el_tight_Up","",3,x);
			hlepSFeldnEFT[i*18+j]=		new TH1F(wc[i]+wcp[j]+"_lep_sf_el_tight_Down","",3,x);
			hlepSFmuupEFT[i*18+j]=		new TH1F(wc[i]+wcp[j]+"_lep_sf_mu_tight_Up","",3,x);
			hlepSFmudnEFT[i*18+j]=		new TH1F(wc[i]+wcp[j]+"_lep_sf_mu_tight_Down","",3,x);
			hbtagSFupEFT[i*18+j]=		new TH1F(wc[i]+wcp[j]+"_btag_sf_medium_Up","",3,x);
			hbtagSFdnEFT[i*18+j]=		new TH1F(wc[i]+wcp[j]+"_btag_sf_medium_Down","",3,x);
			hjesupEFT[i*18+j]=		new TH1F(wc[i]+wcp[j]+"_jes_Up","",3,x);
			hjesdnEFT[i*18+j]=		new TH1F(wc[i]+wcp[j]+"_jes_Down","",3,x);
			hjerupEFT[i*18+j]=		new TH1F(wc[i]+wcp[j]+"_jer_Up","",3,x);
			hjerdnEFT[i*18+j]=		new TH1F(wc[i]+wcp[j]+"_jer_Down","",3,x);
		}
	}

}

HistCollection::~HistCollection(){
	delete htotal;
	delete hpuWeightup;
        delete hpuWeightdn;
        delete hprefireWeightup;
        delete hprefireWeightdn;
        delete htriggerWeightup;
        delete htriggerWeightdn;
        delete hlepSFelup;
        delete hlepSFeldn;
        delete hlepSFmuup;
        delete hlepSFmudn;
        delete hbtagSFup;
        delete hbtagSFdn;

	delete hjesup;
        delete hjesdn;
        delete hjerup;
        delete hjerdn;

	for(unsigned int i=0;i<6;i++){
		delete hqcdscale[i];
	}
	for(unsigned int i=0;i<100;i++){
		delete hpdf[i];
	}

	for(unsigned int i=0;i<216;i++){
		delete htotalEFT[i];
		delete hpuWeightupEFT[i];
		delete hpuWeightdnEFT[i];
		delete hprefireWeightupEFT[i];
		delete hprefireWeightdnEFT[i];
		delete htriggerWeightupEFT[i];
		delete htriggerWeightdnEFT[i];
		delete hlepSFelupEFT[i];
		delete hlepSFeldnEFT[i];
		delete hlepSFmuupEFT[i];
		delete hlepSFmudnEFT[i];
		delete hbtagSFupEFT[i];
		delete hbtagSFdnEFT[i];

		delete hjesupEFT[i];
		delete hjesdnEFT[i];
		delete hjerupEFT[i];
		delete hjerdnEFT[i];
	}

}

bool HistCollection::Cut(const Event &evt){
        return true;
}

bool HistCollection::LepID(const LepInfo &lep){
	return true;
}

bool HistCollection::FatJetID(const FatJetInfo &fatjet){
	return true;
}

bool HistCollection::JetID(const JetInfo &jet){
	return true;
}

bool HistCollection::Fill(const float lumi,const float cross_section,const char* path){
    TChain *chain=new TChain("t");
	Event evt;
    vector<LV> *lepp4=0,*jetp4=0,*fatjetp4=0,*genp4=0;
	LV *met=0,*metjesup=0,*metjesdn=0,*metjerup=0,*metjerdn=0;
	vector<float> *fatjetpt=0;
	vector<bool> *jetbloose=0,*jetbmedium=0;
	vector<float> *fatjetsdm=0,*tau21=0,*tau2=0,*tau1=0,*deepMDW=0,*deepW=0,*deepMDZ=0,*deepZ=0,*deepT=0,*deepMDbb=0,*EFT=0,*fatjetptjesup=0,*fatjetptjesdn=0,*fatjetptjerup=0,*fatjetptjerdn=0,*jetptjesup=0,*jetptjesdn=0,*jetptjerup=0,*jetptjerdn=0,*qcd=0,*pdf=0;
	vector<int> *leptight=0,*WPid=0,*WPida=0,*pdgid=0,*genpdgid=0,*genmotheridx=0,*genmotherid=0;
	float puWeight,puWeightup,puWeightdn;
	float prefireWeight,prefireWeightup,prefireWeightdn;
	float triggerWeight,triggerWeightup,triggerWeightdn;
	float lepSF,lepSFelup,lepSFeldn,lepSFmuup,lepSFmudn;
	float btagSF,btagSFup,btagSFdn; 
	float genWeight;
	int nbmedium,nbloose,run,isData;
	bool pass_duplicate_removal;
	unsigned long long event;
	chain->SetBranchAddress("Common_lep_tight",		&leptight);
    chain->SetBranchAddress("Common_lep_p4",                &lepp4);
	chain->SetBranchAddress("Common_lep_pdgid",		&pdgid);
    chain->SetBranchAddress("Common_jet_p4",                &jetp4);
	chain->SetBranchAddress("Common_jet_passBloose",	&jetbloose);
	chain->SetBranchAddress("Common_jet_passBmedium",       &jetbmedium);
	chain->SetBranchAddress("Common_nb_loose",		&nbloose);
	chain->SetBranchAddress("Common_nb_medium",             &nbmedium);
    chain->SetBranchAddress("Common_fatjet_p4",             &fatjetp4);
    chain->SetBranchAddress("Common_fatjet_msoftdrop",      &fatjetsdm);
    chain->SetBranchAddress("Common_fatjet_tau21",          &tau21);
    chain->SetBranchAddress("Common_fatjet_tau2",           &tau2);
    chain->SetBranchAddress("Common_fatjet_tau1",           &tau1);
	chain->SetBranchAddress("Common_fatjet_particleNetMD_W",       &deepMDW);
    chain->SetBranchAddress("Common_fatjet_deep_W",         &deepW);
    chain->SetBranchAddress("Common_fatjet_deepMD_Z",       &deepMDZ);
	chain->SetBranchAddress("Common_fatjet_deep_Z",         &deepZ);
    chain->SetBranchAddress("Common_fatjet_deep_T",         &deepT);
    chain->SetBranchAddress("Common_fatjet_deepMD_bb",      &deepMDbb);
	chain->SetBranchAddress("Common_fatjet_WP_MD",		&WPid);
    chain->SetBranchAddress("Common_met_p4",		&met);
	chain->SetBranchAddress("Common_run",			&run);
	chain->SetBranchAddress("Common_evt",			&event);
	chain->SetBranchAddress("Common_isData",		&isData);
	chain->SetBranchAddress("Common_LHEWeight_mg_reweighting", &EFT);
    chain->SetBranchAddress("Common_genWeight",             &genWeight);
    chain->SetBranchAddress("SS2jet_raw_gen_mother_idx",	&genmotheridx);
    chain->SetBranchAddress("SS2jet_raw_gen_mother_id",	&genmotherid);
	chain->SetBranchAddress("SS2jet_raw_gen_pdgid",		&genpdgid);
	chain->SetBranchAddress("SS2jet_raw_gen_p4s",		&genp4);

	//scale factors

	chain->SetBranchAddress("Common_event_puWeight",	&puWeight);
	chain->SetBranchAddress("Common_event_puWeightup",	&puWeightup);
	chain->SetBranchAddress("Common_event_puWeightdn",	&puWeightdn);
	chain->SetBranchAddress("Common_event_prefireWeight",	&prefireWeight);
	chain->SetBranchAddress("Common_event_prefireWeightup",	&prefireWeightup);
	chain->SetBranchAddress("Common_event_prefireWeightdn",	&prefireWeightdn);
	chain->SetBranchAddress("Common_event_triggerWeight",	&triggerWeight);
	chain->SetBranchAddress("Common_event_triggerWeightup",	&triggerWeightup);
	chain->SetBranchAddress("Common_event_triggerWeightdn",	&triggerWeightdn);
	chain->SetBranchAddress("Common_event_lepSFTight",	&lepSF);
	chain->SetBranchAddress("Common_event_lepSFelupTight",	&lepSFelup);
	chain->SetBranchAddress("Common_event_lepSFeldnTight",	&lepSFeldn);
	chain->SetBranchAddress("Common_event_lepSFmuupTight",	&lepSFmuup);
    chain->SetBranchAddress("Common_event_lepSFmudnTight",	&lepSFmudn);
	chain->SetBranchAddress("Common_event_mediumBtagSF",	&btagSF);
	chain->SetBranchAddress("Common_event_mediumBtagSFup",	&btagSFup);
	chain->SetBranchAddress("Common_event_mediumBtagSFdn",	&btagSFdn);
	chain->SetBranchAddress("Common_event_QCDScale",	&qcd);
	chain->SetBranchAddress("Common_event_PDF",		&pdf);

	chain->SetBranchAddress("Common_met_p4_jesup",&metjesup);	
	chain->SetBranchAddress("Common_met_p4_jesdn",&metjesdn);
    chain->SetBranchAddress("Common_met_p4_jerup",&metjerup);
    chain->SetBranchAddress("Common_met_p4_jerdn",&metjerdn);
	chain->SetBranchAddress("Common_fatjet_pt_jesup",&fatjetptjesup);
    chain->SetBranchAddress("Common_fatjet_pt_jesdn",&fatjetptjesdn);
    chain->SetBranchAddress("Common_fatjet_pt_jerup",&fatjetptjerup);
    chain->SetBranchAddress("Common_fatjet_pt_jerdn",&fatjetptjerdn);
	chain->SetBranchAddress("Common_jet_pt_jesup",&jetptjesup);
    chain->SetBranchAddress("Common_jet_pt_jesdn",&jetptjesdn);
    chain->SetBranchAddress("Common_jet_pt_jerup",&jetptjerup);
	chain->SetBranchAddress("Common_jet_pt_jerdn",&jetptjerdn);



	chain->Add(path);
	TFile *file=TFile::Open(path);
	float scale;
	if(!file->IsOpen()) return false;
	scale = cross_section * 1000 * lumi / ((TH1F*) file->Get("Wgt__h_nevents"))->Integral();
    for(unsigned int i=0;i<chain->GetEntries();i++){
        chain->GetEntry(i);
		if(EFT->size()!=0)      scale = cross_section * 1000 * lumi / ((TH1F*) file->Get("Root__h_Common_LHEWeight_mg_reweighting_times_genWeight"))->GetBinContent(1) * EFT->at(0);

		if(isData) evt.scale=1;
		else evt.scale =  ( (genWeight>0) - (genWeight<0) ) * scale * puWeight * prefireWeight * triggerWeight * lepSF * btagSF;
		//build the event
		evt.lep.clear();
		evt.fatjet.clear();
		evt.jet.clear();
		evt.st=0;
		evt.stjesup=0;
        evt.stjesdn=0;
        evt.stjerup=0;
        evt.stjerdn=0;

		evt.nbloose=0;
		evt.nbmedium=0;
		evt.njet_overlapremoved=0;
		evt.run=run;
		evt.event=event;

		for(unsigned int ilep=0;ilep<lepp4->size();ilep++){
			LepInfo temp;
			temp.pt=lepp4->at(ilep).pt();
			temp.eta=lepp4->at(ilep).eta();
			temp.phi=lepp4->at(ilep).phi();
			temp.tight=leptight->at(ilep);
			temp.pdgid=pdgid->at(ilep);
			float phi1 = lepp4->at(ilep).Phi();
		    float phi2 = met->Phi();
			float Et1  = lepp4->at(ilep).Et();
			float Et2  = met->Et();
			temp.mt=sqrt(2*Et1*Et2*(1.0 - cos(phi1-phi2)));
			if(!LepID(temp)) continue;
			evt.lep.push_back(temp);
			evt.st+=lepp4->at(ilep).pt();
			evt.stjesup+=lepp4->at(ilep).pt();
            evt.stjesdn+=lepp4->at(ilep).pt();
            evt.stjerup+=lepp4->at(ilep).pt();
            evt.stjerdn+=lepp4->at(ilep).pt();

		}

		evt.deltaRll=ROOT::Math::VectorUtil::DeltaR(lepp4->at(0),lepp4->at(1));
		evt.mll=(lepp4->at(0)+lepp4->at(1)).M();

		for(unsigned int ifatjet=0;ifatjet<fatjetp4->size();ifatjet++){
			FatJetInfo temp;
			temp.sdmass=fatjetsdm->at(ifatjet);
			temp.pt=fatjetp4->at(ifatjet).pt();
			temp.eta=fatjetp4->at(ifatjet).eta();
			temp.phi=fatjetp4->at(ifatjet).phi();
			temp.tau21=tau21->at(ifatjet);
			temp.tau2=tau2->at(ifatjet);
			temp.tau1=tau1->at(ifatjet);
			temp.deepMDW=deepMDW->at(ifatjet);
			temp.deepW=deepW->at(ifatjet);
			temp.deepMDZ=deepZ->at(ifatjet);
			temp.deepZ=deepZ->at(ifatjet);
			temp.deepT=deepT->at(ifatjet);
			temp.deepMDbb=deepMDbb->at(ifatjet);
			temp.WPid=WPid->at(ifatjet);
			if(!FatJetID(temp)) continue;	
			
			evt.st+=fatjetp4->at(ifatjet).pt();
			if(fatjetptjesup->size()>0){
				evt.stjesup+=fatjetptjesup->at(ifatjet);
        	                evt.stjesdn+=fatjetptjesdn->at(ifatjet);
	                        evt.stjerup+=fatjetptjerup->at(ifatjet);
                        	evt.stjerdn+=fatjetptjerdn->at(ifatjet);
			}
			evt.fatjet.push_back(temp);

		}
		for(unsigned int ijet=0;ijet<jetp4->size();ijet++){
			JetInfo temp;
			temp.pt=jetp4->at(ijet).pt();
			temp.eta=jetp4->at(ijet).eta();
			temp.phi=jetp4->at(ijet).phi();
			temp.passbloose=jetbloose->at(ijet);
			temp.passbmedium=jetbmedium->at(ijet);
			bool overlap=false;
			for(unsigned int ifatjet=0;ifatjet<fatjetp4->size();ifatjet++){
				if(ROOT::Math::VectorUtil::DeltaR(jetp4->at(ijet),fatjetp4->at(ifatjet))<0.8){
					overlap=true;
					break;
				}				
			}
			temp.overlapwithfatjet=overlap;
			if(!JetID(temp)) continue;
			evt.jet.push_back(temp);
			if(!overlap){
				evt.st+=jetp4->at(ijet).pt();
				if(jetptjesup->size()>0){
					evt.stjesup+=jetptjesup->at(ijet);
					evt.stjesdn+=jetptjesdn->at(ijet);
					evt.stjerup+=jetptjerup->at(ijet);
					evt.stjerdn+=jetptjerdn->at(ijet);
				}
				evt.ht+=jetp4->at(ijet).pt();
				evt.njet_overlapremoved++;
			}
		}
		evt.nbmedium=nbmedium;
		evt.nbloose=nbloose;
		evt.met=met->pt();
		evt.st+=met->pt();
		evt.stjesup+=metjesup->pt();
		evt.stjesdn+=metjesdn->pt();
		evt.stjerup+=metjerup->pt();
		evt.stjerdn+=metjerdn->pt();
                
		evt.st=evt.st>1500?1499.5:evt.st;
        	evt.stjesup=evt.stjesup>1500?1499.5:evt.stjesup;
		evt.stjesdn=evt.stjesdn>1500?1499.5:evt.stjesdn;
        	evt.stjerup=evt.stjerup>1500?1499.5:evt.stjerup;
        	evt.stjerdn=evt.stjerdn>1500?1499.5:evt.stjerdn;

		//build the event
		//fill in histograms
		if(Cut(evt)){//addition cut
			htotal->Fill(evt.st,				evt.scale);
			hpuWeightup->Fill(evt.st,			evt.scale/puWeight*puWeightup);
			hpuWeightdn->Fill(evt.st,			evt.scale/puWeight*puWeightdn);
			hprefireWeightup->Fill(evt.st,			evt.scale/prefireWeight*prefireWeightup);
			hprefireWeightdn->Fill(evt.st,			evt.scale/prefireWeight*prefireWeightdn);
			htriggerWeightup->Fill(evt.st,			evt.scale/triggerWeight*triggerWeightup);
			htriggerWeightdn->Fill(evt.st,			evt.scale/triggerWeight*triggerWeightdn);
			hlepSFelup->Fill(evt.st,			evt.scale/lepSF*lepSFelup);
			hlepSFeldn->Fill(evt.st,			evt.scale/lepSF*lepSFeldn);
			hlepSFmuup->Fill(evt.st,			evt.scale/lepSF*lepSFmuup);
			hlepSFmudn->Fill(evt.st,			evt.scale/lepSF*lepSFmudn);
			hbtagSFup->Fill(evt.st,				evt.scale/btagSF*btagSFup);
			hbtagSFdn->Fill(evt.st,				evt.scale/btagSF*btagSFdn);
			hjesup->Fill(evt.stjesup,			evt.scale);
			hjesdn->Fill(evt.stjesdn,                      	evt.scale);
			hjerup->Fill(evt.stjerup,                       evt.scale);
			hjerdn->Fill(evt.stjerdn,                       evt.scale);
			if(qcd->size()==9){
				hqcdscale[0]->Fill(evt.st,		evt.scale*qcd->at(0));
				hqcdscale[1]->Fill(evt.st,		evt.scale*qcd->at(1));
				hqcdscale[2]->Fill(evt.st,		evt.scale*qcd->at(3));
				hqcdscale[3]->Fill(evt.st,		evt.scale*qcd->at(5));	
				hqcdscale[4]->Fill(evt.st,		evt.scale*qcd->at(7));
				hqcdscale[5]->Fill(evt.st,		evt.scale*qcd->at(8));
			}
			else{
				hqcdscale[0]->Fill(evt.st,		evt.scale*qcd->at(0));
				hqcdscale[1]->Fill(evt.st,		evt.scale*qcd->at(1));
				hqcdscale[2]->Fill(evt.st,		evt.scale*qcd->at(3));
				hqcdscale[3]->Fill(evt.st,		evt.scale*qcd->at(4));	
				hqcdscale[4]->Fill(evt.st,		evt.scale*qcd->at(6));
				hqcdscale[5]->Fill(evt.st,		evt.scale*qcd->at(7));
			}
			for(unsigned int i=0;i<100;i++){
				hpdf[i]->Fill(evt.st,			evt.scale*pdf->at(i+1));
			}

			if(EFT->size()!=0){
				for(unsigned i=0;i<216;i++){
					float scale=evt.scale/EFT->at(0)*EFT->at(i+1);
					htotalEFT[i]->Fill(evt.st,				scale);
					hpuWeightupEFT[i]->Fill(evt.st,				scale/puWeight*puWeightup);
					hpuWeightdnEFT[i]->Fill(evt.st,				scale/puWeight*puWeightdn);
					hprefireWeightupEFT[i]->Fill(evt.st,			scale/prefireWeight*prefireWeightup);
					hprefireWeightdnEFT[i]->Fill(evt.st,			scale/prefireWeight*prefireWeightdn);
					htriggerWeightupEFT[i]->Fill(evt.st,			scale/triggerWeight*triggerWeightup);
					htriggerWeightdnEFT[i]->Fill(evt.st,			scale/triggerWeight*triggerWeightdn);
					hlepSFelupEFT[i]->Fill(evt.st,				scale/lepSF*lepSFelup);
					hlepSFeldnEFT[i]->Fill(evt.st,				scale/lepSF*lepSFeldn);
					hlepSFmuupEFT[i]->Fill(evt.st,				scale/lepSF*lepSFmuup);
					hlepSFmudnEFT[i]->Fill(evt.st,				scale/lepSF*lepSFmudn);
					hbtagSFupEFT[i]->Fill(evt.st,				scale/btagSF*btagSFup);
					hbtagSFdnEFT[i]->Fill(evt.st,				scale/btagSF*btagSFdn);
					hjesupEFT[i]->Fill(evt.stjesup,				scale);
					hjesdnEFT[i]->Fill(evt.stjesdn,                       	scale);
					hjerupEFT[i]->Fill(evt.stjerup,                       	scale);
					hjerdnEFT[i]->Fill(evt.stjerdn,                       	scale);
				}
			}			
			//fill in histograms
		}
	}
	return true;
}

bool HistCollection::Write(const char* path){
	TFile *file=TFile::Open(path,"RECREATE");
	htotal->Write();
	hpuWeightup->Write();
	hpuWeightdn->Write();
	hprefireWeightup->Write();
	hprefireWeightdn->Write();
	htriggerWeightup->Write();
	htriggerWeightdn->Write();
	hlepSFelup->Write();
	hlepSFeldn->Write();
	hlepSFmuup->Write();
	hlepSFmudn->Write();
	hbtagSFup->Write();
	hbtagSFdn->Write();

	hjesup->Write();
	hjesdn->Write();
	hjerup->Write();
	hjerdn->Write();

	for(unsigned int i=0;i<6;i++){
		hqcdscale[i]->Write();
	}

	for(unsigned int i=0;i<100;i++){
		hpdf[i]->Write();
	}

	for(unsigned i=0;i<216;i++){
		htotalEFT[i]->Write();
		hpuWeightupEFT[i]->Write();
		hpuWeightdnEFT[i]->Write();
		hprefireWeightupEFT[i]->Write();
		hprefireWeightdnEFT[i]->Write();
		htriggerWeightupEFT[i]->Write();
		htriggerWeightdnEFT[i]->Write();
		hlepSFelupEFT[i]->Write();
		hlepSFeldnEFT[i]->Write();
		hlepSFmuupEFT[i]->Write();
		hlepSFmudnEFT[i]->Write();
		hbtagSFupEFT[i]->Write();
		hbtagSFdnEFT[i]->Write();

		hjesupEFT[i]->Write();
		hjesdnEFT[i]->Write();
		hjerupEFT[i]->Write();
		hjerdnEFT[i]->Write();
	}


	file->Close();
	return true;
}

class SignalRegion:public HistCollection{
	protected:
		bool LepID(const LepInfo &lep){
			return (lep.tight==1);
		}

                bool FatJetID(const FatJetInfo &fatjet){
//			return true;
		       	return (fatjet.WPid>=2 && fatjet.sdmass<=105. && fatjet.sdmass>=65. && fatjet.pt>200.);
//			return (fatjet.deepMDW>0.82 && fatjet.sdmass<=105. && fatjet.sdmass>=65. && fatjet.pt>200.);
                }
                bool Cut(const Event &evt){
                        if(!(evt.lep.size()==2)) return false;
                        if(!(evt.fatjet.size()>=1)) return false;
			//if(!(evt.lep.at(0).pdgid * evt.lep.at(1).pdgid==143)) return false;
                        if(!(evt.nbmedium==0)) return false;
                        if(evt.lep.at(0).pdgid * evt.lep.at(1).pdgid==121 && abs(evt.mll-mZ)<20) return false;
                        if(!(evt.lep.at(0).pt>40 && evt.lep.at(1).pt>30)) return false;
			if(!(evt.deltaRll>1.2)) return false;
			//if(!(evt.met>60)) return false;
                        return true;
                }
};



