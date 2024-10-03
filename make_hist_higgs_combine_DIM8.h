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

#include "FakeFatjetScaleFactorProvider_new.h"

#define LUMI2016APV 19.52
#define LUMI2016 16.81
#define LUMI2017 41.48
#define LUMI2018 59.83
//#define LUMI2018 137.65
#define mZ 91

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > LV;

struct LepInfo{
	LV p4;
	int tight;
	int pdgid;
};

struct FatJetInfo{
	LV p4;
	float sdmass;
	float WPid;
};

struct JetInfo{
	LV p4;
	bool overlapwithfatjet;
};

struct Event{
	float scale;
	vector<LepInfo> lep;
	vector<FatJetInfo> fatjet;
	vector<JetInfo> jet;
	int run;
	unsigned long long event;
	float met;
	float mll;
	float st;
	float deltaRll;
	int nbmedium;
	int WFatJet;
};


class HistCollection{
	protected:
	TH1F* hcentral;
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
	TH1F* hfatjetSFup;
	TH1F* hfatjetSFdn;
	TH1F* hjesup[27];
	TH1F* hjesdn[27];
	TH1F* hjerup;
	TH1F* hjerdn;
	TH1F* halphasup;
	TH1F* halphasdn;
	TH1F* hqcdscale[6];
	TH1F* hpdf[100];

	//TH2F* hcentral_vs_EFTCoefficient;
	TH1F* hcentralEFT[252];
	TH1F* hpuWeightupEFT[252];
	TH1F* hpuWeightdnEFT[252];
	TH1F* hprefireWeightupEFT[252];
	TH1F* hprefireWeightdnEFT[252];
	TH1F* htriggerWeightupEFT[252];
	TH1F* htriggerWeightdnEFT[252];
	TH1F* hlepSFelupEFT[252];
	TH1F* hlepSFeldnEFT[252];
	TH1F* hlepSFmuupEFT[252];
	TH1F* hlepSFmudnEFT[252];
	TH1F* hbtagSFupEFT[252];
	TH1F* hbtagSFdnEFT[252];
	TH1F* hfatjetSFupEFT[252];
	TH1F* hfatjetSFdnEFT[252];
	TH1F* hjesupEFT[27][252];
	TH1F* hjesdnEFT[27][252];
	TH1F* hjerupEFT[252];
	TH1F* hjerdnEFT[252];
	TH1F* halphasupEFT[252];
	TH1F* halphasdnEFT[252];
	TH1F* hqcdscaleEFT[6][252];
	TH1F* hpdfEFT[100][252];

	TString jes[27]={
		"AbsoluteStat",
		"AbsoluteScale",
		"AbsoluteMPFBias",
		"Fragmentation",
		"SinglePionECAL",
		"SinglePionHCAL",
		"FlavorQCD",
		"TimePtEta",
		"RelativeJEREC1",
		"RelativeJEREC2",
		"RelativeJERHF",
		"RelativePtBB",
		"RelativePtEC1",
		"RelativePtEC2",
		"RelativePtHF",
		"RelativeBal",
		"RelativeSample",
		"RelativeFSR",
		"RelativeStatFSR",
		"RelativeStatEC",
		"RelativeStatHF",
		"PileUpDataMC",
		"PileUpPtRef",
		"PileUpPtBB",
		"PileUpPtEC1",
		"PileUpPtEC2",
		"PileUpPtHF"
	};
	FakeFatjetScaleFactorProvider fakefatjetsfprovider;
	TString wcpoint_large[12]={"_m3000p0","_m1500p0","_m1000p0","_m800p0","_m400p0","_m200p0","_200p0","_400p0","_800p0","_1000p0","_1500p0","_3000p0"};
	TString wcpoint_medium[12]={"_m30p0","_m15p0","_m10p0","_m8p0","_m4p0","_m2p0","_2p0","_4p0","_8p0","_10p0","_15p0","_30p0"};
	TString wcpoint_small[12]={"_m3p0","_m1p5","_m1p0","_m0p8","_m0p4","_m0p2","_0p2","_0p4","_0p8","_1p0","_1p5","_3p0"};
	int WWW_index[21]={1,13,25,37,49,0 ,0 ,0 ,0 ,61 ,73 ,85 ,97 ,109,121,133,0  ,0  ,0  ,0  ,0  };
	int WVZ_index[21]={1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,0  ,0  };
	int ZZZ_index[21]={1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241};
	TString wc[21]={
		"FS0",
		"FS1",
		"FS2",
		"FM0",
		"FM1",
		"FM2",
		"FM3",
		"FM4",
		"FM5",
		"FM6",
		"FM7",
		"FT0",
		"FT1",
		"FT2",
		"FT3",
		"FT4",
		"FT5",
		"FT6",
		"FT7",
		"FT8",
		"FT9"
	};

	vector<unsigned long long> checkDuplicates;
	virtual bool Cut(const Event &evt);
	virtual bool LepID(const LepInfo &lep);
	virtual bool FatJetID(const FatJetInfo &fatjet);
	virtual bool WFatJetID(const FatJetInfo &fatjet);
	virtual bool JetID(const JetInfo &jet);
	float Variable_of_Interest(const Event &evt);
	public:
	HistCollection();
	~HistCollection();
	bool Write(const char* path);
	bool Fill(const float lumi,const float cross_section,const char* path, bool pdf_replica=false);
};

HistCollection::HistCollection(){
	checkDuplicates.clear();
	double x[4]={0,1000,1600,2000};
	//3,x
	//18,200,2000
	hcentral=		new TH1F("hcentral",""				,3,x);
	hpuWeightup=		new TH1F("pu_weight_Up",""			,3,x);
	hpuWeightdn=		new TH1F("pu_weight_Down",""			,3,x);
	hprefireWeightup=	new TH1F("prefire_weight_Up",""			,3,x);
	hprefireWeightdn=	new TH1F("prefire_weight_Down",""		,3,x);
	htriggerWeightup=	new TH1F("trigger_weight_Up",""			,3,x);
	htriggerWeightdn=	new TH1F("trigger_weight_Down",""		,3,x);
	hlepSFelup=		new TH1F("lep_sf_el_tight_Up",""		,3,x);
	hlepSFeldn=		new TH1F("lep_sf_el_tight_Down",""		,3,x);
	hlepSFmuup=		new TH1F("lep_sf_mu_tight_Up",""		,3,x);
	hlepSFmudn=		new TH1F("lep_sf_mu_tight_Down",""		,3,x);
	hbtagSFup=		new TH1F("btag_sf_medium_Up",""			,3,x);
	hbtagSFdn=		new TH1F("btag_sf_medium_Down",""		,3,x);
	hfatjetSFup=		new TH1F("fatjet_sf_medium_Up",""		,3,x);
	hfatjetSFdn=		new TH1F("fatjet_sf_medium_Down",""		,3,x);
	for(unsigned int i=0;i<27;i++){
		hjesup[i]=	new TH1F("jes"+jes[i]+"_Up",""			,3,x);
		hjesdn[i]=	new TH1F("jes"+jes[i]+"_Down",""		,3,x);
	}
	hjerup=			new TH1F("jer_Up",""				,3,x);
	hjerdn=			new TH1F("jer_Down",""				,3,x);
	halphasup=		new TH1F("alphas_Up",""				,3,x);
	halphasdn=		new TH1F("alphas_Down",""			,3,x);
	//hcentral_vs_EFTCoefficient= new TH2F("central_vs_EFTCoefficient","",3,x,91,-0.5,90.5);
	for(unsigned int i=0;i<6;i++){
		hqcdscale[i]=	new TH1F(TString::Format("qcd_%d",i),""		,3,x);
	}
	for(unsigned int i=0;i<100;i++){
		hpdf[i]=		new TH1F(TString::Format("pdf_%d",i),""	,3,x);
	}

	for(unsigned int i=0;i<21;i++){
		TString *wcp=wcpoint_medium;
		if(i==0 || i==1 || i==2){
			wcp=wcpoint_large;
		}
		if(i==11 || i==12 || i==13){
			wcp=wcpoint_small;
		}
		for(unsigned int j=0;j<12;j++){
			hcentralEFT[i*12+j]=		new TH1F("hcentral_"+wc[i]+wcp[j],""			,3,x);
			hpuWeightupEFT[i*12+j]=		new TH1F(wc[i]+wcp[j]+"_pu_weight_Up",""		,3,x);
			hpuWeightdnEFT[i*12+j]=		new TH1F(wc[i]+wcp[j]+"_pu_weight_Down",""		,3,x);
			hprefireWeightupEFT[i*12+j]=	new TH1F(wc[i]+wcp[j]+"_prefire_weight_Up",""		,3,x);
			hprefireWeightdnEFT[i*12+j]=	new TH1F(wc[i]+wcp[j]+"_prefire_weight_Down",""		,3,x);
			htriggerWeightupEFT[i*12+j]=	new TH1F(wc[i]+wcp[j]+"_trigger_weight_Up",""		,3,x);
			htriggerWeightdnEFT[i*12+j]=	new TH1F(wc[i]+wcp[j]+"_trigger_weight_Down",""		,3,x);
			hlepSFelupEFT[i*12+j]=		new TH1F(wc[i]+wcp[j]+"_lep_sf_el_tight_Up",""		,3,x);
			hlepSFeldnEFT[i*12+j]=		new TH1F(wc[i]+wcp[j]+"_lep_sf_el_tight_Down",""	,3,x);
			hlepSFmuupEFT[i*12+j]=		new TH1F(wc[i]+wcp[j]+"_lep_sf_mu_tight_Up",""		,3,x);
			hlepSFmudnEFT[i*12+j]=		new TH1F(wc[i]+wcp[j]+"_lep_sf_mu_tight_Down",""	,3,x);
			hbtagSFupEFT[i*12+j]=		new TH1F(wc[i]+wcp[j]+"_btag_sf_medium_Up",""		,3,x);
			hbtagSFdnEFT[i*12+j]=		new TH1F(wc[i]+wcp[j]+"_btag_sf_medium_Down",""		,3,x);
			hfatjetSFupEFT[i*12+j]=		new TH1F(wc[i]+wcp[j]+"_fatjet_sf_medium_Up",""		,3,x);
			hfatjetSFdnEFT[i*12+j]=		new TH1F(wc[i]+wcp[j]+"_fatjet_sf_medium_Down",""	,3,x);
			for(unsigned int k=0;k<27;k++){
				hjesupEFT[k][i*12+j]=	new TH1F(wc[i]+wcp[j]+"_jes"+jes[k]+"_Up",""		,3,x);
				hjesdnEFT[k][i*12+j]=	new TH1F(wc[i]+wcp[j]+"_jes"+jes[k]+"_Down",""		,3,x);
			}			
			hjerupEFT[i*12+j]=		new TH1F(wc[i]+wcp[j]+"_jer_Up",""			,3,x);
			hjerdnEFT[i*12+j]=		new TH1F(wc[i]+wcp[j]+"_jer_Down",""			,3,x);
			halphasupEFT[i*12+j]=		new TH1F(wc[i]+wcp[j]+"_alphas_Up",""			,3,x);
			halphasdnEFT[i*12+j]=		new TH1F(wc[i]+wcp[j]+"_alphas_Down",""			,3,x);
			for(unsigned int k=0;k<6;k++){
				hqcdscaleEFT[k][i*12+j]=new TH1F(TString::Format(wc[i]+wcp[j]+"_qcd_%d",k),""	,3,x);
			}
			for(unsigned int k=0;k<100;k++){
				hpdfEFT[k][i*12+j]=	new TH1F(TString::Format(wc[i]+wcp[j]+"_pdf_%d",k),""	,3,x);
			}
		}
	}


}

HistCollection::~HistCollection(){
	delete hcentral;
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
	delete hfatjetSFup;
	delete hfatjetSFdn;
	delete halphasup;
	delete halphasdn;
	//delete hcentral_vs_EFTCoefficient;
	for(unsigned int i=0;i<27;i++){
		delete hjesup[i];
		delete hjesdn[i];
	}
	delete hjerup;
	delete hjerdn;
	for(unsigned int i=0;i<6;i++){
		delete hqcdscale[i];
	}
	for(unsigned int i=0;i<100;i++){
		delete hpdf[i];
	}

	for(unsigned int i=0;i<252;i++){
		delete hcentralEFT[i];
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
		delete hfatjetSFupEFT[i];
		delete hfatjetSFdnEFT[i];
		for(unsigned int j=0;j<27;j++){
			delete hjesupEFT[j][i];
			delete hjesdnEFT[j][i];
		}
		delete hjerupEFT[i];
		delete hjerdnEFT[i];
		delete halphasupEFT[i];
		delete halphasdnEFT[i];
		for(unsigned int j=0;j<6;j++){
			delete hqcdscaleEFT[j][i];
		}
		for(unsigned int j=0;j<100;j++){
			delete hpdfEFT[j][i];
		}
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

bool HistCollection::WFatJetID(const FatJetInfo &fatjet){
	return true;
}

bool HistCollection::JetID(const JetInfo &jet){
	return true;
}

float HistCollection::Variable_of_Interest(const Event &evt){
	//For SS+1fatjet channel, varialbe of interest is stmet
	//overlap with lepton is already handle at looper level
	float rtn=0;
	for(unsigned int ilep=0;ilep<evt.lep.size();ilep++){
		rtn+= evt.lep.at(ilep).p4.pt();
	}
	for(unsigned int ifatjet=0;ifatjet<evt.fatjet.size();ifatjet++){
		rtn+= evt.fatjet.at(ifatjet).p4.pt();
	}
	for(unsigned int ijet=0;ijet<evt.jet.size();ijet++){
		bool overlap=false;
		for(unsigned int ifatjet=0;ifatjet<evt.fatjet.size();ifatjet++){
			if(ROOT::Math::VectorUtil::DeltaR(evt.jet.at(ijet).p4,evt.fatjet.at(ifatjet).p4)<0.8){
				overlap=true;
				break;
			}
		}
		if(overlap) continue;
		rtn+= evt.jet.at(ijet).p4.pt();
	}
	rtn+= evt.met;
	return rtn;
}
bool HistCollection::Fill(const float lumi,const float cross_section,const char* path, bool pdf_replica=false){
	TChain *chain=new TChain("t");
	Event evt;
	Event evtvar[56];
	vector<LV> *lepp4=0,*jetp4=0,*fatjetp4=0,*genp4=0;
	vector<LV> *jetp4var[56],*fatjetp4var[56];
	vector<float> *fatjetsdmvar[56]; 
	LV *met=0;
	LV *metp4var[56];
	vector<float> *fatjetpt=0;
	vector<bool> *jetbmedium=0;
	vector<float> *fatjetsdm=0,*EFT=0,*qcd=0,*pdf=0,*EFTCoef=0;
	vector<int> *leptight=0,*WPid=0,*WPida=0,*pdgid=0,*genpdgid=0,*genmotheridx=0,*genmotherid=0;
	float puWeight,puWeightup,puWeightdn;
	float prefireWeight,prefireWeightup,prefireWeightdn;
	float triggerWeight,triggerWeightup,triggerWeightdn;
	float lepSF,lepSFelup,lepSFeldn,lepSFmuup,lepSFmudn;
	float btagSF,btagSFup,btagSFdn; 
	vector<float> *fjSF=0,*fjSFup=0,*fjSFdn=0;
	float genWeight;
	int nbmedium,run,isData;
	int nbmediumvar[56];
	bool pass_duplicate_removal;
	unsigned long long event;
	chain->SetBranchAddress("Common_lep_tight",			&leptight);
	chain->SetBranchAddress("Common_lep_p4",			&lepp4);
	chain->SetBranchAddress("Common_lep_pdgid",			&pdgid);
	chain->SetBranchAddress("Common_jet_p4",			&jetp4);
	chain->SetBranchAddress("Common_nb_medium",			&nbmedium);
	chain->SetBranchAddress("Common_fatjet_p4",			&fatjetp4);
	chain->SetBranchAddress("Common_fatjet_msoftdrop",		&fatjetsdm);
	chain->SetBranchAddress("Common_fatjet_WP_MD",			&WPid);
	chain->SetBranchAddress("Common_met_p4",			&met);
	chain->SetBranchAddress("Common_run",				&run);
	chain->SetBranchAddress("Common_evt",				&event);
	chain->SetBranchAddress("Common_isData",			&isData);
	chain->SetBranchAddress("Common_LHEWeight_mg_reweighting",	&EFT);
	chain->SetBranchAddress("Common_genWeight",			&genWeight);
	//chain->SetBranchAddress("Common_EFTfitCoefficient",		&EFTCoef);
	//jes and jer
	for(unsigned int i=0;i<27;i++){
		jetp4var[i]=0;jetp4var[i+28]=0;fatjetp4var[i]=0;fatjetp4var[i+28]=0;metp4var[i]=0;metp4var[i+28]=0;fatjetsdmvar[i]=0;fatjetsdmvar[i+28]=0;
		chain->SetBranchAddress("Common_jet_p4_jes"+jes[i]+"up",	&jetp4var[i]);
		chain->SetBranchAddress("Common_fatjet_p4_jes"+jes[i]+"up",	&fatjetp4var[i]);
		chain->SetBranchAddress("Common_fatjet_msoftdrop_jes"+jes[i]+"up", &fatjetsdmvar[i]);
		chain->SetBranchAddress("Common_met_p4_jes"+jes[i]+"up",	&metp4var[i]);
		chain->SetBranchAddress("Common_nb_medium_jes"+jes[i]+"up",	&nbmediumvar[i]);

		chain->SetBranchAddress("Common_jet_p4_jes"+jes[i]+"dn",	&jetp4var[i+28]);
		chain->SetBranchAddress("Common_fatjet_p4_jes"+jes[i]+"dn",	&fatjetp4var[i+28]);
		chain->SetBranchAddress("Common_fatjet_msoftdrop_jes"+jes[i]+"dn", &fatjetsdmvar[i+28]);
		chain->SetBranchAddress("Common_met_p4_jes"+jes[i]+"dn",	&metp4var[i+28]);	
		chain->SetBranchAddress("Common_nb_medium_jes"+jes[i]+"dn",	&nbmediumvar[i+28]);
	}
	jetp4var[27]=0;jetp4var[55]=0;fatjetp4var[27]=0;fatjetp4var[55]=0;metp4var[27]=0;metp4var[55]=0;fatjetsdmvar[27]=0;fatjetsdmvar[55]=0;
	chain->SetBranchAddress("Common_jet_p4_jerup",		&jetp4var[27]);
	chain->SetBranchAddress("Common_fatjet_p4_jerup",	&fatjetp4var[27]);
	chain->SetBranchAddress("Common_fatjet_msoftdrop_jerup", &fatjetsdmvar[27]);
	chain->SetBranchAddress("Common_met_p4_jerup",		&metp4var[27]);
	chain->SetBranchAddress("Common_nb_medium_jerup",	&nbmediumvar[27]);

	chain->SetBranchAddress("Common_jet_p4_jerdn",		&jetp4var[55]);
	chain->SetBranchAddress("Common_fatjet_p4_jerdn",	&fatjetp4var[55]);
	chain->SetBranchAddress("Common_fatjet_msoftdrop_jerdn", &fatjetsdmvar[55]);
	chain->SetBranchAddress("Common_met_p4_jerdn",		&metp4var[55]);	
	chain->SetBranchAddress("Common_nb_medium_jerdn",	&nbmediumvar[55]);

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
	chain->SetBranchAddress("Common_fatjet_MD_SFMedium",	&fjSF);
	chain->SetBranchAddress("Common_fatjet_MD_SFupMedium",	&fjSFup);
	chain->SetBranchAddress("Common_fatjet_MD_SFdnMedium",	&fjSFdn);
	chain->SetBranchAddress("Common_event_QCDScale",	&qcd);
	chain->SetBranchAddress("Common_event_PDF",		&pdf);

	chain->SetBranchAddress("Common_gen_mother_idx",	&genmotheridx);
	chain->SetBranchAddress("Common_gen_mother_id",		&genmotherid);
	chain->SetBranchAddress("Common_gen_pdgid",		&genpdgid);
	chain->SetBranchAddress("Common_gen_p4s",		&genp4);

	chain->Add(path);
	TFile *file=TFile::Open(path);
	float scale;
	if(!file->IsOpen()) return false;
	scale = cross_section * 1000 * lumi / ((TH1D*) file->Get("Wgt__h_nevents"))->Integral();
	for(unsigned int i=0;i<chain->GetEntries();i++){
        chain->GetEntry(i);
		if(EFT->size()!=0)      scale = cross_section * 1000 * lumi / ((TH1F*) file->Get("Root__h_Common_LHEWeight_mg_reweighting_times_genWeight"))->GetBinContent(1) * EFT->at(0);
		int *index;
		if(EFT->size()==253){
			index=ZZZ_index;
		}
		else if(EFT->size()==229){
			index=WVZ_index;
		}
		else if(EFT->size()==145){
			index=WWW_index;
		}
		else{
			return EFT->size();
		}

		//else if (EFTCoef->size()!=0)  scale = cross_section * 1000 * lumi / ((TH1F*) file->Get("Root__h_Common_LHEWeight_mg_reweighting_times_genWeight"))->GetBinContent(1) * EFTCoef->at(0);
		if(isData) evt.scale=1;
		else evt.scale =  ( (genWeight>0) - (genWeight<0) ) * scale * puWeight * prefireWeight * triggerWeight * lepSF * btagSF;
		//build the event
		evt.lep.clear();
		evt.fatjet.clear();
		evt.jet.clear();
		evt.st=0;
		evt.nbmedium=0;
		evt.WFatJet=0;
		evt.run=run;
		evt.event=event;
		float fatjetSF=1,fatjetSFup=1,fatjetSFdn=1;

		for(unsigned int ivar=0;ivar<56;ivar++){
			evtvar[ivar].scale =  ( (genWeight>0) - (genWeight<0) ) * scale * puWeight * prefireWeight * triggerWeight * lepSF * btagSF;
			evtvar[ivar].lep.clear();
			evtvar[ivar].fatjet.clear();
			evtvar[ivar].jet.clear();
			evtvar[ivar].st=0;
			evtvar[ivar].nbmedium=0;
			evtvar[ivar].WFatJet=0;
			evtvar[ivar].run=run;
			evtvar[ivar].event=event;
		}
		for(unsigned int ilep=0;ilep<lepp4->size();ilep++){
			LepInfo temp;
			temp.p4=lepp4->at(ilep);
			temp.tight=leptight->at(ilep);
			temp.pdgid=pdgid->at(ilep);
			if(!LepID(temp)) continue;
			evt.lep.push_back(temp);
			evt.st+=lepp4->at(ilep).pt();
			for(unsigned int ivar=0;ivar<56;ivar++){
				evtvar[ivar].lep.push_back(temp);
				evtvar[ivar].st+=lepp4->at(ilep).pt();
			}
		}

		evt.deltaRll=ROOT::Math::VectorUtil::DeltaR(lepp4->at(0),lepp4->at(1));
		evt.mll=(lepp4->at(0)+lepp4->at(1)).M();
		for(unsigned int ivar=0;ivar<56;ivar++){
			evtvar[ivar].deltaRll=ROOT::Math::VectorUtil::DeltaR(lepp4->at(0),lepp4->at(1));
			evtvar[ivar].mll=(lepp4->at(0)+lepp4->at(1)).M();
		}
		for(unsigned int ifatjet=0;ifatjet<fatjetp4->size();ifatjet++){
			FatJetInfo temp;
			temp.p4=fatjetp4->at(ifatjet);
			temp.sdmass=fatjetsdm->at(ifatjet);
			temp.WPid=WPid->at(ifatjet);
			if(!FatJetID(temp)) continue;
			evt.st+=fatjetp4->at(ifatjet).pt();
			evt.fatjet.push_back(temp);
			if(WFatJetID(temp)){
				evt.WFatJet++;
				if(evt.WFatJet==1){
					int type = fakefatjetsfprovider.type(*genp4,*genpdgid,*genmotherid,fatjetp4->at(ifatjet));
					if(type==-1){
						fatjetSF=fjSF->at(ifatjet);
						fatjetSFup=fjSFup->at(ifatjet);
						fatjetSFdn=fjSFdn->at(ifatjet);					
					}
					else{
						int year = 10;
						if((lumi>59. && lumi<61.) || (lumi>137. && lumi<138.)){
							year = 2018;
						}
						else if(lumi>41. && lumi<43.){
							year = 2017;
						}
						else if(lumi>16. && lumi<17.){
							year = 2016;
						}
						fatjetSF=fakefatjetsfprovider.eval(year,type,2,fatjetp4->at(ifatjet).pt());
						fatjetSFup=fakefatjetsfprovider.eval_up(year,type,2,fatjetp4->at(ifatjet).pt());
						fatjetSFdn=fakefatjetsfprovider.eval_down(year,type,2,fatjetp4->at(ifatjet).pt());
					}
				}
			} 
		}
		evt.scale*=fatjetSF;
		for(unsigned int ivar=0;ivar<56;ivar++){
			for(unsigned int ifatjet=0;ifatjet<fatjetp4var[ivar]->size();ifatjet++){
				FatJetInfo tempvar;
				tempvar.p4=fatjetp4var[ivar]->at(ifatjet);
				tempvar.sdmass=fatjetsdmvar[ivar]->at(ifatjet);
				tempvar.WPid=WPid->at(ifatjet);
				if(!FatJetID(tempvar)) continue;		
				evtvar[ivar].st+=fatjetp4var[ivar]->at(ifatjet).pt();
				evtvar[ivar].fatjet.push_back(tempvar);
				if(WFatJetID(tempvar)){
					evtvar[ivar].WFatJet++;
				} 
			}
			evtvar[ivar].scale*=fatjetSF;
		}
		for(unsigned int ijet=0;ijet<jetp4->size();ijet++){
			JetInfo temp;
			temp.p4=jetp4->at(ijet);
			bool overlap=false;
			for(unsigned int ifatjet=0;ifatjet<evt.fatjet.size();ifatjet++){
				if(ROOT::Math::VectorUtil::DeltaR(jetp4->at(ijet),evt.fatjet.at(ifatjet).p4)<0.8){
					overlap=true;
					break;
				}				
			}
			temp.overlapwithfatjet=overlap;
			if(!JetID(temp)) continue;
			evt.jet.push_back(temp);
			if(!overlap){
				evt.st+=jetp4->at(ijet).pt();
			}
		}
		for(unsigned int ivar=0;ivar<56;ivar++){
			for(unsigned int ijet=0;ijet<jetp4var[ivar]->size();ijet++){
				JetInfo tempvar;
				tempvar.p4=jetp4var[ivar]->at(ijet);
				bool overlap=false;
				for(unsigned int ifatjet=0;ifatjet<evtvar[ivar].fatjet.size();ifatjet++){
					if(ROOT::Math::VectorUtil::DeltaR(jetp4var[ivar]->at(ijet),evtvar[ivar].fatjet.at(ifatjet).p4)<0.8){
						overlap=true;
						break;
					}				
				}
				tempvar.overlapwithfatjet=overlap;
				if(!JetID(tempvar)) continue;
				evtvar[ivar].jet.push_back(tempvar);
				if(!overlap){
					evtvar[ivar].st+=jetp4var[ivar]->at(ijet).pt();
				}
			}
		}
		evt.nbmedium=nbmedium;
		evt.met=met->pt();
		evt.st+=met->pt();
		evt.st=evt.st>2000?1999.5:evt.st;
		for(unsigned int ivar=0;ivar<56;ivar++){
			evtvar[ivar].nbmedium=nbmedium;
			evtvar[ivar].met=metp4var[ivar]->pt();
			evtvar[ivar].st+=metp4var[ivar]->pt();
			evtvar[ivar].st=evtvar[ivar].st>2000?1999.5:evtvar[ivar].st;
		}

		//build the event
		//fill in histograms
		if(Cut(evt)){//addition cut
			float variable_of_interest=Variable_of_Interest(evt);
			variable_of_interest=variable_of_interest>2000?1999.5:variable_of_interest;
			hcentral->Fill(variable_of_interest,				evt.scale);
			hpuWeightup->Fill(variable_of_interest,			evt.scale/puWeight*puWeightup);
			hpuWeightdn->Fill(variable_of_interest,			evt.scale/puWeight*puWeightdn);
			hprefireWeightup->Fill(variable_of_interest,		evt.scale/prefireWeight*prefireWeightup);
			hprefireWeightdn->Fill(variable_of_interest,		evt.scale/prefireWeight*prefireWeightdn);
			htriggerWeightup->Fill(variable_of_interest,		evt.scale/triggerWeight*triggerWeightup);
			htriggerWeightdn->Fill(variable_of_interest,		evt.scale/triggerWeight*triggerWeightdn);
			hlepSFelup->Fill(variable_of_interest,			evt.scale/lepSF*lepSFelup);
			hlepSFeldn->Fill(variable_of_interest,			evt.scale/lepSF*lepSFeldn);
			hlepSFmuup->Fill(variable_of_interest,			evt.scale/lepSF*lepSFmuup);
			hlepSFmudn->Fill(variable_of_interest,			evt.scale/lepSF*lepSFmudn);
			hbtagSFup->Fill(variable_of_interest,				evt.scale/btagSF*btagSFup);
			hbtagSFdn->Fill(variable_of_interest,				evt.scale/btagSF*btagSFdn);
			hfatjetSFup->Fill(variable_of_interest,			evt.scale/fatjetSF*fatjetSFup);
			hfatjetSFdn->Fill(variable_of_interest,			evt.scale/fatjetSF*fatjetSFdn);
			if(qcd->size()==9){
				hqcdscale[0]->Fill(variable_of_interest,		evt.scale*qcd->at(0));
				hqcdscale[1]->Fill(variable_of_interest,		evt.scale*qcd->at(1));
				hqcdscale[2]->Fill(variable_of_interest,		evt.scale*qcd->at(3));
				hqcdscale[3]->Fill(variable_of_interest,		evt.scale*qcd->at(5));	
				hqcdscale[4]->Fill(variable_of_interest,		evt.scale*qcd->at(7));
				hqcdscale[5]->Fill(variable_of_interest,		evt.scale*qcd->at(8));
			}
			else if(qcd->size()==8){
				hqcdscale[0]->Fill(variable_of_interest,		evt.scale*qcd->at(0));
				hqcdscale[1]->Fill(variable_of_interest,		evt.scale*qcd->at(1));
				hqcdscale[2]->Fill(variable_of_interest,		evt.scale*qcd->at(3));
				hqcdscale[3]->Fill(variable_of_interest,		evt.scale*qcd->at(4));	
				hqcdscale[4]->Fill(variable_of_interest,		evt.scale*qcd->at(6));
				hqcdscale[5]->Fill(variable_of_interest,		evt.scale*qcd->at(7));
			}
			if(pdf->size()<100){
				for(unsigned int i=0;i<100;i++){
					hpdf[i]->Fill(variable_of_interest,			evt.scale);
				}
			}
			else{
				if(pdf_replica){
					for(unsigned int i=0;i<100;i++){
						hpdf[i]->Fill(variable_of_interest,			evt.scale*((pdf->at(i+1)-1)/sqrt(99)+1));
					}
				}
				else{
					for(unsigned int i=0;i<100;i++){
						hpdf[i]->Fill(variable_of_interest,			evt.scale*pdf->at(i+1));
					}
				}
			}

			float alphas_up=1,alphas_dn=1;
			if(pdf->size()==103){
				alphas_dn=pdf->at(101);
				alphas_up=pdf->at(102);
			}
			halphasup->Fill(variable_of_interest,				evt.scale*alphas_up);
			halphasdn->Fill(variable_of_interest,				evt.scale*alphas_dn);
			/*
			if(EFTCoef->size()!=0){
				for(unsigned int i=0;i<EFTCoef->size();i++){
					hcentral_vs_EFTCoefficient->Fill(variable_of_interest,i,evt.scale / EFTCoef->at(0) * EFTCoef->at(i));
				}
			}
			*/
			if(EFT->size()!=0){
				for(unsigned iWC=0;iWC<21;iWC++){
					for(unsigned iPoint=0;iPoint<12;iPoint++){
						int iEFT=iWC*12+iPoint;
						float scale=evt.scale;
						if(index[iWC]!=0){
							scale=evt.scale/EFT->at(0)*EFT->at(index[iWC]+iPoint);	
						}
						hcentralEFT[iEFT]->Fill(variable_of_interest,				scale);
						hpuWeightupEFT[iEFT]->Fill(variable_of_interest,				scale/puWeight*puWeightup);
						hpuWeightdnEFT[iEFT]->Fill(variable_of_interest,				scale/puWeight*puWeightdn);
						hprefireWeightupEFT[iEFT]->Fill(variable_of_interest,			scale/prefireWeight*prefireWeightup);
						hprefireWeightdnEFT[iEFT]->Fill(variable_of_interest,			scale/prefireWeight*prefireWeightdn);
						htriggerWeightupEFT[iEFT]->Fill(variable_of_interest,			scale/triggerWeight*triggerWeightup);
						htriggerWeightdnEFT[iEFT]->Fill(variable_of_interest,			scale/triggerWeight*triggerWeightdn);
						hlepSFelupEFT[iEFT]->Fill(variable_of_interest,				scale/lepSF*lepSFelup);
						hlepSFeldnEFT[iEFT]->Fill(variable_of_interest,				scale/lepSF*lepSFeldn);
						hlepSFmuupEFT[iEFT]->Fill(variable_of_interest,				scale/lepSF*lepSFmuup);
						hlepSFmudnEFT[iEFT]->Fill(variable_of_interest,				scale/lepSF*lepSFmudn);
						hbtagSFupEFT[iEFT]->Fill(variable_of_interest,				scale/btagSF*btagSFup);
						hbtagSFdnEFT[iEFT]->Fill(variable_of_interest,				scale/btagSF*btagSFdn);
						hfatjetSFupEFT[iEFT]->Fill(variable_of_interest,				scale/fatjetSF*fatjetSFup);
						hfatjetSFdnEFT[iEFT]->Fill(variable_of_interest,				scale/fatjetSF*fatjetSFdn);
						if(qcd->size()==9){
							hqcdscaleEFT[0][iEFT]->Fill(variable_of_interest,		scale*qcd->at(0));
							hqcdscaleEFT[1][iEFT]->Fill(variable_of_interest,		scale*qcd->at(1));
							hqcdscaleEFT[2][iEFT]->Fill(variable_of_interest,		scale*qcd->at(3));
							hqcdscaleEFT[3][iEFT]->Fill(variable_of_interest,		scale*qcd->at(5));	
							hqcdscaleEFT[4][iEFT]->Fill(variable_of_interest,		scale*qcd->at(7));
							hqcdscaleEFT[5][iEFT]->Fill(variable_of_interest,		scale*qcd->at(8));
						}
						else{
							hqcdscaleEFT[0][iEFT]->Fill(variable_of_interest,		scale*qcd->at(0));
							hqcdscaleEFT[1][iEFT]->Fill(variable_of_interest,		scale*qcd->at(1));
							hqcdscaleEFT[2][iEFT]->Fill(variable_of_interest,		scale*qcd->at(3));
							hqcdscaleEFT[3][iEFT]->Fill(variable_of_interest,		scale*qcd->at(4));	
							hqcdscaleEFT[4][iEFT]->Fill(variable_of_interest,		scale*qcd->at(6));
							hqcdscaleEFT[5][iEFT]->Fill(variable_of_interest,		scale*qcd->at(7));
						}
						for(unsigned int i=0;i<100;i++){
							hpdfEFT[i][iEFT]->Fill(variable_of_interest,			scale*pdf->at(i+1));
						}
						halphasupEFT[iEFT]->Fill(variable_of_interest,                           	scale*alphas_up);
						halphasdnEFT[iEFT]->Fill(variable_of_interest,                           	scale*alphas_dn);
					}
				}
			}			
			//fill in histograms
		}
		for(unsigned int ivar=0;ivar<27;ivar++){
			if(Cut(evtvar[ivar])){
				float variable_of_interest=Variable_of_Interest(evtvar[ivar]);
				variable_of_interest=variable_of_interest>2000?1999.5:variable_of_interest;	
				hjesup[ivar]->Fill(variable_of_interest,				evtvar[ivar].scale);
				if(EFT->size()!=0){
					for(unsigned iWC=0;iWC<21;iWC++){
						for(unsigned iPoint=0;iPoint<12;iPoint++){
							int iEFT=iWC*12+iPoint;
							float scale=evt.scale;
							if(index[iWC]!=0){
								scale=evt.scale/EFT->at(0)*EFT->at(index[iWC]+iPoint);	
							}
							hjesupEFT[ivar][iEFT]->Fill(variable_of_interest,				scale);
						}
					}
				}
			}
			if(Cut(evtvar[ivar+28])){
				float variable_of_interest=Variable_of_Interest(evtvar[ivar+28]);
				variable_of_interest=variable_of_interest>2000?1999.5:variable_of_interest;
				hjesdn[ivar]->Fill(variable_of_interest,			evtvar[ivar+28].scale);
				if(EFT->size()!=0){
					for(unsigned iWC=0;iWC<21;iWC++){
						for(unsigned iPoint=0;iPoint<12;iPoint++){
							int iEFT=iWC*12+iPoint;
							float scale=evt.scale;
							if(index[iWC]!=0){
								scale=evt.scale/EFT->at(0)*EFT->at(index[iWC]+iPoint);	
							}
							hjesdnEFT[ivar][iEFT]->Fill(variable_of_interest,			scale);
						}
					}
				}
			}		
		}

		if(Cut(evtvar[27])){
			float variable_of_interest=Variable_of_Interest(evtvar[27]);
			variable_of_interest=variable_of_interest>2000?1999.5:variable_of_interest;
			hjerup->Fill(variable_of_interest,				evtvar[27].scale);
			if(EFT->size()!=0){
				for(unsigned iWC=0;iWC<21;iWC++){
					for(unsigned iPoint=0;iPoint<12;iPoint++){
						int iEFT=iWC*12+iPoint;
						float scale=evt.scale;
						if(index[iWC]!=0){
							scale=evt.scale/EFT->at(0)*EFT->at(index[iWC]+iPoint);	
						}
						hjerupEFT[iEFT]->Fill(variable_of_interest,				scale);
					}
				}
			}
		}
		if(Cut(evtvar[55])){
			float variable_of_interest=Variable_of_Interest(evtvar[55]);
			variable_of_interest=variable_of_interest>2000?1999.5:variable_of_interest;
			hjerdn->Fill(variable_of_interest,				evtvar[55].scale);
			if(EFT->size()!=0){
				for(unsigned iWC=0;iWC<21;iWC++){
					for(unsigned iPoint=0;iPoint<12;iPoint++){
						int iEFT=iWC*12+iPoint;
						float scale=evt.scale;
						if(index[iWC]!=0){
							scale=evt.scale/EFT->at(0)*EFT->at(index[iWC]+iPoint);	
						}
						hjerupEFT[iEFT]->Fill(variable_of_interest,				scale);
					}
				}
			}
		}
	}
	return true;
}

bool HistCollection::Write(const char* path){
	TFile *file=TFile::Open(path,"RECREATE");
	//hcentral_vs_EFTCoefficient->Write();
	hcentral->Write();
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
	hfatjetSFup->Write();
	hfatjetSFdn->Write();
	for(unsigned int i=0;i<27;i++){
		hjesup[i]->Write();
		hjesdn[i]->Write();
	}
	hjerup->Write();
	hjerdn->Write();
	halphasup->Write();
	halphasdn->Write();
	for(unsigned int i=0;i<6;i++){
		hqcdscale[i]->Write();
	}
	for(unsigned int i=0;i<100;i++){
		hpdf[i]->Write();
	}
	for(unsigned i=0;i<252;i++){
		hcentralEFT[i]->Write();
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
		hfatjetSFupEFT[i]->Write();
		hfatjetSFdnEFT[i]->Write();
		for(unsigned int j=0;j<27;j++){
			hjesupEFT[j][i]->Write();
			hjesdnEFT[j][i]->Write();
		}
		hjerupEFT[i]->Write();
		hjerdnEFT[i]->Write();
		halphasupEFT[i]->Write();
		halphasdnEFT[i]->Write();
		for(unsigned int j=0;j<6;j++){
			hqcdscaleEFT[j][i]->Write();
		}
		for(unsigned int j=0;j<100;j++){
			hpdfEFT[j][i]->Write();
		}
	}
	file->Close();
	return true;
}

class SignalRegion:public HistCollection{
	protected:
		bool LepID(const LepInfo &lep){
			return (lep.tight==1);
			//return true;
		}

		bool JetID(const JetInfo &jet){
			return (jet.p4.pt()>30.);
		}
		bool FatJetID(const FatJetInfo &fatjet){
			return (fatjet.p4.pt()>200.);
		}
		bool WFatJetID(const FatJetInfo &fatjet){
			return (fatjet.WPid>=2 && fatjet.sdmass<=105. && fatjet.sdmass>=65. && fatjet.p4.pt()>200.);
		}
		bool Cut(const Event &evt){
				if(!(evt.lep.size()==2)) return false;
				if(!(evt.fatjet.size()>=1)) return false;
				if(!(evt.WFatJet>=1)) return false;
				//if(!(evt.lep.at(0).pdgid * evt.lep.at(1).pdgid==143)) return false;
				if(!(evt.nbmedium==0)) return false;
				if(evt.lep.at(0).pdgid * evt.lep.at(1).pdgid==121 && abs(evt.mll-mZ)<20) return false;
				if(!(evt.lep.at(0).p4.pt()>40 && evt.lep.at(1).p4.pt()>30)) return false;
				if(!(evt.deltaRll>1.2)) return false;
				return true;
		}
};

class ControlRegion_ttbar:public HistCollection{
        protected:
                bool LepID(const LepInfo &lep){
                        return (lep.tight==1);
                }

		bool JetID(const JetInfo &jet){
			return (jet.p4.pt()>30.);
		}

		bool FatJetID(const FatJetInfo &fatjet){
                        return (fatjet.p4.pt()>200.);
                }
                bool WFatJetID(const FatJetInfo &fatjet){
                        return (fatjet.WPid>=2 && fatjet.sdmass<=105. && fatjet.sdmass>=65. && fatjet.p4.pt()>200.);
                }

                bool Cut(const Event &evt){
                        if(!(evt.lep.size()==2)) return false;
                        if(!(evt.WFatJet>=1)) return false;
                        if(evt.lep.at(0).pdgid * evt.lep.at(1).pdgid==121 && abs(evt.mll-mZ)<20) return false;
                        if(!(evt.lep.at(0).p4.pt()>40 && evt.lep.at(1).p4.pt()>30)) return false;
                        if(!(evt.deltaRll>1.2)) return false;
                        if(!(evt.nbmedium>=1)) return false;
                        return true;
                }
};


