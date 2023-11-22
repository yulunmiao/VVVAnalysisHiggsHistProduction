int extract_pdf_qcd(){
	TFile* fVVV=new TFile("../hist/SR_higgs/all/VVV_1Jet.root");
	TH1F*  hpdf[100];
	TH1F*  hqcd[6];
	TH1F*  hVVV=(TH1F*) fVVV->Get("htotal");
	for(unsigned int i=0;i<100;i++){
		hpdf[i]=(TH1F*) fVVV->Get(TString::Format("pdf_%d",i));
	}
	for(unsigned int i=0;i<6;i++){
		hqcd[i]=(TH1F*) fVVV->Get(TString::Format("qcd_%d",i));
	}
	for(unsigned int i=1;i<=3;i++){
		cout<<hVVV->GetBinContent(i)<<endl;
		float qcdup=0,qcddn=0;
		for(unsigned int j=0;j<6;j++){
                	float temp=hqcd[j]->GetBinContent(i)-hVVV->GetBinContent(i);
			if (temp>qcdup) qcdup=temp;
			if (temp<qcddn) qcddn=temp;
        	}
		cout<<qcdup<<"	"<<qcddn<<endl;
		float pdf=0;
		for(unsigned int j=0;j<100;j++){
			pdf += (hpdf[i]->GetBinContent(i)-hVVV->GetBinContent(i))*(hpdf[i]->GetBinContent(i)-hVVV->GetBinContent(i));
		}
		pdf = sqrt(pdf);
		cout<<pdf<<endl;
	}
	return 0;
}

