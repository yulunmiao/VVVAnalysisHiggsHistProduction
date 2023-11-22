#include "make_hist_higgs_combine.h"

int SR_make_hist_EFT_higgs_combine(){
	//WWW
	SignalRegion *VVVJ_n=new SignalRegion;
	VVVJ_n->Fill(137.64,0.051998,"WWW_1Jet.root");
        VVVJ_n->Fill(137.64,0.032426,"WWZ_1Jet.root");
        VVVJ_n->Fill(137.64,0.010671,"WZZ_1Jet.root");
        VVVJ_n->Fill(137.64,0.003822,"ZZZ_1Jet.root");
	VVVJ_n->Write("../hist/SR_higgs/all/VVV_1Jet.root");
	delete VVVJ_n;
	//WWW
	/*
	//WWZ
	SignalRegion *WWZ=new SignalRegion;
	WWZ->Fill(LUMI2018,0.003448,"./WWZ_EFT.root");
	WWZ->Write("../hist/SR/2018/hWWZ_EFT.root");
	delete WWZ;
	//WWZ
	*/
	return 0;
}
