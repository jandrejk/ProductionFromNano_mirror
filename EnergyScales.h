const struct ES{
	struct Electron{
		float oneProng0p0 = 0.0;
		float oneProng1p0 = 0.0;
		float threeProng0p0 = 0.0;
		float uncertaintyShift = 0.03;
	} Electron;

	struct Muon{
		float oneProng0p0 = 0.0;
		float oneProng1p0 = 0.0;
		float threeProng0p0 = 0.0;
		float uncertaintyShift = 0.03;
	} Muon;

	struct Tau{
		float oneProng0p0 =   0.008;
		float oneProng1p0 =  -0.008;
		float threeProng0p0 = 0.010;
		float uncertaintyShift = 0.03;
	} Tau;
} ES;
