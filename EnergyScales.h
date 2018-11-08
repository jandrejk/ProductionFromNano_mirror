const struct ES{
	struct Electron{
		float oneProng0p0 = 0.01;
		float oneProng1p0 = 0.041;
		float threeProng0p0 = 0.0;
		float uncertaintyShift1p   = 0.014;
		float uncertaintyShift1p1p = 0.018;
		float uncertaintyShift3p   = 0.03;
	} Electron;

	struct Muon{
		float oneProng0p0 = 0.0;
		float oneProng1p0 = 0.0;
		float threeProng0p0 = 0.0;
		float uncertaintyShift1p   = 0.02;
		float uncertaintyShift1p1p = 0.02;
		float uncertaintyShift3p   = 0.02;
	} Muon;

	struct Tau{
		float oneProng0p0 =   0.007;
		float oneProng1p0 =  -0.002;
		float threeProng0p0 = 0.001;
		float uncertaintyShift1p   = 0.008;
		float uncertaintyShift1p1p = 0.008;
		float uncertaintyShift3p   = 0.009;
	} Tau;
} ES;
