namespace LeptonCuts
{


	const struct Baseline
	{
		struct Electron
		{
			float pt  = 25.0;
			float eta = 2.1;
		}Electron;

		struct Muon
		{
			float pt  = 21.0;
			float eta = 2.1;
		}Muon;

	}Baseline;

	const struct Di
	{
		struct Electron
		{
			float pt  = 15.0;
			float eta = 2.5;
		}Electron;

		struct Muon
		{
			float pt  = 15.0;
			float eta = 2.4;
		}Muon;

	}Di;

	const struct Extra
	{
		struct Electron
		{
			float pt  = 10.0;
			float eta = 2.5;
		}Electron;

		struct Muon
		{
			float pt  = 10.0;
			float eta = 2.4;
		}Muon;

	}Extra;

	const struct Additional
	{
		struct Electron
		{
			float pt  = 13.0;
			float eta = 2.1;
		}Electron;

		struct Muon
		{
			float pt  = 5.0;
			float eta = 2.1;
		}Muon;

	}Additional;


}
