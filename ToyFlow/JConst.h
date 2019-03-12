#define NC 6
static double CentBins[NC+3] = {0,5,10,20,30,40,50, 60, 70};

enum DETECTOR{
	D_TPC, //TPC full coverage
	D_TPC_ETAA, //TPC with eta gap
	D_TPC_ETAC,
	D_V0A,
	D_V0C,
	D_V0P, //V0+
	D_COUNT
};
static double cov[D_COUNT][2] = {
	{-1.5,1.5},
	{-1.5,-0.4},
	{0.4,1.5},
	{2.8,5.1},
	{-3.7,-1.7},
	{2.19,5.08}
};

enum RESOLUTION{
	R_V0A,
	R_V0C,
	R_V0P,
	R_COUNT
};
static const char *presn[] = {"V0A","V0C","V0P"};