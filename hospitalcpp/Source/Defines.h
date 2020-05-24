#pragma once

enum ETestOptions // flags
{
	ETO_NOACTION,
	ETO_HEURISTICS = 1,
	ETO_CPLEX = 1 << 1,
};

/////////////////////////////////////////////////////////////////////////////

#define RUN_CALCULATIONS // if commented - generate the results for the paper
#define HEURISTICS_USE_CPLEX_TIME
//#define UMAXVALUE
const ETestOptions eTestOptions = (ETestOptions)(ETO_HEURISTICS | ETO_CPLEX);

/////////////////////////////////////////////////////////////////////////////

#ifdef UMAXVALUE
	const bool c_bUMaxValue = true;
#else
	const bool c_bUMaxValue = false;
#endif

const int c_S = 3; // uncertain problem: amount of scenarios

const int nSituations = 18;
const int nInstances = 20;

const int nHeuristics = 6;

struct SSituation
{
	SSituation(const int n, const int t, const int m, const long nMaxValTime, const long nUMaxValTime) :
	n(n), t(t), m(m), nMaxVHeurTime(nMaxValTime / nHeuristics), nUMaxVHeurTime(nUMaxValTime / nHeuristics) {}

	const int n, t, m;
	const long nMaxVHeurTime; // not used
	const long nUMaxVHeurTime; // not used
};

const SSituation c_arrSituations[nSituations] = // {n, t, m, MaxValTime}
{
	{100, 10, 2, 810, 0},
	{100, 10, 4, 2258, 0},
	{100, 10, 6, 3784, 0},
	{100, 20, 2, 4652, 0},
	{100, 20, 4, 4748, 0},
	{100, 20, 6, 2293, 0},
	{500, 10, 2, 705, 0},
	{500, 10, 4, 1613, 0},
	{500, 10, 6, 2825, 0},
	{500, 20, 2, 1941, 0},
	{500, 20, 4, 6873, 0},
	{500, 20, 6, 33790, 0},
	{1000, 10, 2, 1352, 0},
	{1000, 10, 4, 3216, 0},
	{1000, 10, 6, 12388, 0},
	{1000, 20, 2, 3651, 0},
	{1000, 20, 4, 27446, 0},
	{1000, 20, 6, 69328, 0}
};
