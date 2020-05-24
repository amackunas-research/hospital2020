#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <functional>
#include <string>

#include "Data.h"
#include "Heuristics.h"
#include "Results.h"
#include "Defines.h"

static const char *pszCalcResFile = "Data/calcresults.inl";
static const char *pszHeurResFile = "Data/heuristicsres.inl";
static const char *pszDatFile = "../hospitalcplex/hospital/hospital.dat";
static const char *pszProject = "../hospitalcplex/hospital";
#define pszResultsMaxVal "..\Data\maxvalue.inl"
#define pszResultsUMaxVal "..\Data\umaxvalue.inl"

static void SetRandomSeed(const int nSituation, const int nInstance)
{
	srand(nInstance);
}

static void LaunchHeuristics(const Data &data, const int n, const int t, const int m, const int nSituation, const int nInstance)
{
	HeuristicsResults res(data.n, data.t, data.m);

	double V_D = 0; // best deterministic heuristic net value
	double V_R = 0; // best randomised heuristic net value
	int nIdl = 0;
	std::vector<int> vUdl;

	ofstream hr(pszHeurResFile, std::ios_base::app);
	hr << "//(" << n << "," << t << "," << m << ") #" << nInstance << endl;

#ifdef HEURISTICS_USE_CPLEX_TIME
	double time[nSituations][nInstances];
	{
		#ifdef UMAXVALUE
			double V_s[nSituations][nInstances];
			int udl[c_S][nSituations][nInstances];
			#include pszResultsUMaxVal
		#else
			double V_1[nSituations][nInstances];
			int idl[nSituations][nInstances];
			#include pszResultsMaxVal
		#endif
	}
#endif // HEURISTICS_USE_CPLEX_TIME

	int nHeuristicNr = 0; // just for output of each heuristic results

	auto fLaunchHeuristic = [&](const SHeuristic &heuristic, const bool bRandom, 
		double &dNetValue)
	{
#ifdef HEURISTICS_USE_CPLEX_TIME
		const clock_t timeLimit = (clock_t)time[nSituation][nInstance];
#else
		//const clock_t timeLimit = (clock_t)c_arrSituations[nSituation].nMaxVHeurTime;
		const clock_t timeLimit = (clock_t)5*60*1000;
#endif // HEURISTICS_USE_CPLEX_TIME

		int nTimes = 0;
		double dHeuristicVal = 0;
		const clock_t timeStart = clock();
		do
		{
			heuristic.fHeuristic(data, bRandom, res);
			const double dVal = nmResults::CalcNetValue(data, res);

			if (dVal > dHeuristicVal + 0.00001)
			{
				dHeuristicVal = dVal;

				if (dHeuristicVal > dNetValue + 0.00001)
				{
					dNetValue = dHeuristicVal;
					if (c_bUMaxValue)
						vUdl = nmResults::CalculateUdl(data, res);
					else
						nIdl = nmResults::CalculateIdl(data, res);
				}
			}
			++nTimes;
		}
		while (bRandom && clock() - timeStart < timeLimit);
		//if (bRandom)
		//	cout << "Rand heuristic has been run: " << nTimes << " times" << endl;

		hr << "V_" << (bRandom ? "R_" : "D_") << heuristic.sName <<
			"[" << nSituation << "][" << nInstance << "] = " << dHeuristicVal << ";" << endl;

		//// intermediate results of one of the heuristics
		//if (bRandom)
		//{
		//	hr << "V_sR[" << nSituation << "][" << nInstance << "] = " << V_R << ";" << endl;
		//	for (int s = 0; s < vUdl.size(); ++s)
		//		hr << "hudl[" << s << "][" << nSituation << "][" << nInstance << "] = " << vUdl[s] << ";" << endl;
		//	hr << endl;
		//}

		//cout << "Heuristic " << heuristic.sName << " calculated." << endl;
	};

	nHeuristicNr = 0;
	fLaunchHeuristic(H_OTR, false, V_D);
	fLaunchHeuristic(H_ORT, false, V_D);
	fLaunchHeuristic(H_TRO, false, V_D);
	fLaunchHeuristic(H_TOR, false, V_D);
	fLaunchHeuristic(H_ROT, false, V_D);
	fLaunchHeuristic(H_RTO, false, V_D);

	nHeuristicNr = 0;
	fLaunchHeuristic(H_OTR, true, V_R);
	fLaunchHeuristic(H_ORT, true, V_R);
	fLaunchHeuristic(H_TRO, true, V_R);
	fLaunchHeuristic(H_TOR, true, V_R);
	fLaunchHeuristic(H_ROT, true, V_R);
	fLaunchHeuristic(H_RTO, true, V_R);

	const double V_0 = max(V_D, V_R); // best heuristic net value

	// outputting final heuristics results
	{
		if (!c_bUMaxValue)
		{
			hr << "V_D[" << nSituation << "][" << nInstance << "] = " << V_D << ";" << endl;
			hr << "V_R[" << nSituation << "][" << nInstance << "] = " << V_R << ";" << endl;
			hr << "V_0[" << nSituation << "][" << nInstance << "] = " << V_0 << ";" << endl;
			hr << "hidl[" << nSituation << "][" << nInstance << "] = " << nIdl << ";" << endl;
		}
		else
		{
			hr << "V_sD[" << nSituation << "][" << nInstance << "] = " << V_D << ";" << endl;
			hr << "V_sR[" << nSituation << "][" << nInstance << "] = " << V_R << ";" << endl;
			hr << "V_s0[" << nSituation << "][" << nInstance << "] = " << V_0 << ";" << endl;
			for (int s = 0; s < vUdl.size(); ++s)
				hr << "hudl[" << s << "][" << nSituation << "][" << nInstance << "] = " << vUdl[s] << ";" << endl;
		}

		hr << endl;
		hr.close();
	}
}

static void LaunchCplex(const Data &data, const int nSituation, const int nInstance)
{
	const char *pszConfMaxValue = "MaxValue";
	const char *pszConfUMaxValue = "UMaxValue";

	data.Output(pszDatFile);

	// appending dat file to import it in cplex and print results
	{
		ofstream df(pszDatFile, std::ios_base::app);
		df << "situation = " << nSituation << ";" << endl;
		df << "instance = " << nInstance << ";" << endl;
		df.close();
	}

	string sCommand;
	sCommand = sCommand + "oplrun -p " + pszProject + " " + (c_bUMaxValue ? pszConfUMaxValue : pszConfMaxValue);
	system(sCommand.c_str());
}

static void LaunchTest(const int n, const int t, const int m, const int nSituation, const int nInstance)
{
	SetRandomSeed(nSituation, nInstance);

	Data data;
	data.Generate(n, t, m);
	if (c_bUMaxValue)
		data.GenerateUncertain();

	if (eTestOptions & ETO_HEURISTICS)
		::LaunchHeuristics(data, n, t, m, nSituation, nInstance);

	if (eTestOptions & ETO_CPLEX)
		::LaunchCplex(data, nSituation, nInstance);
}

static void GetMinMax(double arrVals[], int nSize, int &iMin, int &iMax)
{
	iMin = 0, iMax = 0;

	if (nSize < 0)
		return;

	double dMin = arrVals[0];
	double dMax = arrVals[0];

	for (int i = 1; i < nSize; ++i)
	{
		if (arrVals[i] < dMin)
		{
			dMin = arrVals[i];
			iMin = i;
		}

		if (arrVals[i] > dMax)
		{
			dMax = arrVals[i];
			iMax = i;
		}
	}
}

static double CalcTrimmedMean(double arrVal[nInstances], double arrDivisor[nInstances])
{
	double dTrimmedMean = 0;
	double dMinVal = DBL_MAX;
	double dMaxVal = -DBL_MAX;

	for (int i = 0; i < nInstances; ++i)
	{
		const double dVal = arrVal[i] / arrDivisor[i];

		dTrimmedMean += dVal;
		dMinVal = min(dMinVal, dVal);
		dMaxVal = max(dMaxVal, dVal);
	}

	dTrimmedMean -= (dMinVal + dMaxVal);

	dTrimmedMean /= (nInstances - 2);

	return dTrimmedMean;
};

static void OutputMaxValue()
{
#ifndef RUN_CALCULATIONS
	double V_D[nSituations][nInstances];
	double V_R[nSituations][nInstances];
	double V_0[nSituations][nInstances];
	double V_1[nSituations][nInstances];
	double time[nSituations][nInstances];
	int idl[nSituations][nInstances];
	int hidl[nSituations][nInstances];

	double V_D_OTR[nSituations][nInstances];
	double V_D_ORT[nSituations][nInstances];
	double V_D_TRO[nSituations][nInstances];
	double V_D_TOR[nSituations][nInstances];
	double V_D_ROT[nSituations][nInstances];
	double V_D_RTO[nSituations][nInstances];
	double V_R_OTR[nSituations][nInstances];
	double V_R_ORT[nSituations][nInstances];
	double V_R_TRO[nSituations][nInstances];
	double V_R_TOR[nSituations][nInstances];
	double V_R_ROT[nSituations][nInstances];
	double V_R_RTO[nSituations][nInstances];

	#include pszResultsMaxVal
	#include "..\Data\hmaxvalue.inl"

	ofstream rf1("../results/maxvalue/results1.txt");
	ofstream rf2("../results/maxvalue/results2.txt");

	rf1 << fixed;
	rf2 << fixed;

	for (int i = 0; i < nSituations; ++i)
	{
		double dTimeTM = 0;
		double dIdlTM = 0;
		double dHIdlTM = 0;
		double dVDV_TM = 0;
		double dVRV_TM = 0;
		double dV0V_TM = 0;

		double dV0V_bad = 0;

		const int To = c_arrSituations[i].m * c_arrSituations[i].t * 240;

		int jMinTime, jMaxTime; // these indexes will be skipped during trimmed mean calculation
		::GetMinMax(time[i], nInstances, jMinTime, jMaxTime);

		const double dTimeMin = time[i][jMinTime];
		const double dTimeMax = time[i][jMaxTime];

		double dIdlMin = DBL_MAX, dIdlMax = -DBL_MAX;
		double dHIdlMin = DBL_MAX, dHIdlMax = -DBL_MAX;
		double dVD1Min = DBL_MAX, dVD1Max = -DBL_MAX;
		double dVR1Min = DBL_MAX, dVR1Max = -DBL_MAX;
		double dV01Min = DBL_MAX, dV01Max = -DBL_MAX;

		for (int j = 0; j < nInstances; ++j)
		{
			const double dIdlTo = (double)idl[i][j] / To;
			const double dHIdlTo = (double)hidl[i][j] / To;
			const double dVD1 = V_D[i][j] / V_1[i][j];
			const double dVR1 = V_R[i][j] / V_1[i][j];
			const double dV01 = V_0[i][j] / V_1[i][j];

			dTimeTM += time[i][j];

			dIdlTM += dIdlTo;
			dIdlMin = min(dIdlMin, dIdlTo);
			dIdlMax = max(dIdlMax, dIdlTo);

			dHIdlTM += dHIdlTo;
			dHIdlMin = min(dHIdlMin, dHIdlTo);
			dHIdlMax = max(dHIdlMax, dHIdlTo);

			dVDV_TM += dVD1;
			dVD1Min = min(dVD1Min, dVD1);
			dVD1Max = max(dVD1Max, dVD1);

			dVRV_TM += dVR1;
			dVR1Min = min(dVR1Min, dVR1);
			dVR1Max = max(dVR1Max, dVR1);

			dV0V_TM += dV01;
			dV01Min = min(dV01Min, dV01);
			dV01Max = max(dV01Max, dV01);

			if (j == jMaxTime)
				dV0V_bad = dV01;
		}

		dTimeTM -= (dTimeMin + dTimeMax);
		dIdlTM -= (dIdlMin + dIdlMax);
		dHIdlTM -= (dHIdlMin + dHIdlMax);
		dVDV_TM -= (dVD1Min + dVD1Max);
		dVRV_TM -= (dVR1Min + dVR1Max);
		dV0V_TM -= (dV01Min + dV01Max);

		dTimeTM /= (nInstances - 2);
		dIdlTM /= (nInstances - 2);
		dHIdlTM /= (nInstances - 2);
		dVDV_TM /= (nInstances - 2);
		dVRV_TM /= (nInstances - 2);
		dV0V_TM /= (nInstances - 2);

		const double dVD_OTR_TM = CalcTrimmedMean(V_D_OTR[i], V_1[i]);
		const double dVD_ORT_TM = CalcTrimmedMean(V_D_ORT[i], V_1[i]);
		const double dVD_TRO_TM = CalcTrimmedMean(V_D_TRO[i], V_1[i]);
		const double dVD_TOR_TM = CalcTrimmedMean(V_D_TOR[i], V_1[i]);
		const double dVD_ROT_TM = CalcTrimmedMean(V_D_ROT[i], V_1[i]);
		const double dVD_RTO_TM = CalcTrimmedMean(V_D_RTO[i], V_1[i]);
		const double dVR_OTR_TM = CalcTrimmedMean(V_R_OTR[i], V_1[i]);
		const double dVR_ORT_TM = CalcTrimmedMean(V_R_ORT[i], V_1[i]);
		const double dVR_TRO_TM = CalcTrimmedMean(V_R_TRO[i], V_1[i]);
		const double dVR_TOR_TM = CalcTrimmedMean(V_R_TOR[i], V_1[i]);
		const double dVR_ROT_TM = CalcTrimmedMean(V_R_ROT[i], V_1[i]);
		const double dVR_RTO_TM = CalcTrimmedMean(V_R_RTO[i], V_1[i]);

		rf1 << setprecision(0);
		rf1 << "(" << c_arrSituations[i].n << "," << c_arrSituations[i].t << "," << c_arrSituations[i].m << ")" << endl;
		rf1 << "Trimmed mean time: " << dTimeTM << endl;
		rf1 << "Worst time: " << dTimeMax << endl;
		rf1 << setprecision(4);
		rf1 << "(Idl/To)t.m.: " << dIdlTM << endl;
		rf1 << "(Idl/To)wor: " << dIdlMax << endl;
		rf1 << endl;

		rf2 << setprecision(4);
		rf2 << "(" << c_arrSituations[i].n << "," << c_arrSituations[i].t << "," << c_arrSituations[i].m << ")" << endl;
		rf2 << "(VD/V*)t.m.: " << dVDV_TM << endl;
		rf2 << "(VR/V*)t.m.: " << dVRV_TM << endl;
		rf2 << "(V0/V*)t.m.: " << dV0V_TM << endl;
		rf2 << "(V0/V*)wor: " << dV01Min << endl;
		rf2 << "V0bad/V*bad: " << dV0V_bad << endl;
		rf2 << "(Idl/To)t.m.: " << dHIdlTM << endl;
		rf2 << "(Idl/To)wor: " << dHIdlMax << endl;
		//rf2 << "----- Heuristics quality -----" << endl;
		//rf2 << "(D_OTR)t.m.: " << dVD_OTR_TM << endl;
		//rf2 << "(D_ORT)t.m.: " << dVD_ORT_TM << endl;
		//rf2 << "(D_TRO)t.m.: " << dVD_TRO_TM << endl;
		//rf2 << "(D_TOR)t.m.: " << dVD_TOR_TM << endl;
		//rf2 << "(D_ROT)t.m.: " << dVD_ROT_TM << endl;
		//rf2 << "(D_RTO)t.m.: " << dVD_RTO_TM << endl;
		//rf2 << "(R_OTR)t.m.: " << dVR_OTR_TM << endl;
		//rf2 << "(R_ORT)t.m.: " << dVR_ORT_TM << endl;
		//rf2 << "(R_TRO)t.m.: " << dVR_TRO_TM << endl;
		//rf2 << "(R_TOR)t.m.: " << dVR_TOR_TM << endl;
		//rf2 << "(R_ROT)t.m.: " << dVR_ROT_TM << endl;
		//rf2 << "(R_RTO)t.m.: " << dVR_RTO_TM << endl;
		rf2 << endl;
	}
#endif
}

static void OutputUMaxValue()
{
#ifndef RUN_CALCULATIONS
	double V_1[nSituations][nInstances];	// maxvalue optimal
	double V_sD[nSituations][nInstances];	// umaxvalue deterministic heuristic
	double V_sR[nSituations][nInstances];	// umaxvalue randomized heuristic
	double V_s0[nSituations][nInstances];	// umaxvalue heuristic best
	double V_s[nSituations][nInstances];	// umaxvalue optimal
	double time[nSituations][nInstances];
	int udl[c_S][nSituations][nInstances];
	int hudl[c_S][nSituations][nInstances];

	double V_D_OTR[nSituations][nInstances];
	double V_D_ORT[nSituations][nInstances];
	double V_D_TRO[nSituations][nInstances];
	double V_D_TOR[nSituations][nInstances];
	double V_D_ROT[nSituations][nInstances];
	double V_D_RTO[nSituations][nInstances];
	double V_R_OTR[nSituations][nInstances];
	double V_R_ORT[nSituations][nInstances];
	double V_R_TRO[nSituations][nInstances];
	double V_R_TOR[nSituations][nInstances];
	double V_R_ROT[nSituations][nInstances];
	double V_R_RTO[nSituations][nInstances];

	{
		int idl[nSituations][nInstances];

		#include pszResultsMaxVal
		#include pszResultsUMaxVal
		#include "..\Data\humaxvalue.inl"
	}

	ofstream rf1("../results/umaxvalue/uresults1.txt");
	ofstream rf2("../results/umaxvalue/uresults2.txt");

	rf1 << fixed;
	rf2 << fixed;

	for (int i = 0; i < nSituations; ++i)
	{
		double dTimeTM = 0;
		double dUdlTM = 0;
		double dHUdlTM = 0;
		double dVp_Vps_TM = 0;
		double dU0_Us_TM = 0;

		double dU0_Us_bad = 0;

		const int To = c_arrSituations[i].m * c_arrSituations[i].t * 240;

		int jMinTime, jMaxTime; // these indexes will be skipped during trimmed mean calculation
		::GetMinMax(time[i], nInstances, jMinTime, jMaxTime);

		const double dTimeMin = time[i][jMinTime], dTimeMax = time[i][jMaxTime];
		double dUdlMin = DBL_MAX, dUdlMax = -DBL_MAX;
		double dHUdlMin = DBL_MAX, dHUdlMax = -DBL_MAX;
		double dVpVpsMin = DBL_MAX, dVpVpsMax = -DBL_MAX;
		double dU0UsMin = DBL_MAX, dU0UsMax= -DBL_MAX;

		for (int j = 0; j < nInstances; ++j)
		{
			dTimeTM += time[i][j];

			for (int s = 0; s < c_S; ++s)
			{
				const double dUdlTo = (double)udl[s][i][j] / To;
				const double dHUdlTo = (double)hudl[s][i][j] / To;

				dUdlTM += dUdlTo;
				dUdlMin = min(dUdlMin, dUdlTo);
				dUdlMax = max(dUdlMax, dUdlTo);

				dHUdlTM += dHUdlTo;
				dHUdlMin = min(dHUdlMin, dHUdlTo);
				dHUdlMax = max(dHUdlMax, dHUdlTo);
			}

			const double dVps_Vp = V_s[i][j] / V_1[i][j];

			dVp_Vps_TM += dVps_Vp;
			dVpVpsMin = min(dVpVpsMin, dVps_Vp);
			dVpVpsMax = max(dVpVpsMax, dVps_Vp);

			const double dU0_Us = V_s0[i][j] / V_s[i][j];

			dU0_Us_TM += dU0_Us;
			dU0UsMin = min(dU0UsMin, dU0_Us);
			dU0UsMax = max(dU0UsMax, dU0_Us);

			if (j == jMaxTime)
				dU0_Us_bad = dU0_Us;
		}

		dTimeTM -= (dTimeMin + dTimeMax);
		dUdlTM -= (dUdlMin + dUdlMax);
		dHUdlTM -= (dHUdlMin + dHUdlMax);
		dVp_Vps_TM -= (dVpVpsMin + dVpVpsMax);
		dU0_Us_TM -= (dU0UsMin + dU0UsMax);

		dTimeTM /= (nInstances - 2);
		dUdlTM = dUdlTM / (nInstances - 2) / c_S;
		dHUdlTM = dHUdlTM / (nInstances - 2) / c_S;
		dVp_Vps_TM /= (nInstances - 2);
		dU0_Us_TM /= (nInstances - 2);

		rf1 << setprecision(0);
		rf1 << "(" << c_arrSituations[i].n << "," << c_arrSituations[i].t << "," << c_arrSituations[i].m << ")" << endl;
		rf1 << "Trimmed mean time: " << dTimeTM << endl;
		rf1 << "Worst time: " << dTimeMax << endl;
		rf1 << setprecision(4);
		rf1 << "(Vpb/Vp*)t.m.: " << dVp_Vps_TM << endl;
		rf1 << "(Vpb/Vp*)wor: " << dVpVpsMin << endl;
		rf1 << "(Udl/To)t.m.: " << dUdlTM << endl;
		rf1 << "(Udl/To)wor: " << dUdlMax << endl;
		rf1 << endl;

		rf2 << setprecision(4);
		rf2 << "(" << c_arrSituations[i].n << "," << c_arrSituations[i].t << "," << c_arrSituations[i].m << ")" << endl;
		rf2 << "(U0/U*)t.m.: " << dU0_Us_TM << endl;
		rf2 << "(U0/U*)wor: " << dU0UsMin << endl;
		rf2 << "U0bad/U*bad: " << dU0_Us_bad << endl;
		rf2 << "(Udl/To)t.m.: " << dHUdlTM << endl;
		rf2 << "(Udl/To)wor: " << dHUdlMax << endl;
		rf2 << endl;
	}
#endif
}

static void GenerateResults()
{
#ifdef RUN_CALCULATIONS
	cout << "0 - exit, 1 - continue, 2 - regenerate:" << endl;
	int nGenerate = 0;
	cin >> nGenerate;

	int nContinueSit = 0;
	int nContinueInst = 0;

	switch (nGenerate)
	{
		case 0:
			return;
		case 1:
			cout << "Situation and instance (0-based): ";
			cin >> nContinueSit >> nContinueInst;
			break;
		case 2: // clearing files
		{
			ofstream rf;
			rf.open(pszCalcResFile); rf.close();
			rf.open(pszHeurResFile); rf.close();
			break;
		}
	}

	for (int i = nContinueSit; i < nSituations; ++i)
	{
		for (int j = 0; j < nInstances; ++j)
		{
			if (i == nContinueSit && j < nContinueInst)
				continue;

			const int n = c_arrSituations[i].n;
			const int t = c_arrSituations[i].t;
			const int m = c_arrSituations[i].m;

			cout << "Generating sit " << i << "(" << n << "," << t << "," << m << ") #" << j << endl;

			LaunchTest(n, t, m, i, j);
		}
	}

	cout << "Results have been calculated. Close the program." << endl;
#endif // RUN_CALCULATIONS
}

int main()
{
	::GenerateResults();

	if (c_bUMaxValue)
		::OutputUMaxValue();
	else
		::OutputMaxValue();

	return 0;
}