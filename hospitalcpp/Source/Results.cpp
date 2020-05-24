#include <algorithm>

#include "Results.h"

static int CalculateIdl(const Data &data, const Vector2<Vector2<int>> &L_it)
{
	int nIdl = 0;

	for (int i = 1; i <= data.m; ++i)
	{
		for (int k = 1; k <= data.t; ++k)
		{
			const int nRoomIdl = std::max(0, data.D_t[k] - L_it[i][k]);
			nIdl += nRoomIdl;
		}
	}

	return nIdl;
}

double nmResults::CalcNetValue(const Data &data, const HeuristicsResults &res)
{
	double dRevenue = 0;
	for (int j : data.N)
	{
		for (int k : data.T_j[j])
		{
			for (int i : data.M_jt[j][k])
			{
				dRevenue += (data.v_j[j] - data.c_jit[j][i][k]) * res.x_jit[j][i][k];
			}
		}
	}

	double dCosts = 0;

	for (int k : data.T)
	{
		for (int i : data.M_t[k])
		{
			dCosts += data.e_it[i][k] * res.z_it[i][k];
		}
	}

	return dRevenue - dCosts;
}

int nmResults::CalculateIdl(const Data &data, const HeuristicsResults &res)
{
	return ::CalculateIdl(data, res.L_it);
}

std::vector<int> nmResults::CalculateUdl(const Data &data, const HeuristicsResults &res)
{
	std::vector<int> vUdl;

	for (int s = 1; s <= c_S; ++s)
		vUdl.push_back(::CalculateIdl(data, res.L_sit[s]));

	return std::move(vUdl);
}