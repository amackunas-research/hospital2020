#include "Heuristics.h"

#include <algorithm>

#include "Data.h"
#include "Results.h"
#include "Defines.h"

static bool Rule1(const Data &data, const int i, const int j, const int k, HeuristicsResults &res)
{
	const double z_it = (double)max(0, res.L_it[i][k] + data.p_j[j] - data.D_t[k]); // extension time

	if (res.L_it[i][k] + data.p_j[j] <= data.D_t[k] + data.delta_t[k] &&
		data.v_j[j] > data.c_jit[j][i][k] + data.e_it[i][k] * z_it)
	{
		res.x_jit[j][i][k] = 1;
		res.z_it[i][k] = z_it;
		res.L_it[i][k] += data.p_j[j];

		return true;
	}

	return false;
}

static bool Rule2(const Data &data, const int i, const int j, const int k, HeuristicsResults &res)
{
	for (int s = 1; s <= c_S; ++s)
	{
		if (!(res.L_sit[s][i][k] + data.b_sj[s][j] <= data.D_t[k] + data.delta_t[k]))
			return false;
	}

	const double z_it = [&]() // extension time
	{
		int nMaxZ_it = 0;
		for (int s = 1; s <= c_S; ++s)
			nMaxZ_it = max(nMaxZ_it, max(0, res.L_sit[s][i][k] + data.p_j[j] - data.D_t[k]));
		return nMaxZ_it;
	}();

	if (data.v_j[j] > data.c_jit[j][i][k] + data.e_it[i][k] * z_it)
	{
		res.x_jit[j][i][k] = 1;
		res.z_it[i][k] = z_it;

		for (int s = 1; s <= c_S; ++s)
			res.L_sit[s][i][k] += data.b_sj[s][j];

		return true;
	}

	return false;
}

static auto ApplyRule = c_bUMaxValue ? Rule2 : Rule1;

static void RenumberAlphaMaxVal(const Data &data, Vector2<int> &N)
{
	sort(N.begin(), N.end(), [&data](const int j1, const int j2)
	{
		return (data.p_j[j1] < data.p_j[j2]);
	});
}

static void RenumberAlphaUMaxVal(const Data &data, Vector2<int> &N)
{
	sort(N.begin(), N.end(), [&data](const int j1, const int j2)
	{
		int sum_bj1 = 0;
		int sum_bj2 = 0;

		for (int s = 1; s <= c_S; ++s)
		{
			sum_bj1 += data.b_sj[s][j1];
			sum_bj2 += data.b_sj[s][j2];
		}

		return (sum_bj1 < sum_bj2);
	});
}

static void RenumberBeta(const Data &data, Vector2<int> &T)
{
	sort(T.begin(), T.end(), [&data](const int k1, const int k2)
	{
		return k1 < k2;
	});
}

static void RenumberGamma(const Data &data, Vector2<int> &M)
{
	sort(M.begin(), M.end(), [&](const int i1, const int i2)
	{
		return (data.N_i[i1].size() < data.N_i[i2].size());
	});
}

static void RenumberRandom(Vector2<int> &v)
{
	std::random_shuffle(v.begin(), v.end());
}

static auto RenumberAlpha = c_bUMaxValue ? RenumberAlphaUMaxVal : RenumberAlphaMaxVal;

void nmHeuristics::OTR(const Data &data, const bool bRandom, HeuristicsResults &res)
{
	res.Clear();

	Vector2<int> N = data.N;

	if (bRandom)
		RenumberRandom(N);
	else
		RenumberAlpha(data, N);

	for (const int j : N)
	{
		bool bAssigned = false;

		Vector2<int> T_j = data.T_j[j];
		if (bRandom)
			RenumberRandom(T_j);
		else
			RenumberBeta(data, T_j);

		for (const int k : T_j)
		{
			Vector2<int> M_jt = data.M_jt[j][k];
			if (bRandom)
				RenumberRandom(M_jt);
			else
				RenumberGamma(data, M_jt);

			for (const int i : M_jt)
			{
				if (::ApplyRule(data, i, j, k, res))
				{
					bAssigned = true;
					break;
				}
			}

			if (bAssigned)
				break;
		}
	}
}

void nmHeuristics::OTRall(const Data &data, HeuristicsResults &res)
{
	res.Clear();

	Vector2<int> N = data.N;
	std::sort(N.begin(), N.end());

	do
	{
		for (const int j : N)
		{
			bool bAssigned = false;

			Vector2<int> T_j = data.T_j[j];
			std::sort(T_j.begin(), T_j.end());

			do
			{
				for (const int k : T_j)
				{
					Vector2<int> M_jt = data.M_jt[j][k];
					std::sort(M_jt.begin(), M_jt.end());

					do
					{
						for (const int i : M_jt)
						{
							if (::ApplyRule(data, i, j, k, res))
							{
								bAssigned = true;
								break;
							}

							if (bAssigned)
								break;
						}
					} while (std::next_permutation(M_jt.begin(), M_jt.end()));

					if (bAssigned)
						break;
				}

				if (bAssigned)
						break;
			} while (std::next_permutation(T_j.begin(), T_j.end()));
		}
	}
	while (std::next_permutation(N.begin(), N.end()));
}

void nmHeuristics::ORT(const Data &data, const bool bRandom, HeuristicsResults &res)
{
	res.Clear();

	Vector2<int> N = data.N;

	if (bRandom)
		RenumberRandom(N);
	else
		RenumberAlpha(data, N);

	for (const int j : N)
	{
		bool bAssigned = false;

		Vector2<int> M_j = data.M_j[j];

		if (bRandom)
			RenumberRandom(M_j);
		else
			RenumberGamma(data, M_j);

		for (const int i : M_j)
		{
			Vector2<int> T_ji = data.T_ji[j][i];

			if (bRandom)
				RenumberRandom(T_ji);
			else
				RenumberBeta(data, T_ji);

			for (const int k : T_ji)
			{
				if (::ApplyRule(data, i, j, k, res))
				{
					bAssigned = true;
					break;
				}
			}

			if (bAssigned)
				break;
		}
	}
}

void nmHeuristics::TRO(const Data &data, const bool bRandom, HeuristicsResults &res)
{
	res.Clear();

	vector<bool> Asgn_j(data.n + 1, false); // is operation j assigned

	Vector2<int> T = data.T;

	if (bRandom)
		RenumberRandom(T);
	else
		RenumberBeta(data, T);

	for (const int k : T)
	{
		Vector2<int> M_t = data.M_t[k];

		if (bRandom)
			RenumberRandom(M_t);
		else
			RenumberGamma(data, M_t);

		for (const int i : M_t)
		{
			Vector2<int> N_it = data.N_it[i][k];

			if (bRandom)
				RenumberRandom(N_it);
			else
				RenumberAlpha(data, N_it);

			for (const int j : N_it)
			{
				if (!Asgn_j[j] && ::ApplyRule(data, i, j, k, res))
					Asgn_j[j] = true;
			}
		}
	}
}

void nmHeuristics::TOR(const Data &data, const bool bRandom, HeuristicsResults &res)
{
	res.Clear();

	vector<bool> Asgn_j(data.n + 1, false); // is operation j assigned

	Vector2<int> T = data.T;

	if (bRandom)
		RenumberRandom(T);
	else
		RenumberBeta(data, T);

	for (const int k : T)
	{
		Vector2<int> N_t = data.N_t[k];

		if (bRandom)
			RenumberRandom(N_t);
		else
			RenumberAlpha(data, N_t);

		for (const int j : N_t)
		{
			if (Asgn_j[j])
				continue;

			Vector2<int> M_jt = data.M_jt[j][k];

			if (bRandom)
				RenumberRandom(M_jt);
			else
				RenumberGamma(data, M_jt);

			for (const int i : M_jt)
			{
				if (::ApplyRule(data, i, j, k, res))
				{
					Asgn_j[j] = true;
					break;
				}
			}
		}
	}
}

void nmHeuristics::ROT(const Data &data, const bool bRandom, HeuristicsResults &res)
{
	res.Clear();

	vector<bool> Asgn_j(data.n + 1, false); // is operation j assigned

	Vector2<int> M = data.M;

	if (bRandom)
		RenumberRandom(M);
	else
		RenumberGamma(data, M);

	for (const int i : M)
	{
		Vector2<int> N_i = data.N_i[i];

		if (bRandom)
			RenumberRandom(N_i);
		else
			RenumberAlpha(data, N_i);

		for (const int j : N_i)
		{
			if (Asgn_j[j])
				continue;

			Vector2<int> T_ji = data.T_ji[j][i];

			if (bRandom)
				RenumberRandom(T_ji);
			else
				RenumberBeta(data, T_ji);

			for (const int k : T_ji)
			{
				if (::ApplyRule(data, i, j, k, res))
				{
					Asgn_j[j] = true;
					break;
				}
			}
		}
	}
}

void nmHeuristics::RTO(const Data &data, const bool bRandom, HeuristicsResults &res)
{
	res.Clear();

	vector<bool> Asgn_j(data.n + 1, false); // is operation j assigned

	Vector2<int> M = data.M;

	if (bRandom)
		RenumberRandom(M);
	else
		RenumberGamma(data, M);

	for (const int i : M)
	{
		Vector2<int> T_i = data.T_i[i];

		if (bRandom)
			RenumberRandom(T_i);
		else
			RenumberBeta(data, T_i);

		for (const int k : T_i)
		{
			Vector2<int> N_it = data.N_it[i][k];

			if (bRandom)
				RenumberRandom(N_it);
			else
				RenumberAlpha(data, N_it);

			for (const int j : N_it)
			{
				if (!Asgn_j[j] && ::ApplyRule(data, i, j, k, res))
					Asgn_j[j] = true;
			}
		}
	}
}