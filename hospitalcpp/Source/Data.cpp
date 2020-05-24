#include "Data.h"

#include <algorithm>
#include <fstream>

#include "Defines.h"

static int Rand(const int nMin, const int nMax)
{
	const int nRange = nMax - nMin + 1;
	return (nMin + rand() % nRange);
}

template<unsigned int nSize>
static int Rand(const int(&arr)[nSize])
{
	return arr[rand() % nSize];
}

static void Intersection(Vector2<int> v1, Vector2<int> v2, Vector2<int> &vRez)
{
	vRez.clear();

	std::sort(v1.begin(), v1.end());
	std::sort(v2.begin(), v2.end());

	vRez.resize(min(v1.size(), v2.size()));

	auto it = std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), vRez.begin());
									
	vRez.resize(it - vRez.begin());
}

static const int P_J_UPPERBOUND = 240;

void Data::Generate(const int _n, const int _t, const int _m)
{
	// initialization
	n = _n;
	t = _t;
	m = _m;

	// help data
	Vector2<int> v_j0;
	for (int j = 1; j <= n; ++j)
	{
		const int arrRand[] = {0, 50, 100, 200};
		v_j0.push_back(Rand(arrRand));
	}

	Vector2<int> c_i0;
	for (int i = 1; i <= m; ++i)
	{
		const int arrRand[] = {0, 50, 100, 150};
		c_i0.push_back(Rand(arrRand));
	}

	Vector2<int> c_t1;
	for (int k = 1; k <= t; ++k)
	{
		const int arrRand[] = {0, 100};
		c_t1.push_back(Rand(arrRand));
	}

	// initialization
	for (int j = 1; j <= n; ++j)
		N.push_back(j);

	for (int k = 1; k <= t; ++k)
		T.push_back(k);

	for (int i = 1; i <= m; ++i)
		M.push_back(i);

	for (int k = 1; k <= t; ++k)
	{
		D_t.push_back(240);
		delta_t.push_back(60);
		M_t.push_back(M);
	}

	b_sj.resize(c_S);

	for (int j = 1; j <= n; ++j)
	{
		Vector2<int> vRandRooms = M;
		random_shuffle(vRandRooms.begin(), vRandRooms.end());
		M_j.push_back(vector<int>(vRandRooms.begin(), vRandRooms.begin() + 2));

		r_j.push_back(Rand(1, (int)ceil(4 * t / 5)));
		d_j.push_back(Rand(r_j.back(), t));

		p_j.push_back(Rand(50, P_J_UPPERBOUND));
		v_j.push_back(20 * (p_j.back() - 50) + v_j0[j]);

		Vector2<Vector2<double>> v1;
		for (int i = 1; i <= m; ++i)
		{
			Vector2<double> v2;
			for (int k = 1; k <= t; ++k)
			{
				v2.push_back(0.5 * v_j.back() + c_i0[i] + c_t1[k]);
			}
			v1.push_back(v2);
		}
		c_jit.push_back(v1);
	}

	for (int i = 1; i <= m; ++i)
	{
		Vector2<double> v1;
		for (int k = 1; k <= t; ++k)
		{
			v1.push_back(10 + 0.1 * (c_i0[i] + c_t1[k]));
		}
		e_it.push_back(v1);
	}

	M_jt.resize(n);
	for (int j = 1; j <= n; ++j)
	{
		M_jt[j].resize(t);
		for (int k = 1; k <= t; ++k)
		{
			Intersection(M_j[j], M_t[k], M_jt[j][k]);
		}
	}

	N_i.resize(m);
	for (int j = 1; j <= n; ++j)
	{
		for (auto i : M_j[j])
		{
			N_i[i].push_back(j);
		}
	}

	N_t.resize(t);
	for (int j = 1; j <= n; ++j)
	{
		for (int k = 1; k <= t; ++k)
		{
			if (r_j[j] <= k && k <= d_j[j])
			{
				N_t[k].push_back(j);
			}
		}
	}

	N_it.resize(m);
	for (int i = 1; i <= m; ++i)
	{
		N_it[i].resize(t);
		for (int k = 1; k <= t; ++k)
		{
			Intersection(N_i[i], N_t[k], N_it[i][k]);
		}
	}

	T_j.resize(n);
	for (int k = 1; k <= t; ++k)
	{
		for (auto j : N_t[k])
		{
			T_j[j].push_back(k);
		}
	}

	T_i.resize(m);
	for (int k = 1; k <= t; ++k)
	{
		for (auto i : M_t[k])
		{
			T_i[i].push_back(k);
		}
	}

	T_ji.resize(n);
	for (int j = 1; j <= n; ++j)
	{
		T_ji[j].resize(m);
		for (int i = 1; i <= m; ++i)
		{
			Intersection(T_j[j], T_i[i], T_ji[j][i]);
		}
	}
}

void Data::GenerateUncertain()
{
	// initializing b_sj for the uncertain problem
	for (int j = 1; j <= n; ++j)
	{
		const int nBjLB = p_j[j];
		const int nBjUB = min(P_J_UPPERBOUND, (int)ceil(1.05 * nBjLB));
		for (int s = 1; s <= c_S; ++s)
			b_sj[s].push_back(Rand(nBjLB, nBjUB));
	}
}

void Data::GenerateSimple()
{
	// initialization
	n = 8;
	t = 2;
	m = 3;

	// initialization
	for (int j = 1; j <= n; ++j)
		N.push_back(j);

	for (int k = 1; k <= t; ++k)
		T.push_back(k);

	for (int i = 1; i <= m; ++i)
		M.push_back(i);

	for (int k = 1; k <= t; ++k)
	{
		D_t.push_back(240);
		delta_t.push_back(60);
	}

	Vector2<int> c_i0 = vector<int> {50, 50, 50};
	Vector2<int> c_t1 = vector<int> {0, 100};

	p_j = vector<int>{120, 120, 120, 120, 250, 200, 200, 200};
	v_j = vector<int>{1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000};

	r_j = vector<int>{1, 1, 1, 1, 1, 1, 1, 1};
	d_j = vector<int>{1, 1, 1, 1, 2, 2, 2, 2};

	M_t.resize(2);
	M_t[1] = vector<int> {1, 2};
	M_t[2] = vector<int> {1, 2, 3};

	for (int j = 1; j <= n; ++j)
	{
		M_j.push_back(M);

		Vector2<Vector2<double>> v1;
		for (int i = 1; i <= m; ++i)
		{
			Vector2<double> v2;
			for (int k = 1; k <= t; ++k)
			{
				v2.push_back(0.5 * v_j.back() + c_i0[i] + c_t1[k]);
			}
			v1.push_back(v2);
		}
		c_jit.push_back(v1);
	}

	for (int i = 1; i <= m; ++i)
	{
		Vector2<double> v1;
		for (int k = 1; k <= t; ++k)
		{
			v1.push_back(10 + 0.1 * (c_i0[i] + c_t1[k]));
		}
		e_it.push_back(v1);
	}

	M_jt.resize(n);
	for (int j = 1; j <= n; ++j)
	{
		M_jt[j].resize(t);
		for (int k = 1; k <= t; ++k)
		{
			Intersection(M_j[j], M_t[k], M_jt[j][k]);
		}
	}

	N_i.resize(m);
	for (int j = 1; j <= n; ++j)
	{
		for (auto i : M_j[j])
		{
			N_i[i].push_back(j);
		}
	}

	N_t.resize(t);
	for (int j = 1; j <= n; ++j)
	{
		for (int k = 1; k <= t; ++k)
		{
			if (r_j[j] <= k && k <= d_j[j])
			{
				N_t[k].push_back(j);
			}
		}
	}

	N_it.resize(m);
	for (int i = 1; i <= m; ++i)
	{
		N_it[i].resize(t);
		for (int k = 1; k <= t; ++k)
		{
			Intersection(N_i[i], N_t[k], N_it[i][k]);
		}
	}

	T_j.resize(n);
	for (int k = 1; k <= t; ++k)
	{
		for (auto j : N_t[k])
		{
			T_j[j].push_back(k);
		}
	}

	T_i.resize(m);
	for (int k = 1; k <= t; ++k)
	{
		for (auto i: M_t[k])
		{
			T_i[i].push_back(k);
		}
	}

	T_ji.resize(n);
	for (int j = 1; j <= n; ++j)
	{
		T_ji[j].resize(m);
		for (int i = 1; i <= m; ++i)
		{
			Intersection(T_j[j], T_i[i], T_ji[j][i]);
		}
	}
}

void Data::Output(const char *pszFile) const
{
	ofstream rf(pszFile);

	auto fOutNumbers = [&](const auto &v)
	{
		for (size_t i = 1; i <= v.size(); ++i)
		{
			rf << v[i];
			if (i != v.size())
				rf << " ";
		}
	};

	auto fOutVector = [&](const auto &v)
	{
		rf << "[";
		fOutNumbers(v);
		rf << "]";
	};

	auto fOutSet = [&](const auto &v)
	{
		rf << "{";
		fOutNumbers(v);
		rf << "}";
	};

	auto fOutVectVect = [&](const auto &v)
	{
		rf << "[";
		for (size_t i = 1; i <= v.size(); ++i)
		{
			fOutVector(v[i]);
			if (i != v.size())
				rf << endl;
		}
		rf << "]";
	};

	auto fOutVectSet = [&](const auto &v)
	{
		rf << "[";
		for (size_t i = 1; i <= v.size(); ++i)
		{
			fOutSet(v[i]);
			if (i != v.size())
				rf << " ";
		}
		rf << "]";
	};

	auto fOutVectVectVect = [&](const auto &v)
	{
		rf << "[";
		for (size_t i = 1; i <= v.size(); ++i)
		{
			fOutVectVect(v[i]);
			if (i != v.size())
				rf << endl;
		}
		rf << "]";
	};

	auto fOutVectVectSet = [&](const auto &v)
	{
		rf << "[";
		for (size_t i = 1; i <= v.size(); ++i)
		{
			fOutVectSet(v[i]);
			if (i != v.size())
				rf << endl;
		}
		rf << "]";
	};

	rf << "n = " << n << ";" << endl;
	rf << "t = " << t << ";" << endl;
	rf << "m = " << m << ";" << endl;

	if (c_bUMaxValue)
		rf << "S = " << c_S << ";" << endl;

	if (c_bUMaxValue)
	{
		rf << "b_sj = ";
		fOutVectVect(b_sj);
		rf << ";" << endl;
	}
	else
	{
		rf << "p_j = ";
		fOutVector(p_j);
		rf << ";" << endl;
	}

	rf << "v_j = ";
	fOutVector(v_j);
	rf << ";" << endl;

	rf << "r_j = ";
	fOutVector(r_j);
	rf << ";" << endl;

	rf << "d_j = ";
	fOutVector(d_j);
	rf << ";" << endl;

	rf << "D_t = ";
	fOutVector(D_t);
	rf << ";" << endl;

	rf << "delta_t = ";
	fOutVector(delta_t);
	rf << ";" << endl;

	rf << "M_t = ";
	fOutVectSet(M_t);
	rf << ";" << endl;

	rf << "M_j = ";
	fOutVectSet(M_j);
	rf << ";" << endl;
	
	rf << "M_jt = ";
	fOutVectVectSet(M_jt);
	rf << ";" << endl;

	rf << "N_i = ";
	fOutVectSet(N_i);
	rf << ";" << endl;

	rf << "N_t = ";
	fOutVectSet(N_t);
	rf << ";" << endl;

	rf << "N_it = ";
	fOutVectVectSet(N_it);
	rf << ";" << endl;

	rf << "T_j = ";
	fOutVectSet(T_j);
	rf << ";" << endl;

	rf << "T_i = ";
	fOutVectSet(T_i);
	rf << ";" << endl;

	rf << "T_ji = ";
	fOutVectVectSet(T_ji);
	rf << ";" << endl;

	rf << "c_jit = ";
	fOutVectVectVect(c_jit);
	rf << ";" << endl;

	rf << "e_it = ";
	fOutVectVect(e_it);
	rf << ";" << endl;
}
