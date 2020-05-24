#pragma once

#include "Data.h"
#include "Defines.h"

struct HeuristicsResults
{
	HeuristicsResults(const int n, const int t, const int m)
	{
		z_it.resize(m);
		for (auto &z_t : z_it)
			z_t.resize(t);

		x_jit.resize(n);
		for (auto &x_it : x_jit)
		{
			x_it.resize(m);
			for (auto &x_t : x_it)
				x_t.resize(t);
		}

		L_it.resize(m);
		for (auto &L_t : L_it)
			L_t.resize(t);

		L_sit.resize(c_S);
		for (auto &L_it : L_sit)
		{
			L_it.resize(m);
			for (auto &L_t : L_it)
				L_t.resize(t);
		}
	}

	void Clear()
	{
		for (auto &z_t : z_it)
			for (auto &z : z_t)
				z = 0;

		for (auto &x_it : x_jit)
			for (auto &x_t : x_it)
				for (auto &x : x_t)
					x = 0;

		for (auto &L_t : L_it)
			for (auto &L : L_t)
				L = 0;

		for (auto &L_it : L_sit)
			for (auto &L_t : L_it)
				for (auto &L : L_t)
					L = 0;
	}

	Vector2<Vector2<double>> z_it;				// the extension value of the time slot t in operating room i
	Vector2<Vector2<Vector2<uint8_t>>> x_jit;	// patient is assigned to operating room in the time slot
	Vector2<Vector2<int>> L_it;					// total time load of operating room i in time slot t
	Vector2<Vector2<Vector2<int>>> L_sit;		// umaxvalue: total time load of operating room i in time slot t
};

namespace nmResults
{
	double CalcNetValue(const Data &data, const HeuristicsResults &res);
	int CalculateIdl(const Data &data, const HeuristicsResults &res);
	std::vector<int> CalculateUdl(const Data &data, const HeuristicsResults &res);
}