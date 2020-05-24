#pragma once

#include <vector>

using namespace std;

template <typename T>
class Vector2 : public vector<T>
{
public:
	Vector2() { }
	Vector2(const vector &v) : vector<T>(v) { }

	T& operator[](const size_t n)
	{
		return (T&)vector<T>::operator[](n - 1);
	}

	const T& operator[](const size_t n) const
	{
		return (T&)vector<T>::operator[](n - 1);
	}
};

struct Data
{
	void Generate(const int _n, const int _t, const int _m);
	void GenerateUncertain();
	void GenerateSimple();
	void Output(const char *pszFile) const;

	int n; // number of patients of a set N
	int t; // number of non-intersecting time slots of a set T
	int m; // number of operating rooms of a set M

	Vector2<int> N;
	Vector2<int> T;
	Vector2<int> M;

	Vector2<int> D_t; // duration of time slot [minutes]
	Vector2<int> delta_t; // maximum time slot extension value [minutes]

	Vector2<Vector2<int>> M_t; // operating rooms which can be used simultaneously in time slots
	Vector2<Vector2<int>> M_j; // set of operating rooms in which operation j can be performed
	Vector2<Vector2<Vector2<int>>> M_jt; // set of eligible operating rooms in time slot

	Vector2<Vector2<int>> N_i; // set of eligible operations for room i
	Vector2<Vector2<int>> N_t; // set of eligible operations for time slot
	Vector2<Vector2<Vector2<int>>> N_it; // set of eligible operations for room i in time slot

	Vector2<Vector2<int>> T_j; // set of time slots eligible for operation j
	Vector2<Vector2<int>> T_i; // set of time slots in which room i is available
	Vector2<Vector2<Vector2<int>>> T_ji; // set of time slots eligible for performing operation j in room i

	Vector2<int> r_j; // ready time slot
	Vector2<int> d_j; // due time slot
	Vector2<int> p_j; // time of the operation
	Vector2<int> v_j; // value of the operation
	Vector2<Vector2<int>> b_sj; // for umaxvalue
	
	Vector2<Vector2<Vector2<double>>> c_jit; // cost
	Vector2<Vector2<double>> e_it; // extension expense per time unit
};
