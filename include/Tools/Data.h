#pragma once

#include <vector>
#include <iostream>
#include <cassert>
#include <iomanip>
#include <cmath>

using namespace std;

typedef long double Real;
typedef std::vector<Real> Vector;
typedef Real(*Function)(const Vector & x);

//const Real COMPARE_EPS = 0.0000000000000001L;


/**
 *
 * @struct IterationData
 * @brief Структура данных итерации метода
 * @var IterationData::x_prev
 * Предыдущий рассчитанный минимум
 * @var IterationData::f_prev
 * Значение функции в предыдущей рассчитанной точке
 * @var IterationData::x_curr
 * Текущий рассчитанный минимум
 * @var IterationData::f_curr
 * Значение функции в текущей точке
 * @var IterationData::iter_counter
 * Число отработанных итераций на текущий момент
*/

/*
typedef struct IterationData {
    Vector x_prev;
	Real   f_prev;
	Vector x_curr;
	Real   f_curr;
	int    iter_counter;

	void next(const Vector& x_next, const Real f_next);

	IterationData();
} ITERDATA;

typedef ITERDATA(*Method)(Function f, Vector startPoint, Vector parameters, Real grad_accuracy, int iter_lim);
*/