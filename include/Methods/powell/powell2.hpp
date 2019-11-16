#pragma once


#include <fstream>
#include <complex>
#include <iostream>
#include <cmath>
#include "Tools/math.hpp"
#include "Tools/StopCondition.hpp"

IterationData powell2(Function func, Vector p, const StopCondition& stop_condition = default_stop_condition);

Real brent2(const Real ax, const Real bx, const Real cx, Real f(const Real, int, Vector *, Vector *, Function),
	const Real tol, int ncom, Vector *pcom_p, Vector *xicom_p, Function);

Real f1dim2(const Real x, int ncom, Vector *pcom_p, Vector *xicom_p, Function func);

std::pair<Vector, Vector> linmin2(Vector &pInit, Vector &xiInit, Function func);

std::pair<Vector, Vector> mnbrak2(Real axInit, Real bxInit, Real func(const Real, int, Vector *, Vector *, Function), int ncom, Vector *pcom_p, Vector *xicom_p, Function);
// func - указатель на целевую функцию
// p - начальное приближение
// stop_condition - критерий остановы
// Результат работы метода будет лежать в структуре данных о последней итерации
