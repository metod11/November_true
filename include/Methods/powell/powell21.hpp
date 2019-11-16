#pragma once

//#include "math21.hpp"
#include <cmath>
#include <math.h>
#include <limits>
#include "Tools/math.hpp"
#include "Tools/StopCondition.hpp"

IterationData powell21(Function f, Vector startingPoint, const StopCondition& stop_condition = default_stop_condition);

std::pair<Vector, Vector> linmin(Vector &p, Vector &xi, Function f);

Real f1dim(const Real x, int ncom, Vector *pcom_p, Vector *xicom_p, Function f);

//void powell1(Vector &p, Matrix &xi, const double ftol, int &iter,
//	double &fret, double func(Vector &));

//std::tuple<Vector, int, Real> powell1(Vector &p, Real **xi, const Real ftol, &fret, Function f);

std::pair<Real, Real> SWAP(Real &a, Real &b);
