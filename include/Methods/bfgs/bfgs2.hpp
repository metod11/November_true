#pragma once

//#include "math1.hpp"
#include "Tools/math.hpp"
#include "Tools/StopCondition.hpp"

IterationData bfgs2(Function f, Vector start_point, const StopCondition& stop_condition = default_stop_condition);

Matrix out_pr_bfgs2(Vector& x, Vector& y);

Matrix hes_upd_bfgs2(Function f, Matrix& B, Vector& x_cur, Vector& x_prv);

Real search_alpha_bfgs2(Function f, Vector& x, Vector& p, int iter_limit);

typedef std::pair<Vector, int>(*method)(Function, Vector, int);
std::pair<Vector, int> bfgs2(Function f, Vector start_point, int iter_limit = 100);
