#pragma once

/* Tools */
//#include "Data.h"
#include "Tools/math.hpp"
#include "Tools/StopCondition.hpp"

/**
 *
 * @struct AMSGrad
 * @brief AMSGrad that uses the maximum of past squared gradients v_t
 * rather than the exponential average to update the parameters.
 * v_t is defined the same as in Adam
 * @param AMSGrad::startPoint
 * Точка старта минимизации
 * @param AMSGrad::grad_accuracy
 * Точность расчета градиента в методе
 * @param AMSGrad::iter_limit
 * Максимальное число итераций метода
 *
*/


//IterationData AMSGrad(Function f, Vector startPoint, Vector parameters, Real grad_accuracy, int iter_lim);
IterationData AMSGrad(Function f, Vector start_point, const StopCondition& stop_condition = default_stop_condition);
