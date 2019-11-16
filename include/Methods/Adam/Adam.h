#pragma once

/* Tools */
//#include "Data.h"
#include "Tools/math.hpp"
#include "Tools/StopCondition.hpp"


/**
 *
 * @struct Adam
 * @brief Adaptive Moment Estimation method
 * @param Adam::f
 * Минимизируемая функция
 * @param Adam::startPoint
 * Точка старта минимизации
 * @param Adam::parameters
 * Параметры метода
 * @param Adam::grad_accuracy
 * Точность расчета градиента в методе
 * @param Adam::iter_limit
 * Максимальное число итераций метода
 *
*/

//IterationData Adam(Function f, Vector startPoint, Vector parameters, Real grad_accuracy, int iter_lim);
IterationData Adam(Function f, Vector startPoint, const StopCondition& stop_condition = default_stop_condition);
