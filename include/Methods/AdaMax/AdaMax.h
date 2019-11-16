#pragma once

/* Tools */
//#include "Data.h"
#include "Tools/math.hpp"
#include "Tools/StopCondition.hpp"


/**
 *
 * @struct AdaMax
 * @brief Adaptive Moment Estimation method
 * @param AdaMax::f
 * Минимизируемая функция
 * @param AdaMax::startPoint
 * Точка старта минимизации
 * @param AdaMax::parameters
 * Параметры метода
 * @param AdaMax::grad_accuracy
 * Точность расчета градиента в методе
 * @param AdaMax::iter_limit
 * Максимальное число итераций метода
 *
*/

//IterationData AdaMax(Function f, Vector startPoint, Vector parameters, Real grad_accuracy, int iter_lim);
IterationData AdaMax(Function f, Vector start_point, const StopCondition& stop_condition = default_stop_condition);
