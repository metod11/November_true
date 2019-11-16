#pragma once

/* Tools */
//#include "Data3.h"
#include "Tools/math.hpp"
#include "Tools/StopCondition.hpp"

/**
 *
 * @struct Nadam
 * @brief Nesterov-accelerated Adaptive Moment Estimation (Nadam)
 * @param Nadam::f
 * Минимизируемая функция
 * @param Nadam::startPoint
 * Точка старта минимизации
 * @param Nadam::parameters
 * Параметры метода
 * @param Nadam::grad_accuracy
 * Точность расчета градиента в методе
 * @param Nadam::iter_limit
 * Максимальное число итераций метода
 *
*/


IterationData Nadam(Function f, Vector startPoint, const StopCondition& stop_condition = default_stop_condition);
