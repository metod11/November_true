#pragma once

/* Tools */
//#include "Data3.h"
#include "Tools/math.hpp"
#include "Tools/StopCondition.hpp"

/**
*
* @struct Adadelta
* @brief Adaptive Moment Estimation method
* @param Adadelta::f
* Минимизируемая функция
* @param Adadelta::startPoint
* Точка старта минимизации
* @param Adadelta::parameters
* Параметры метода
* @param Adadelta::grad_accuracy
* Точность расчета градиента в методе
* @param Adadelta::iter_limit
* Максимальное число итераций метода
*
*/


IterationData Adadelta(Function f, Vector startPoint, const StopCondition& stop_condition = default_stop_condition);
