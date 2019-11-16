/* Adaptive Gradient */
#include "Methods/Adagrad/Adagrad.h"

//#include "Parameters.hpp"

IterationData Adagrad(Function f, Vector startPoint, const StopCondition& stop_condition) {

  IterationData data;
  data.x_curr = startPoint;
  data.f_curr = f(startPoint);
  data.iter_counter = 0;
  data.method_title = "Adardad";


  Vector x_next;
  // Вектор-градиент и его покомпонентный квадрат
  Vector g, gd;
  // Момент
  Vector v;

  Real f_next;

  // Параметр метода
  Real beta;
  Real grad_accuracy;
  //beta = parameters[0];

  beta = 1.0;  //0.9
  grad_accuracy = 0.00000001;

  g = grad(f, data.x_curr, grad_accuracy);
  gd = g * g;
  v = gd;
  // Итеративная формула
  x_next = data.x_curr - beta * (g / (sqrt(notNull(v))));
  f_next = f(x_next);

  do{
    data.next(x_next, f_next);
    g = grad(f, data.x_curr, grad_accuracy);
    gd = g * g;
    v = v + gd;
    // Итеративная формула
    x_next = data.x_curr - beta * (g / sqrt(notNull(v)));
    f_next = f(x_next);
  }while(!stop_condition(data));

  return data;

}
