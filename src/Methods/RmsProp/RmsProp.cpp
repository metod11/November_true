#include "Methods/RmsProp/RmsProp.h"
//#include "Tools/Parameters.h"

//IterationData RmsProp(Function f, Vector startPoint, Vector parameters, Real grad_accuracy, int iter_lim) {
IterationData RmsProp(Function f, Vector startPoint, const StopCondition& stop_condition) {

  IterationData data;
  data.x_curr = startPoint;
  data.f_curr = f(startPoint);
  data.iter_counter = 0;

  Vector x_next;

  // Градиент и его покомпонентный квадрат
  Vector g, gd;
  // Mомент
  Vector v;

  Real f_next;
  // Параметр сохранения момента
  Real beta;
  // Параметр метода
  Real alpha;
  Real grad_accuracy;

  beta = 0.93005474;//parameters[0];
  alpha = 0.06376147;//parameters[1];
  grad_accuracy = 0.00000001;

  g = grad(f, data.x_curr, grad_accuracy);
  gd = g * g;
  v = (1 - beta) * gd;
  // Итеративная формула
  x_next = data.x_curr - alpha * (g / (sqrt(notNull(v))));
  f_next = f(x_next);

  //while (data.iter_counter < iter_lim) {
  while(!stop_condition(data)) {

    data.next(x_next, f_next);
    g = grad(f, data.x_curr, grad_accuracy);
    gd = g * g;
    v = beta * v + (1 - beta) * gd;
    // Итеративная формула
    x_next = data.x_curr - alpha * (g / (sqrt(notNull(v))));
    f_next = f(x_next);
  }

  return data;
}
