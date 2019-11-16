#include "Methods/AdaMax/AdaMax.h"
//#include "Tools/Parameters.h"

//IterationData AdaMax(Function f, Vector startPoint, Vector parameters, Real grad_accuracy, int iter_lim) {
IterationData AdaMax(Function f, Vector startPoint, const StopCondition& stop_condition) {

  Real grad_accuracy = 0.00000001;
  //int iter_lim = 300;

  IterationData data;
  data.x_curr = startPoint;
  data.f_curr = f(startPoint);
  data.iter_counter = 0;

  Vector x_next;

  // Вектор-градиент
  Vector g;

  Real f_next;


  // Параметры сохранения момента и бесконечной нормы v
  Real beta1, beta2;

  // Параметр метода
  Real gamma;

  // Первый момент
  Vector m;

  // Бесконечная норма градиента в пространстве Lp
  Real v = 0;

  // Норма градиента
  Real nrm;

  /* Присвоение соответствующих параметров из структуры параметров */

  /*//before: beta1 = parameters[0];beta2 = parameters[1];gamma = parameters[2];*/

  /* Авторы Adam предлагаю эти значения по умолчанию*/
  beta1 = 0.92096770;//AdaMax_beta1 = 0.92096770;	//0.9
  beta2 = 0.95563197;//AdaMax_beta2 = 0.95563197;	//0.999
  gamma = 0.18211471;//AdaMax_gamma = 0.18211471;	//0.00000001

  g = grad(f, data.x_curr, grad_accuracy);
  m = (1 - beta1) * g;
  nrm = norm(g);

  //Считаем норму по формуле v(t) = max(beta2 * v(t-1), |g(t)|)
  v = beta2 * v > nrm ? beta2 * v : nrm;
  if (!v)
    v = v + COMPARE_EPS;

  // Итеративная формула
  x_next = data.x_curr - (gamma / (1 - beta1)) * (m * (1 / v));
  f_next = f(x_next);

  //while (data.iter_counter < iter_lim) {
  while (!stop_condition(data)) {

    data.next(x_next, f_next);
    g = grad(f, data.x_curr, grad_accuracy);
    m = beta1 * m + (1 - beta1) * g;
    nrm = norm(g);

    //Считаем норму по формуле v(t) = max(beta2 * v(t-1), |g(t)|)
    v = beta2 * v > nrm ? beta2 * v : nrm;

    // Итеративная формула
    x_next = data.x_curr - (gamma / (1 - pow(beta1, data.iter_counter))) * (m * (1 / v));
    f_next = f(x_next);
  }

  return data;

}

