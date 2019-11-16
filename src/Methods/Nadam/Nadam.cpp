#include "Methods/Nadam/Nadam.h"
//#include "Tools/Parameters.h"

IterationData Nadam(Function f, Vector startPoint, const StopCondition& stop_condition) {
  Real grad_accuracy = 0.00000001;
  //int iter_lim = 3000;
  IterationData data;
  data.x_curr = startPoint;
  data.f_curr = f(startPoint);
  data.iter_counter = 0;

  Vector x_next;
  // Вектор-градиент и его покомпонентный квадрат
  Vector g, gd;

  Real f_next;

  // Параметры сохранения моментов
  Real beta1, beta2;

  // Параметр метода
  Real gamma;

  // Первый и второй моменты
  Vector m, m_t, v, v_t;

  // Временная переменная для расчетов
  Vector tmp;

  /* Присвоение соответствующих параметров из структуры параметров */

  beta1 = 0.88441989;
  beta2 = 0.88441989;
  gamma = 0.08381000;


  g = grad(f, data.x_curr, grad_accuracy);
  gd = g * g;
  m = (1 - beta1) * g;
  v = (1 - beta2) * gd;
  m_t = m * (1 / (1 - beta1));
  v_t = v * (1 / (1 - beta2));

  tmp = beta1 * m_t + ((1 - beta1) * g) * (1 / ( 1 - beta1 ) );

  // Итеративная формула
  x_next = data.x_curr - gamma * (tmp / (sqrt(notNull(v_t))));
  f_next = f(x_next);

  while (!stop_condition(data)) {
    data.next(x_next, f_next);
    g = grad(f, data.x_curr, grad_accuracy);
    gd = g * g;
    m = beta1 * m + (1 - beta1) * g;
    v = beta2 * v + (1 - beta2) * gd;

    m_t = m * (1 / (1 - pow(beta1, data.iter_counter)));
    v_t = v * (1 / (1 - pow(beta2, data.iter_counter)));

    tmp = beta1 * m_t + ((1 - beta1) * g) * (1 / ( 1 - pow(beta1, data.iter_counter) ) );

    // Итеративная формула
    x_next = data.x_curr - gamma * (tmp / (sqrt(notNull(v_t))));
    f_next = f(x_next);
  }

  return data;
}
