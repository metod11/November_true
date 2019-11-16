#include "Methods/powell/powell21.hpp"
//#include "math21.hpp"
#include <cmath>
#include "Tools/math.hpp"
#include <limits>
#include <tuple>

Real f1dim21(const Real x, int ncom, Vector* pcom_p, Vector* xicom_p, Function nrfunc)
{
    Vector xt(ncom);
    Vector& pcom = *pcom_p, & xicom = *xicom_p;
    for (int j = 0; j < ncom; j++)
        xt[j] = pcom[j] + x * xicom[j];
    return nrfunc(xt);
}





std::pair<Real, Real> SWAP(Real &a, Real &b)
{
	Real dum = a; a = b; b = dum;
	return { a, b };
}

std::pair<Vector, Vector>  linmin21(Vector &p, Vector &xi, Function f)
{
	const Real TOL = 1.0e-8;
	int n = p.size();
	int ncom = n;
	Vector *pcom_p = new Vector(n);
	Vector *xicom_p = new Vector(n);
	Function nrfunc = f;
	Vector &pcom = *pcom_p, &xicom = *xicom_p;
	for (int j = 0; j<n; j++) {
		pcom[j] = p[j];
		xicom[j] = xi[j];
	}
	Real ax = 0.0;
	Real bx = 1.0;
	Real  xmin, fb, fa;
	//mnbrak
	//_______________________________________________________________
	const Real GOLD = 1.618034, GLIMIT = 100.0, TINY = 1.0e-20;
	Real ulim, u, r, q, fu;
	//Real te; // ïîä ñâàï
	fa = f1dim21(ax, ncom, pcom_p, xicom_p, nrfunc);
	fb = f1dim21(bx, ncom, pcom_p, xicom_p, nrfunc);
	
	Real cx = bx + GOLD * (bx - ax);
	Real fc = f1dim21(cx, ncom, pcom_p, xicom_p, nrfunc);
	while (fb > fc) {
		r = (bx - ax)*(fb - fc);
		q = (bx - cx)*(fb - fa);
		u = bx - ((bx - cx)*q - (bx - ax)*r) /
			(2.0*SIGN(MAX(fabs(q - r), TINY), q - r));
		ulim = bx + GLIMIT * (cx - bx);
		if ((bx - u)*(u - cx) > 0.0) {
			fu = f1dim21(u, ncom, pcom_p, xicom_p, nrfunc);
			if (fu < fc) {
				ax = bx;
				bx = u;
				fa = fb;
				fb = fu;
			//	return;
			}
			else if (fu > fb) {
				cx = u;
				fc = fu;
			//	return;
			}
			u = cx + GOLD * (cx - bx);
			fu = f1dim21(u, ncom, pcom_p, xicom_p, nrfunc);
		}
		else if ((cx - u)*(u - ulim) > 0.0) {
			fu = f1dim21(u, ncom, pcom_p, xicom_p, nrfunc);
			if (fu < fc) {
				bx = cx; cx = u; u = u + GOLD * (u - cx);
				fb = fc; fc = fu; fu = f1dim21(u, ncom, pcom_p, xicom_p, nrfunc);
			}
		}
		else if ((u - ulim)*(ulim - cx) >= 0.0) {
			u = ulim;
			fu = f1dim21(u, ncom, pcom_p, xicom_p, nrfunc);
		}
		else {
			u = cx + GOLD * (cx - bx);
			fu = f1dim21(u, ncom, pcom_p, xicom_p, nrfunc);
		}
		ax = bx; bx = cx; cx = u;
		fa = fb; fb = fc; fc = fu;
	}

	const int ITMAX = 100;
	const Real CGOLD = 0.3819660;
	const Real ZEPS = exp(-6);
	int iter;
	Real d = 0.0, etemp, fu1, fv, fx1;
	Real p1, q1, r1, tol1, tol2, u1, v, w, xm;     // íó òóò ðèàë õç êàê óïðîñòèòü
	Real e = 0.0;

	Real a = (ax < cx ? ax : cx);
	Real b = (ax > cx ? ax : cx);
	Real x = w = v = bx;
	Real fw = fv = fx1 = f1dim21(x, ncom, pcom_p, xicom_p, nrfunc);
	for (iter = 0; iter<ITMAX; iter++) {
		xm = 0.5*(a + b);
		tol2 = 2.0*(tol1 = TOL * fabs(x) + ZEPS);
		if (fabs(x - xm) <= (tol2 - 0.5*(b - a))) {
			xmin = x;
			//fret = fx1;
			break;
		}
		if (fabs(e) > tol1) {
			r1 = (x - w)*(fx1 - fv);
			q1 = (x - v)*(fx1 - fw);
			p1 = (x - v)*q1 - (x - w)*r1;
			q1 = 2.0*(q1 - r1);
			if (q1 > 0.0) p1 = -p1;
			q1 = fabs(q1);
			etemp = e;
			e = d;
			if (fabs(p1) >= fabs(0.5*q1*etemp) || p1 <= q1 * (a - x) || p1 >= q1 * (b - x))
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			else {
				d = p1 / q1;
				u1 = x + d;
				if (u1 - a < tol2 || b - u1 < tol2)
					d = SIGN(tol1, xm - x);
			}
		}
		else {
			d = CGOLD * (e = (x >= xm ? a - x : b - x));
		}
		u1 = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
		fu1 = f1dim21(u1, ncom, pcom_p, xicom_p, nrfunc);
		if (fu1 <= fx1) {
			if (u1 >= x) a = x; else b = x;

			v = w; w = x; x = u1;

			fv = fw; fw = fx1; fx1 = fu1;
		}
		else {
			if (u1 < x) a = u1; else b = u1;
			if (fu1 <= fw || w == x) {
				v = w;
				w = u1;
				fv = fw;
				fw = fu1;
			}
			else if (fu1 <= fv || v == x || v == w) {
				v = u1;
				fv = fu1;
			}
		}
	}
	//fret = fx1;
	xmin = x;
	//return fx;
	//_____________________________________________________________
	for (int j = 0; j<n; j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	return { p, xi };
	delete xicom_p;
	delete pcom_p;
}




IterationData powell21(Function func, Vector p, const StopCondition& stop_condition)
{
    ///////////////
    IterationData iter_data;
    Real MIN_VAL = 100000; // ìèíèìóì ôóíêöèè
    int NDIM = (int)p.size(); // îïðåëÿåì ðàçìåðíîñòü
    Vector x_curr = p;
    //const Real ftol = 1.0e-6;
    Vector p_d = p;
    iter_data.x_curr = p_d;
    iter_data.f_curr = func(p_d);
    iter_data.iter_counter = 0;
    iter_data.method_title = "powell21";
    Real** xi = new Real * [NDIM];
    for (int count = 0; count < NDIM; count++)
        xi[count] = new Real[NDIM];

    for (int i = 0; i < NDIM; i++)
        for (int j = 0; j < NDIM; j++)
            xi[i][j] = (i == j ? 1.0 : 0.0);

    Real fret;
    ///////////////////////////


	//const int ITMAX = 200;
	//int iter;
	//const Real TINY = 1.0e-25;
	Real  fp, fptt, t;
	int n = p.size();
	Vector pt(n), ptt(n), xit(n);
	fret = func(p);
	for (int j = 0; j<n; j++) pt[j] = p[j];
    do {
    //for (iter = 0; iter<200; iter++) {
		//5
        fp = fret;
		int ibig = 0;
		Real del = 0.0;

		for (int i = 0; i<n; i++) {
			for (int j = 0; j<n; j++) xit[j] = xi[j][i];
			fptt = fret;
			std::pair<Vector, Vector> qwe = linmin21(p, xit, func);

            p = qwe.first;
			xit = qwe.second;
			if (fptt - fret > del) {
				del = fptt - fret;
				ibig = i + 1;
			}
		}

        // 2 новых if
		//if (2.0*(fp - fret) <= ftol * (fabs(fp) + fabs(fret)) + TINY) {
		//	return iter_data;
		//}
		//if (iter == ITMAX) return iter_data;
        ///Та же часть что и в 5
		for (int j = 0; j<n; j++) {
			ptt[j] = 2.0*p[j] - pt[j];
			xit[j] = p[j] - pt[j];
			pt[j] = p[j];
		}
		fptt = func(ptt);
        iter_data.next(ptt, fptt);
		if (fptt < fp) {
			t = 2.0*(fp - 2.0*fret + fptt)*SQR(fp - fret - del) - del * SQR(fp - fptt);
			if (t < 0.0) {
				std::pair<Vector, Vector> qwe = linmin21(p, xit, func);
				p = qwe.first;
				xit = qwe.second;
				for (int j = 0; j<n; j++) {
					xi[j][ibig - 1] = xi[j][n - 1];
					xi[j][n - 1] = xit[j];
				}
			}
		}

	} while (!stop_condition(iter_data)); // Проверяем условие останова

    /////////////////////////////////////
    //x_curr = std::get<0>(ans);
    //x_cur = ans.first;
    //Real iter = ans.second;
    //int iter = std::get<1>(ans) = 1;
    //fret = ans.third;
    //fret = std::get<2>(ans);
    //x_cur = p;
    if (fret < MIN_VAL) { MIN_VAL = fret; }
    return iter_data;
}

