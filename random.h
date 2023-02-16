#ifndef _RANDOM_
#define _RANDOM_

/*
 * Module : random.c   
 */
long rnd32 (long);
double uniform (double, double, long *);
double negexp (double, long *);
int poisson (double, long *);
int geometric0 (double, long *);
int geometric1 (double, long *);
int geometric_trunc1 (double, int, long *);
int trunc_exp (double, long, long *);
double eval_gauss_sample (long *, double, double);
int randint (int, int, long *);
double weibull (double a, double b, long *);
double iperexp (double alpha, double mu_1, double mu_2, long *);
double pareto (double a, long *);
double erlang (double a, double M, long *);
double ipererl (double alpha, double m_1, double a_1, double m_2, double a_2,
		long *);


#endif
