/*
 * **  Module  : random.c **  Version : 1.1 **  Date    : March 13, 1993 
 */

#include <stdio.h>
#include <math.h>

#define MODULE  2147483647
#define A       16807
#define LASTXN  127773
#define UPTOMOD -2836
#define RATIO   0.46566128e-9	/*
				 * 1/MODULE 
				 */


/*
 * **  Function  : long rnd32(long seed) **  Return    : the updated value
 * of 'seed'. **  Remarks   : congruential generator of pseudorandom
 * sequences of numbers **              uniformly distributed between 1 and
 * 2147483646, using the **              congruential relation: Xn+1 = 16807 
 * * Xn  mod  2147483647 . **              The input to the routine is the
 * integer Xn, while the returned **              integer is the number
 * Xn+1. 
 */
long
rnd32 (long seed)
{
  long times, rest, prod1, prod2;

  times = seed / LASTXN;
  rest = seed - times * LASTXN;
  prod1 = times * UPTOMOD;
  prod2 = rest * A;
  seed = prod1 + prod2;
  if (seed < 0)
    seed = seed + MODULE;
  return seed;
}


/*
 * Function : int randint(int min,int max,long *seed) **  Return   : an
 * integer number uniformly distributed between MIN and **             MAX.
 * **  Remarks  : The value of '*seed' is changed 
 */

int
randint (int min, int max, long *seed)
{
  double temp;
  *seed = rnd32 (*seed);
  temp = (double) (*seed) - 1.0;
  temp = temp * (max - min + 1.0) / (MODULE - 1);
  if (((int) temp + min) > max)
    printf ("Errore in rnd_int fuori range!!! temp = %f \n", temp);
  return ((int) temp + min);
}

/*
 * **  Function  : double uniform(double a, double b, long *seed) **  Return
 *   : a value uniformly distributed between 'a' and 'b' **  Remarks   : the
 * value of '*seed' is changed. 
 */
double
uniform (double a, double b, long *seed)
{
  double u;
  *seed = rnd32 (*seed);
  u = (*seed) * RATIO;
  u = a + u * (b - a);
  return u;
}


/*
 * **  Function  : double negexp(double mean, long *seed) **  Return    : a
 * value exponentially distributed with mean 'mean'. **  Remarks   : the
 * value of '*seed' is changed. 
 */
double
negexp (double mean, long *seed)
{
  double u;
  *seed = rnd32 (*seed);
  u = (*seed) * RATIO;
  return (-mean * log (u));
}


/*
 * **  Function  : int poisson(double alpha,long *seed) **  Return    : the
 * number of users arrived, according to a poisson process **
 * with rate 'alpha' user/slot, within a slot. **  Remarks   : the value of
 * '*seed' is changed. 
 */
int
poisson (double alpha, long *seed)
{
  int n = 0;
  double pn, lim;
  double prob;

  lim = pn = exp (-alpha);
  prob = uniform (0.0, 1.0, seed);
  while (prob > lim)
    {
      n++;
      pn *= alpha / n;
      lim += pn;
    }
  return n;
}


/*
 * **  Function  : int geometric0(double mean,long *seed) **  Return    : a
 * random value distributed geometrically with average 'mean', **
 *  starting from 0 (0 w.p. 1-p, 1 w.p. p(1-p), etc.). **  Remarks   : the
 * value of '*seed' is changed. 
 */
int
geometric0 (double mean, long *seed)
{
  int n;
  double prob, status;

  status = mean / (1.0 + mean);	/*
				 * E[X] = p/(1-p) -> p =
				 * E[X]/(1+E[X])  
				 */
  prob = uniform (0.0, 1.0, seed);	/*
					 * 1-p = prob. di avere n = 0
					 *   
					 */
  n = (int) floor (log (1 - prob) / log (status));
  return n;
}


/*
 * **  Function: int geometric1(double mean,long *seed) **  Return  : a
 * random value distributed geometrically with average 'mean', **
 * starting from 1 (1 w.p. 1-p, 2 w.p. p(1-p), etc.). **  Remarks : the value 
 * of '*seed' is changed. 
 */
int
geometric1 (double mean, long *seed)
{
  int n;
  double prob, status;

  status = (mean - 1) / mean;	/*
				 * E[X] = 1/(1-p) -> p = (E[X]-1)/E[X]  
				 */
  prob = uniform (0.0, 1.0, seed);	/*
					 * 1-p = prob. di avere n = 1
					 *   
					 */
  n = 1 + (int) floor (log (1 - prob) / log (status));
  return n;
}


/*
 * **  Function  : int geometric_trunc1(double mean,int max_len,long *seed)
 * **  Return    : a random value distributed geometrically with average
 * 'mean', **              starting from 1. **              The distribution
 * is truncated at the value 'max_len'. **  Remarks   : the value of '*seed'
 * is changed. 
 */
int
geometric_trunc1 (double mean, int max_len, long *seed)
{
  /*
   * These function returns a number distributed quasi-geometrically with
   * 
   */
  /*
   * average mean and maximum value 'max_len'.
   * 
   */
  /*
   * There are some problems with the calculation. Here we explain the way
   * 
   */
  /*
   * the numbers are calculated.
   * 
   */
  /*
   * The mean value of the random variable is:
   * 
   */
  /*
   */
  /*
   * Sum(i*p^(i-1),i=1..N)                              
   */
  /*
   * E[x] = --------------------- = m                          
   */
  /*
   * Sum(p^(i-1),i=1..N)                               
   */
  /*
   * i.e.
   * 
   */
  /*
   * p^N ( Np - N - 1) + 1                              
   */
  /*
   * m = ---------------------          (1)                 
   */
  /*
   * (1-p)(1-p^N)                                   
   */
  /*
   */
  /*
   * where p is the transition probability in the Markov chain of the
   * model. 
   */
  /*
   */
  /*
   * We need the value of p as a function of m and N. The only solution is
   * 
   */
  /*
   * to solve the equation (1) in the variable p using the Newton method,
   * 
   */
  /*
   * i.e.
   * 
   */
  /*
   * p' = p - f(p)/f'(p)                                          
   */
  /*
   * being p' the value of p at the step i+1, p the value at the step i,
   * 
   */
  /*
   * f(p) is (1) and f'(p) is df(p)/dp.
   * 
   */
  /*
   * In our calculations, we use:
   * 
   */
  /*
   */
  /*
   * f(p)  = p^N * ((m-N)p + N - m + 1) - mp + m -1                      
   */
  /*
   * f'(p) = (m-N) p^N + N p^(N-1)((m-N)p + N - m + 1) - m               
   */
  /*
   */
  /*
   * and the value  p = (m-1)/m  as starting point. This is the value of
   * 
   */
  /*
   * p when N tends to infinity.
   * 
   */
  /*
   */
  /*
   * This value of p is used to find the number n to be returned.  A random 
   */
  /*
   * variable q uniformly distributed in (0,1) is extracted, so if
   * 
   */
  /*
   */
  /*
   * sum(p^(i-1),i=1..n)    1 - p^n                           
   */
  /*
   * q = --------------------- = -------                           
   */
  /*
   * sum(p^(i-1),i=1..N)    1 - p^N                           
   */
  /*
   */
  /*
   * we found that
   * 
   */
  /*
   */
  /*
   * |~  log(p^N * q - q - 1) ~|                               
   */
  /*
   * n = |   --------------------  |                               
   */
  /*
   * |         log(p)          |                               
   */
  /*
   */
  /*
   * In order to avoid large computations, the previous values of 'mean'
   * 
   */
  /*
   * and 'max_len' are recorded, so if the function is called twice or more 
   */
  /*
   * times consecutively with the same parameters, the previously computed
   * 
   */
  /*
   * value of p can be used.
   * 
   */
  /*
   */
  /*
   * In the code, there is the corrispondence:
   * 
   */
  /*
   * p     -> status                                             
   */
  /*
   * m     -> mean                                               
   */
  /*
   * N     -> max_len                                            
   */
  /*
   * q     -> prob                                               
   */
  /*
   * f(p)  -> f_p                                                
   */
  /*
   * f'(p) -> df_p                                               
   */
  /*
   * between the symbols used in this comment and the variables names.
   * 
   */

  int n;
  double prob, f_p, df_p, temp_status, temp_res, len;
  static double status = 0.0, old_mean = 0.0, status_N = 0.0;
  static int old_max = 0;

  if (mean >= (double) max_len)
    {
      printf ("Error Calling Geometric_Trunc1()\n");
      return 1;
    }
  if (fabs (old_mean - mean) > 1e-5 || old_max != max_len)
    {
      len = (double) max_len;
      temp_status = (mean - 1) / mean;
      do
	{
	  status = temp_status;
	  status_N = pow (status, len);
	  temp_res = (mean - len) * status + len - mean + 1;
	  f_p = status_N * temp_res - mean * status + mean - 1;
	  df_p =
	    (mean - len) * status_N +
	    len * status_N * temp_res / status - mean;
	  temp_status = status - f_p / df_p;
	}
      while (fabs (temp_status - status) > 1e-9);
      status = temp_status;
      status_N = pow (status, len);
      old_mean = mean;
      old_max = max_len;
    }

  prob = uniform (0.0, 1.0, seed);
  n = 1 + (int) floor (log (1 - prob + prob * status_N) / log (status));
  return n;
}


/*
 * **  Function  : int trunc_exp(double mean,long length,long *seed) **
 * Return    : a value extracted from a truncated exponential density **
 *          function. **  Remarks   : mean and length are expressed in
 * bytes. **              The value of '*seed' is changed. 
 */
int
trunc_exp (double mean, long length, long *seed)
{
  double len, prob;

  *seed = rnd32 (*seed);
  prob = (*seed) * RATIO;
  /*
   * len =  - 8*mean*(log(*seed)-21.4875626); 
   */
  len = -8 * mean * log (prob);
  len = (len > length * 8.) ? length : len / 8.;
  return ((int) len == 0 ? 1 : (int) len);
}

/*
 * **  Function  : double eval_gauss_sample(long *seed, double mean, **
 *                                 double variance) **  Return    : a value
 * extracted from a gaussian density function. **  Remarks   : the sum of 12
 * independent variables uniformly distributed **              between -0.5
 * and 0.5 is used to obtain a Normally **              distributed variable. 
 */
double
eval_gauss_sample (long *seed, double mean, double variance)
{
  double acc_x = 0.0;
  int i = 0;

  for (; i < 12; i++)
    acc_x += uniform (-0.5, 0.5, seed);

  return (acc_x * sqrt (variance) + mean);
}

/*        
**  Function: double weibull(long *seed, double a, double b) 
**  Return  : a random variable with a Weibull p.d.f.:
**			 
**		      b * x^(b-1)          
**	      f(x) = ------------ * exp(-(x/a)^b)
**		         a^b		 
**
**  "a" is the "scale" parameter;
**  "b" is the "shape" parameter.
**
**  If 	b = 3.602, the Weibull distr. is very close to a normal;
**  for b > 3.602, it has a long left tail; 
**  for b < 3.602, it has a long right tail.
**  
**  If 	b <= 1 the Weibull p.d.f. is "L-shaped", while
**  if	b >  1 it is "bell-shaped".
**
**  For large values of b the Weibull p.d.f. has a sharp peak
**  at the mode.                       
**  The mean of the r.v. is: 
**
**		 a
**  		--- * Gamma(1/b),     
**  		 b
**		     	  
**  and the variance is: 
**
** 		 a^2	
**		----- * (2b*Gamma(2/b)-(Gamma(1/b))^2).
**		 b^2
**
**  (Gamma(x) is the Gamma Eulero's function).
** 
**  Remarks : The value of '*seed' is changed
*/
double
weibull (double a, double b, long *seed)
{
  double u;
  u = uniform (0.0, 1.0, seed);
  return (a * pow (log (1 / u), 1 / b));
}



/*        
**  Function: double iperexp(long *seed, double alpha, double mu_1, double mu_2) 
**  Return  : a random variable with a 2th order 
**  hyperexponential p.d.f.:
**
**	f(x) = alpha*mu_1*exp(-mu_1*x)+(1-alpha)*mu_2*exp(-mu_2*x)
**
**  Remarks : The value of '*seed' is changed
*/
double
iperexp (double alpha, double mu_1, double mu_2, long *seed)
{
  double u;
  u = uniform (0.0, 1.0, seed);
  if (u <= alpha)
    return (negexp (mu_1, seed));	/* ?? DUBBIO ?? avendo gia' scelto, ci 
					   vuole ancora la moltiplicazione
					   per alpha?                        */
  else				/* u > alpha */
    return (negexp (mu_2, seed));	/* Stesso DUBBIO di sopra!              */
}


/*        
**  Function: double pareto(long *seed, double a) 
**  Return  : a random variable with a Pareto p.d.f.:
**
**	f(x) =  a * x^(-(a+1))
**
**  The mean of the r.v. is:
**
**	   a	
**	-------		for a > 1
**	 a - 1
**
** while its variance is:
**
**		  a
**	-----------------------		for a > 2
**	 ((a - 1)^2) * (a - 2)	
**
** ???? (chiarire i dubbi: per a<1  la media DIVERGE???
** Idem per la varianza!!!)
**
**  Remarks : The value of '*seed' is changed
**
*/


double
pareto (double a, long *seed)
{
  double u;
  u = uniform (0.0, 1.0, seed);
  return (pow (1 / u, 1 / a));
}

/* 
** Return : a random variable with Erlang p.d.f.:
**
**	f(x) = x^(M-1)*exp(-x/a)/((a^M)*(M-1)!)
** Remarks : The value of '*seed' is changed
**
*/

double
erlang (double a, double M, long *seed)
{
  double u;
  double p;
  int i;
  for (p = 1, i = 0; i < M; i++)
    {
      u = uniform (0.0, 1.0, seed);
      p *= u;
    }
  return (-a * log (p));
}


/*         
**  Function: double ipererl(long *seed, double alpha, double m_1, double a_1,double m_2, double a_2)  
**  Return  : a random variable with a 2th order  
**  hypererlang p.d.f.: 
** 
**	f(x) = 
** alpha*((x^(m_1-1)*exp(-x/a_1))/((a_1^m_1)*(m_1-1)!))+(1-alpha)*((x^(m_2-1)*exp(-x/a_2))/((a_2^m_2)*(m_2-1)!))
**  Remarks : The value of '*seed' is changed 
*/
double
ipererl (double alpha, double m_1, double a_1, double m_2, double a_2,
	 long *seed)
{
  double u;
  u = uniform (0.0, 1.0, seed);
  if (u <= alpha)
    return (erlang (a_1, m_1, seed));	/* ?? DUBBIO ?? avendo gia' scelto, ci  
					   vuole ancora la moltiplicazione 
					   per alpha?                     */
  else				/* u > alpha */
    return (erlang (a_2, m_2, seed));	/* Stesso DUBBIO di sopra!              */
}
