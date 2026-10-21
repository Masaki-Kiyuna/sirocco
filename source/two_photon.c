/**
 * @file two_photon.c
 * @brief H I two-photon continuum process and Cloudy-C17 shape adapter.
 */

#include <math.h>
#include <stdio.h>

#include "atomic.h"
#include "sirocco.h"

typedef struct
{
  double a[8], b[8], c[8], d[8];
  double i0[9], i1[9];
  double norm;
  int initialized;
}
TPCloudyShape;

static TPCloudyShape h1_shape;

static void
tp_half_integrals (s, x, i0, i1)
     const TPCloudyShape *s;
     double x;
     double *i0, *i1;
{
  const double h = 1.0 / 16.0;
  int k;
  double t, v0, v1;

  if (x <= 0.0)
  {
    *i0 = *i1 = 0.0;
    return;
  }
  if (x >= 0.5)
  {
    *i0 = s->i0[8];
    *i1 = s->i1[8];
    return;
  }
  k = (int) (16.0 * x);
  t = x - k * h;
  v0 = t * (s->a[k] + t * (s->b[k] / 2.0 + t * (s->c[k] / 3.0 + t * s->d[k] / 4.0)));
  v1 = t * t * (s->a[k] / 2.0 + t * (s->b[k] / 3.0 + t * (s->c[k] / 4.0 + t * s->d[k] / 5.0)));
  *i0 = s->i0[k] + v0;
  *i1 = s->i1[k] + k * h * v0 + v1;
}

static void
tp_shape_init ()
{
  static const double value[9] = { 0.0, 2.03239, 3.15792, 3.84452, 4.28435, 4.56958, 4.74855, 4.84732, 4.87892 };
  const double h = 1.0 / 16.0;
  double second[9] = { 0 }, cp[9] = { 0 }, dp[9] = { 0 };
  double rhs, denominator, z0, z1;
  int k;

  if (h1_shape.initialized)
    return;

  for (k = 1; k < 8; k++)
  {
    rhs = 6.0 * (value[k + 1] - 2.0 * value[k] + value[k - 1]) / (h * h);
    denominator = 4.0 - cp[k - 1];
    cp[k] = (k < 7 ? 1.0 / denominator : 0.0);
    dp[k] = (rhs - dp[k - 1]) / denominator;
  }
  for (k = 7; k > 0; k--)
    second[k] = dp[k] - cp[k] * second[k + 1];

  h1_shape.i0[0] = h1_shape.i1[0] = 0.0;
  for (k = 0; k < 8; k++)
  {
    h1_shape.a[k] = value[k];
    h1_shape.b[k] = (value[k + 1] - value[k]) / h - h * (2.0 * second[k] + second[k + 1]) / 6.0;
    h1_shape.c[k] = second[k] / 2.0;
    h1_shape.d[k] = (second[k + 1] - second[k]) / (6.0 * h);
    z0 = h * (h1_shape.a[k] + h * (h1_shape.b[k] / 2.0 + h * (h1_shape.c[k] / 3.0 + h * h1_shape.d[k] / 4.0)));
    z1 = h * h * (h1_shape.a[k] / 2.0 + h * (h1_shape.b[k] / 3.0 + h * (h1_shape.c[k] / 4.0 + h * h1_shape.d[k] / 5.0)));
    h1_shape.i0[k + 1] = h1_shape.i0[k] + z0;
    h1_shape.i1[k + 1] = h1_shape.i1[k] + k * h * z0 + z1;
  }
  h1_shape.norm = 2.0 * h1_shape.i0[8];
  h1_shape.initialized = TRUE;
}

static double
tp_cdf (y)
     double y;
{
  double a, b;
  tp_shape_init ();
  if (y <= 0.0)
    return 0.0;
  if (y >= 1.0)
    return 1.0;
  if (y <= 0.5)
  {
    tp_half_integrals (&h1_shape, y, &a, &b);
    return 2.0 * b / h1_shape.norm;
  }
  tp_half_integrals (&h1_shape, 1.0 - y, &a, &b);
  return 1.0 - 2.0 * (a - b) / h1_shape.norm;
}

static double
tp_inverse_cdf (target, lo, hi)
     double target;
     double lo, hi;
{
  double mid;
  int k;

  for (k = 0; k < 60; k++)
  {
    mid = lo + 0.5 * (hi - lo);
    if (mid == lo || mid == hi)
      break;
    if (tp_cdf (mid) < target)
      lo = mid;
    else
      hi = mid;
  }
  return lo + 0.5 * (hi - lo);
}

int
two_photon_process_for_upper (upper)
     int upper;
{
  if (!geo.h1_two_photon_emission || upper < 0 || upper >= nlevels)
    return -1;
  return two_photon_from_config[upper];
}

double
two_photon_rate (iprocess)
     int iprocess;
{
  if (!geo.h1_two_photon_emission || iprocess < 0 || iprocess >= n_two_photon)
    return 0.0;
  return two_photon[iprocess].a_spont * two_photon[iprocess].collapsed_frac;
}

double
two_photon_nu0 (iprocess)
     int iprocess;
{
  if (iprocess < 0 || iprocess >= n_two_photon)
    return 0.0;
  return (xconfig[two_photon[iprocess].nconfigu].ex - xconfig[two_photon[iprocess].nconfigl].ex) / PLANCK;
}

double
two_photon_energy_weight (upper)
     int upper;
{
  int iprocess = two_photon_process_for_upper (upper);
  if (iprocess < 0)
    return 0.0;
  return two_photon_rate (iprocess) * (xconfig[upper].ex - xconfig[two_photon[iprocess].nconfigl].ex);
}

double
two_photon_band_fraction (iprocess, freq_min, freq_max)
     int iprocess;
     double freq_min, freq_max;
{
  double nu0, answer;
  nu0 = two_photon_nu0 (iprocess);
  if (!(nu0 > 0.0) || !(freq_max > freq_min))
    return 0.0;
  answer = tp_cdf (freq_max / nu0) - tp_cdf (freq_min / nu0);
  if (answer < 0.0)
    answer = 0.0;
  if (answer > 1.0)
    answer = 1.0;
  return answer;
}

double
two_photon_sample_freq (iprocess, freq_min, freq_max)
     int iprocess;
     double freq_min, freq_max;
{
  double nu0, lo, hi, cdf_lo, cdf_hi, target, u, nu;

  nu0 = two_photon_nu0 (iprocess);
  if (!(nu0 > 0.0) || !(freq_max > freq_min))
    return NAN;
  lo = fmax (0.0, freq_min / nu0);
  hi = fmin (1.0, freq_max / nu0);
  if (!(hi > lo))
    return NAN;
  cdf_lo = tp_cdf (lo);
  cdf_hi = tp_cdf (hi);
  if (!(cdf_hi > cdf_lo))
    return NAN;
  u = random_number (0.0, 1.0);
  if (u <= 0.0)
    u = nextafter (0.0, 1.0);
  if (u >= 1.0)
    u = nextafter (1.0, 0.0);
  target = cdf_lo + u * (cdf_hi - cdf_lo);
  nu = nu0 * tp_inverse_cdf (target, lo, hi);
  if (nu <= 0.0)
    nu = nextafter (0.0, nu0);
  if (nu >= nu0)
    nu = nextafter (nu0, 0.0);
  return nu;
}

int
two_photon_validate_shape ()
{
  double upper_third_fraction, u, y, previous_y = 0.0, max_roundtrip_error = 0.0, error;
  int k;
  tp_shape_init ();
  upper_third_fraction = tp_cdf (1.0) - tp_cdf (1.0 / 3.0);
  for (k = 1; k < 1000; k++)
  {
    u = k / 1000.0;
    y = tp_inverse_cdf (u, 0.0, 1.0);
    error = fabs (tp_cdf (y) - u);
    if (error > max_roundtrip_error)
      max_roundtrip_error = error;
    if (!(y > previous_y))
      return -1;
    previous_y = y;
  }
  Log ("TwoPhoton shape: Cloudy C17.03 natural spline norm=%.12g energy-fraction(y>1/3)=%.10f inverse-CDF max-error=%.3e\n",
       h1_shape.norm, upper_third_fraction, max_roundtrip_error);
  return (fabs (upper_third_fraction - 0.8839678243) < 1.e-8 && max_roundtrip_error < 1.e-12 ? 0 : -1);
}
