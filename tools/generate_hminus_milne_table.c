/* Generate the H- spontaneous radiative-attachment coefficient from the
 * McLaughlin photodetachment table using the ground-channel Milne relation.
 * Production SIROCCO reads the generated CSV; Python is not required. */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHI_EV 0.7542
#define EMAX_EV 10.95
#define EV_ERG 1.602176634e-12
#define H_PLANCK 6.62607015e-27
#define C_LIGHT 2.99792458e10
#define M_ELECTRON 9.1093837015e-28
#define K_BOLTZMANN 1.380649e-16

typedef struct
{
  size_t n;
  double *energy_ev;
  double *sigma;
} CrossSection;

static const double gl_x[8] = {
  -0.9602898564975363, -0.7966664774136267,
  -0.5255324099163290, -0.1834346424956498,
   0.1834346424956498,  0.5255324099163290,
   0.7966664774136267,  0.9602898564975363
};
static const double gl_w[8] = {
  0.1012285362903763, 0.2223810344533745,
  0.3137066458778873, 0.3626837833783620,
  0.3626837833783620, 0.3137066458778873,
  0.2223810344533745, 0.1012285362903763
};

static int
read_cross_section (const char *path, CrossSection *table)
{
  FILE *stream = fopen (path, "r");
  char line[512];
  size_t capacity = 4096;
  double energy, sigma;

  if (stream == NULL)
  {
    fprintf (stderr, "cannot open %s: %s\n", path, strerror (errno));
    return -1;
  }
  table->energy_ev = calloc (capacity, sizeof (*table->energy_ev));
  table->sigma = calloc (capacity, sizeof (*table->sigma));
  if (table->energy_ev == NULL || table->sigma == NULL)
    return -1;
  while (fgets (line, sizeof (line), stream) != NULL)
  {
    if (sscanf (line, "%lf,%lf", &energy, &sigma) != 2)
      continue;
    if (!(energy > CHI_EV) || !(sigma >= 0.0)
        || (table->n > 0 && energy <= table->energy_ev[table->n - 1]))
    {
      fprintf (stderr, "invalid cross-section row %zu\n", table->n + 2);
      return -1;
    }
    if (table->n == capacity)
    {
      double *new_energy, *new_sigma;
      capacity *= 2;
      new_energy = realloc (table->energy_ev, capacity * sizeof (*table->energy_ev));
      new_sigma = realloc (table->sigma, capacity * sizeof (*table->sigma));
      if (new_energy == NULL || new_sigma == NULL)
        return -1;
      table->energy_ev = new_energy;
      table->sigma = new_sigma;
    }
    table->energy_ev[table->n] = energy;
    table->sigma[table->n] = sigma;
    table->n++;
  }
  fclose (stream);
  return table->n >= 2 ? 0 : -1;
}

static double
sigma_interval (const CrossSection *table, size_t interval, double energy_ev)
{
  double e0, e1, s0, s1, fraction;

  if (interval == 0 && energy_ev < table->energy_ev[0])
  {
    fraction = (energy_ev - CHI_EV) / (table->energy_ev[0] - CHI_EV);
    return table->sigma[0] * pow (fmax (fraction, 0.0), 1.5);
  }
  e0 = table->energy_ev[interval - 1];
  e1 = table->energy_ev[interval];
  s0 = table->sigma[interval - 1];
  s1 = table->sigma[interval];
  fraction = (energy_ev - e0) / (e1 - e0);
  if (s0 > 0.0 && s1 > 0.0)
    return exp (log (s0) + fraction * log (s1 / s0));
  return s0 + fraction * (s1 - s0);
}

static double
milne_coefficient (const CrossSection *table, double temperature)
{
  const double g_hminus = 1.0, g_h_ground = 2.0;
  const double kt_erg = K_BOLTZMANN * temperature;
  const double saha0 = g_hminus / (2.0 * g_h_ground)
    * pow (H_PLANCK * H_PLANCK / (2.0 * M_PI * M_ELECTRON * kt_erg), 1.5);
  const double prefactor = 8.0 * M_PI * saha0
    / (C_LIGHT * C_LIGHT * H_PLANCK * H_PLANCK * H_PLANCK);
  double integral = 0.0;
  size_t interval;

  for (interval = 0; interval < table->n; interval++)
  {
    double lower = interval == 0 ? CHI_EV : table->energy_ev[interval - 1];
    double upper = fmin (table->energy_ev[interval], EMAX_EV);
    double midpoint, halfwidth, subtotal = 0.0;
    int node;
    if (!(upper > lower))
      continue;
    midpoint = 0.5 * (lower + upper);
    halfwidth = 0.5 * (upper - lower);
    for (node = 0; node < 8; node++)
    {
      double energy_ev = midpoint + halfwidth * gl_x[node];
      double energy_erg = energy_ev * EV_ERG;
      double boltzmann = exp (-(energy_ev - CHI_EV) * EV_ERG / kt_erg);
      subtotal += gl_w[node] * energy_erg * energy_erg
        * sigma_interval (table, interval, energy_ev) * boltzmann;
    }
    integral += halfwidth * EV_ERG * subtotal;
    if (table->energy_ev[interval] >= EMAX_EV)
      break;
  }
  return prefactor * integral;
}

int
main (int argc, char **argv)
{
  CrossSection table = { 0, NULL, NULL };
  FILE *temperatures, *output;
  char line[512];
  unsigned int index;
  double temperature;

  if (argc != 4)
  {
    fprintf (stderr, "usage: %s PD.csv temperature_nodes.csv output.csv\n", argv[0]);
    return EXIT_FAILURE;
  }
  if (read_cross_section (argv[1], &table) != 0)
    return EXIT_FAILURE;
  temperatures = fopen (argv[2], "r");
  output = fopen (argv[3], "w");
  if (temperatures == NULL || output == NULL)
  {
    fprintf (stderr, "unable to open temperature input or output\n");
    return EXIT_FAILURE;
  }
  fprintf (output, "T_K,k_RA_Milne_cm3_s\n");
  while (fgets (line, sizeof (line), temperatures) != NULL)
  {
    if (sscanf (line, "%u,%lf", &index, &temperature) != 2)
      continue;
    if (!(temperature > 0.0))
      return EXIT_FAILURE;
    fprintf (output, "%.17g,%.17e\n", temperature,
             milne_coefficient (&table, temperature));
  }
  fclose (temperatures);
  fclose (output);
  free (table.energy_ev);
  free (table.sigma);
  return EXIT_SUCCESS;
}
