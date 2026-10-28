/** @file passive_molecules.c
 *
 * Passive H2/H2+ diagnostics.  These routines deliberately do not alter the
 * active hydrogen, electron, or thermal ledgers.  The first implementation
 * stage owns and validates the closed 70-level H2/H collision data set.
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"

#define H2_H_LIQUE_NT 50
#define H2_H_BOSSION_NT 13
#define H2_H_PAIRS ((H2_PASSIVE_LEVELS * (H2_PASSIVE_LEVELS - 1)) / 2)
#define H2_ELECTRONIC_MAX 7
#define H2_V_MAX 38
#define H2_J_MAX 26
#define H2_UV_UPPER_STATES (H2_ELECTRONIC_MAX * H2_V_MAX * H2_J_MAX)

typedef enum h2_h_rate_source
{
  H2_H_LIQUE2015 = 0,
  H2_H_BOSSION2018 = 1
} H2HRateSource;

typedef struct h2_passive_level_data
{
  int v, j;
  double energy_k, g;
} H2PassiveLevelData;

typedef struct h2_h_rate_data
{
  int upper, lower;
  H2HRateSource source;
  int ntemperature;
  double rate[H2_H_LIQUE_NT];
  double spontaneous_rate;
} H2HRateData;

static H2PassiveLevelData h2_level_data[H2_PASSIVE_LEVELS];
static H2HRateData h2_h_rate_data[H2_H_PAIRS];
static int h2_passive_data_ready = FALSE;

typedef struct h2_uv_line_data
{
  int lower, upper_key;
  double frequency, oscillator_strength, dissociation_probability;
} H2UVLineData;

static H2UVLineData *h2_uv_line = NULL;
static size_t h2_uv_line_count = 0;
static double h2_uv_upper_total[H2_UV_UPPER_STATES];
static double h2_uv_upper_dissociation_probability[H2_UV_UPPER_STATES];
static double *h2_uv_upper_retained = NULL;
static int h2_uv_uses_lte_backsubstitution = FALSE;
static double h2_ad_distribution[7][H2_PASSIVE_LEVELS];
static const double h2_ad_temperature[7] = { 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0, 10000.0 };

static const double lique_temperature[H2_H_LIQUE_NT] = {
  100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
  1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
  2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000,
  3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000,
  4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000
};

static const double bossion_temperature[H2_H_BOSSION_NT] = {
  100, 150, 200, 300, 500, 700, 1000, 1500, 2000, 3000, 5000, 7000, 10000
};

static int h2_load_levels (const char *path);
static int h2_load_h_rates (const char *path);
static int h2_load_spontaneous_rates (const char *path);
static int h2_load_uv_data (const char *line_path, const char *dissociation_path);
static int h2_load_ad_distribution (const char *path);
static double h2_h_downward_rate (const H2HRateData *transition, double temperature);
static int h2_normalized_populations (double temperature, double neutral_h_density, double *population, double *residual);
static int h2_source_response (PlasmaPtr xplasma, double neutral_h_density, const double *source_fraction,
                               const double *destruction_rate, double *population, double *residual);
static int h2_h_collision_selftest (void);
static double passive_h2_photoionization_sigma (double energy_ev);
static double passive_h2plus_photodissociation_sigma (double energy_ev);

int
passive_molecules_init (void)
{
  const char *root;
  char path[LINELENGTH];

  if (!modes.passive_molecules_enabled)
    return 0;
  root = getenv ("SIROCCO");
  if (root == NULL || root[0] == '\0')
  {
    Error ("Passive molecules: SIROCCO is not set; cannot locate H2 data\n");
    return -1;
  }
  snprintf (path, sizeof (path), "%s/data/h2/passive70/hminus_deposit.dat", root);
  if (h2_load_ad_distribution (path) != 0)
    return -1;
  snprintf (path, sizeof (path), "%s/data/h2/passive70/levels70.csv", root);
  if (h2_load_levels (path) != 0)
    return -1;
  snprintf (path, sizeof (path), "%s/data/h2/passive70/downward_rates70.csv", root);
  if (h2_load_h_rates (path) != 0)
    return -1;
  snprintf (path, sizeof (path), "%s/data/h2/passive70/roueff_H2_IR_transitions.csv", root);
  if (h2_load_spontaneous_rates (path) != 0)
    return -1;
  {
    char dissociation_path[LINELENGTH];
    snprintf (path, sizeof (path), "%s/data/h2/passive70/cloudy_C25_H2_UV_lines.csv", root);
    snprintf (dissociation_path, sizeof (dissociation_path),
              "%s/data/h2/passive70/cloudy_C25_H2_UV_dissociation.csv", root);
    if (h2_load_uv_data (path, dissociation_path) != 0)
      return -1;
  }
  if (h2_h_collision_selftest () != 0)
    return -1;
  if (fabs (passive_h2_shielding_factor (0.0, 1.e5) - 1.0) > 1.e-14
      || passive_h2_shielding_factor (1.e22, 1.e5) >= 1.e-4)
  {
    Error ("Passive H2: Draine-Bertoldi shielding limiting-case self-test failed\n");
    return -1;
  }
  Log ("Passive molecular continuum check: sigma_H2_PI(20eV)=%.8e cm2 sigma_H2plus_PD(10.54eV)=%.8e cm2\n",
       passive_h2_photoionization_sigma (20.0), passive_h2plus_photodissociation_sigma (10.54));
  h2_passive_data_ready = TRUE;
  Log ("Passive H2 data: loaded 70 X-state levels and 2415 H-collider pairs; detailed-balance self-test passed\n");
  return 0;
}

static int
h2_load_levels (const char *path)
{
  FILE *stream;
  char *line = NULL, spin[16];
  size_t size = 0;
  int count = 0, index, v, j, meudon;
  double energy, g;

  stream = fopen (path, "r");
  if (stream == NULL)
  {
    Error ("Passive H2: unable to open %s\n", path);
    return -1;
  }
  if (getline (&line, &size, stream) < 0)
    goto malformed;
  while (getline (&line, &size, stream) >= 0)
  {
    if (sscanf (line, "%d,%d,%d,%15[^,],%lf,%lf,%d", &index, &v, &j, spin, &energy, &g, &meudon) != 7
        || index != count + 1 || meudon != index || energy < 0.0 || g <= 0.0)
      goto malformed;
    h2_level_data[count].v = v;
    h2_level_data[count].j = j;
    h2_level_data[count].energy_k = energy;
    h2_level_data[count].g = g;
    count++;
    if (count > H2_PASSIVE_LEVELS)
      goto malformed;
  }
  free (line);
  fclose (stream);
  if (count != H2_PASSIVE_LEVELS)
  {
    Error ("Passive H2: expected %d levels but read %d from %s\n", H2_PASSIVE_LEVELS, count, path);
    return -1;
  }
  return 0;

malformed:
  Error ("Passive H2: malformed level row %d in %s\n", count + 2, path);
  free (line);
  fclose (stream);
  return -1;
}

static int
h2_parse_rate_list (char *text, double *rate, int expected)
{
  int count = 0;
  char *cursor = text, *end;
  while (*cursor != '\0' && *cursor != '\r' && *cursor != '\n' && count < expected)
  {
    rate[count] = strtod (cursor, &end);
    if (end == cursor || rate[count] < 0.0)
      return -1;
    count++;
    cursor = end;
    if (*cursor == ';')
      cursor++;
    else if (*cursor != '\0' && *cursor != '\r' && *cursor != '\n')
      return -1;
  }
  return count == expected && (*cursor == '\0' || *cursor == '\r' || *cursor == '\n') ? 0 : -1;
}

static int
h2_load_h_rates (const char *path)
{
  FILE *stream;
  char *line = NULL, *field[10], *cursor;
  size_t size = 0;
  int count = 0, i, upper, lower;
  unsigned char seen[H2_PASSIVE_LEVELS][H2_PASSIVE_LEVELS] = { { 0 } };

  stream = fopen (path, "r");
  if (stream == NULL)
  {
    Error ("Passive H2: unable to open %s\n", path);
    return -1;
  }
  if (getline (&line, &size, stream) < 0)
    goto malformed;
  while (getline (&line, &size, stream) >= 0)
  {
    cursor = line;
    for (i = 0; i < 9; i++)
    {
      field[i] = cursor;
      cursor = strchr (cursor, ',');
      if (cursor == NULL)
        goto malformed;
      *cursor++ = '\0';
    }
    field[9] = cursor;
    upper = atoi (field[0]);
    lower = atoi (field[1]);
    if (upper <= lower || upper > H2_PASSIVE_LEVELS || lower < 1 || seen[upper - 1][lower - 1])
      goto malformed;
    seen[upper - 1][lower - 1] = 1;
    h2_h_rate_data[count].upper = upper - 1;
    h2_h_rate_data[count].lower = lower - 1;
    if (strcmp (field[7], "Lique2015") == 0)
    {
      h2_h_rate_data[count].source = H2_H_LIQUE2015;
      h2_h_rate_data[count].ntemperature = H2_H_LIQUE_NT;
    }
    else if (strcmp (field[7], "Bossion2018_via_public_Meudon") == 0)
    {
      h2_h_rate_data[count].source = H2_H_BOSSION2018;
      h2_h_rate_data[count].ntemperature = H2_H_BOSSION_NT;
    }
    else
      goto malformed;
    if (h2_parse_rate_list (field[9], h2_h_rate_data[count].rate, h2_h_rate_data[count].ntemperature) != 0)
      goto malformed;
    count++;
    if (count > H2_H_PAIRS)
      goto malformed;
  }
  free (line);
  fclose (stream);
  if (count != H2_H_PAIRS)
  {
    Error ("Passive H2: expected %d H-collider pairs but read %d from %s\n", H2_H_PAIRS, count, path);
    return -1;
  }
  return 0;

malformed:
  Error ("Passive H2: malformed H-collider row %d in %s\n", count + 2, path);
  free (line);
  fclose (stream);
  return -1;
}

static int
h2_find_level (int v, int j)
{
  int level;
  for (level = 0; level < H2_PASSIVE_LEVELS; level++)
    if (h2_level_data[level].v == v && h2_level_data[level].j == j)
      return level;
  return -1;
}

static H2HRateData *
h2_find_pair (int upper, int lower)
{
  int pair;
  if (upper < lower)
  {
    int swap = upper;
    upper = lower;
    lower = swap;
  }
  for (pair = 0; pair < H2_H_PAIRS; pair++)
    if (h2_h_rate_data[pair].upper == upper && h2_h_rate_data[pair].lower == lower)
      return &h2_h_rate_data[pair];
  return NULL;
}

static int
h2_load_spontaneous_rates (const char *path)
{
  FILE *stream;
  char *line = NULL, *field[16], *cursor;
  size_t size = 0;
  int row = 0, retained = 0;

  stream = fopen (path, "r");
  if (stream == NULL)
  {
    Error ("Passive H2: unable to open %s\n", path);
    return -1;
  }
  if (getline (&line, &size, stream) < 0)
    goto malformed;
  while (getline (&line, &size, stream) >= 0)
  {
    int i, upper, lower, v_upper, j_upper, v_lower, j_lower;
    double spontaneous;
    H2HRateData *pair;
    cursor = line;
    for (i = 0; i < 15; i++)
    {
      field[i] = cursor;
      cursor = strchr (cursor, ',');
      if (cursor == NULL)
        goto malformed;
      *cursor++ = '\0';
    }
    field[15] = cursor;
    v_upper = atoi (field[0]);
    j_upper = atoi (field[1]);
    v_lower = atoi (field[2]);
    j_lower = atoi (field[3]);
    spontaneous = strtod (field[8], NULL);
    if (spontaneous < 0.0 || !isfinite (spontaneous))
      goto malformed;
    upper = h2_find_level (v_upper, j_upper);
    lower = h2_find_level (v_lower, j_lower);
    if (upper >= 0 && lower >= 0)
    {
      if (h2_level_data[upper].energy_k <= h2_level_data[lower].energy_k)
        goto malformed;
      pair = h2_find_pair (upper, lower);
      if (pair == NULL)
        goto malformed;
      pair->spontaneous_rate += spontaneous;
      retained++;
    }
    row++;
  }
  free (line);
  fclose (stream);
  if (retained == 0)
  {
    Error ("Passive H2: no Roueff spontaneous transitions connected the retained 70 levels\n");
    return -1;
  }
  Log ("Passive H2 data: retained %d of %d Roueff spontaneous transitions inside the 70-level block\n", retained, row);
  return 0;

malformed:
  Error ("Passive H2: malformed Roueff transition row %d in %s\n", row + 2, path);
  free (line);
  fclose (stream);
  return -1;
}

static int
h2_load_uv_data (const char *line_path, const char *dissociation_path)
{
  FILE *stream;
  char *line = NULL;
  size_t size = 0, capacity = 0;
  double probability[H2_ELECTRONIC_MAX][H2_V_MAX][H2_J_MAX] = { { { 0.0 } } };
  unsigned char probability_present[H2_ELECTRONIC_MAX][H2_V_MAX][H2_J_MAX] = { { { 0 } } };
  int row = 0;

  stream = fopen (dissociation_path, "r");
  if (stream == NULL || getline (&line, &size, stream) < 0)
    goto malformed_dissociation;
  while (getline (&line, &size, stream) >= 0)
  {
    int electronic, v, j;
    double a_diss, a_bound, p_diss, kinetic;
    if (sscanf (line, "%d,%d,%d,%lf,%lf,%lf,%lf", &electronic, &v, &j,
                &a_diss, &a_bound, &p_diss, &kinetic) != 7
        || electronic < 0 || electronic >= H2_ELECTRONIC_MAX || v < 0 || v >= H2_V_MAX
        || j < 0 || j >= H2_J_MAX || p_diss < 0.0 || p_diss > 1.0 || probability_present[electronic][v][j])
      goto malformed_dissociation;
    probability[electronic][v][j] = p_diss;
    h2_uv_upper_total[(electronic * H2_V_MAX + v) * H2_J_MAX + j] = a_diss + a_bound;
    h2_uv_upper_dissociation_probability[(electronic * H2_V_MAX + v) * H2_J_MAX + j] = p_diss;
    probability_present[electronic][v][j] = 1;
    row++;
  }
  free (line);
  line = NULL;
  size = 0;
  fclose (stream);

  row = 0;
  h2_uv_upper_retained = calloc (H2_UV_UPPER_STATES * H2_PASSIVE_LEVELS,
                                 sizeof (*h2_uv_upper_retained));
  if (h2_uv_upper_retained == NULL)
    return -1;
  stream = fopen (line_path, "r");
  if (stream == NULL || getline (&line, &size, stream) < 0)
    goto malformed_line;
  while (getline (&line, &size, stream) >= 0)
  {
    int eu, vu, ju, el, vl, jl, lower;
    double a, wavenumber, wavelength, oscillator;
    if (sscanf (line, "%d,%d,%d,%d,%d,%d,%lf,%lf,%lf,%lf", &eu, &vu, &ju, &el, &vl, &jl,
                &a, &wavenumber, &wavelength, &oscillator) != 10
        || eu < 0 || eu >= H2_ELECTRONIC_MAX || vu < 0 || vu >= H2_V_MAX || ju < 0 || ju >= H2_J_MAX
        || !probability_present[eu][vu][ju] || oscillator < 0.0 || wavenumber <= 0.0)
      goto malformed_line;
    lower = el == 0 ? h2_find_level (vl, jl) : -1;
    if (lower >= 0)
      h2_uv_upper_retained[((eu * H2_V_MAX + vu) * H2_J_MAX + ju) * H2_PASSIVE_LEVELS + lower] += a;
    if (lower >= 0)
    {
      if (h2_uv_line_count == capacity)
      {
        H2UVLineData *new_data;
        capacity = capacity == 0 ? 4096 : 2 * capacity;
        new_data = realloc (h2_uv_line, capacity * sizeof (*h2_uv_line));
        if (new_data == NULL)
          goto malformed_line;
        h2_uv_line = new_data;
      }
      h2_uv_line[h2_uv_line_count].lower = lower;
      h2_uv_line[h2_uv_line_count].upper_key = (eu * H2_V_MAX + vu) * H2_J_MAX + ju;
      h2_uv_line[h2_uv_line_count].frequency = wavenumber * VLIGHT;
      h2_uv_line[h2_uv_line_count].oscillator_strength = oscillator;
      h2_uv_line[h2_uv_line_count].dissociation_probability = probability[eu][vu][ju];
      h2_uv_line_count++;
    }
    row++;
  }
  free (line);
  fclose (stream);
  if (h2_uv_line_count == 0)
    return -1;
  {
    int upper, represented = 0;
    double retained_sum = 0.0, missing_sum = 0.0;
    for (upper = 0; upper < H2_UV_UPPER_STATES; upper++)
      if (h2_uv_upper_total[upper] > 0.0)
      {
        int target;
        double retained = 0.0;
        for (target = 0; target < H2_PASSIVE_LEVELS; target++)
          retained += h2_uv_upper_retained[upper * H2_PASSIVE_LEVELS + target] / h2_uv_upper_total[upper];
        if (retained > 1.0 + 1.e-6)
        {
          Error ("Passive H2: retained UV branch %.8e exceeds unity for upper key %d\n", retained, upper);
          return -1;
        }
        retained_sum += retained;
        missing_sum += fmax (1.0 - h2_uv_upper_dissociation_probability[upper] - retained, 0.0);
        represented++;
      }
    Log ("Passive H2 UV closure: %d upper states mean direct-retained branch %.5f; remainder %.5f uses LTE back-substitution proxy\n",
         represented, retained_sum / fmax (represented, 1), missing_sum / fmax (represented, 1));
    h2_uv_uses_lte_backsubstitution = missing_sum > 0.0;
  }
  Log ("Passive H2 data: retained %zu of %d Cloudy-C25 UV lines from the 70-level block\n", h2_uv_line_count, row);
  return 0;

malformed_dissociation:
  Error ("Passive H2: malformed UV dissociation row %d in %s\n", row + 2, dissociation_path);
  free (line);
  if (stream != NULL)
    fclose (stream);
  return -1;
malformed_line:
  Error ("Passive H2: malformed UV line row %d in %s\n", row + 2, line_path);
  free (line);
  if (stream != NULL)
    fclose (stream);
  return -1;
}

static int
h2_load_ad_distribution (const char *path)
{
  FILE *stream;
  char line[LINELENGTH];
  double sum[7] = { 0.0 };
  int magic1, magic2, magic3, rows = 0;

  stream = fopen (path, "r");
  if (stream == NULL || fgets (line, sizeof (line), stream) == NULL
      || sscanf (line, "%d %d %d", &magic1, &magic2, &magic3) != 3
      || magic1 != 2 || magic2 != 10 || magic3 != 17)
  {
    Error ("Passive H2: invalid associative-detachment distribution header in %s\n", path);
    if (stream != NULL)
      fclose (stream);
    return -1;
  }
  while (fgets (line, sizeof (line), stream) != NULL)
  {
    int v, j, level, itemp;
    double energy, a[7];
    if (line[0] == '#')
      continue;
    if (sscanf (line, "%d %d", &v, &j) == 2 && v < 0)
      break;
    if (sscanf (line, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf", &v, &j, &energy,
                &a[0], &a[1], &a[2], &a[3], &a[4], &a[5], &a[6]) != 10)
    {
      Error ("Passive H2: malformed associative-detachment distribution row %d in %s\n", rows + 1, path);
      fclose (stream);
      return -1;
    }
    level = h2_find_level (v, j);
    if (level >= 0)
    {
      for (itemp = 0; itemp < 7; itemp++)
      {
        h2_ad_distribution[itemp][level] = pow (10.0, -a[itemp]);
        sum[itemp] += h2_ad_distribution[itemp][level];
      }
    }
    rows++;
  }
  fclose (stream);
  for (int itemp = 0; itemp < 7; itemp++)
  {
    int level;
    if (!(sum[itemp] > 0.0))
      return -1;
    for (level = 0; level < H2_PASSIVE_LEVELS; level++)
      h2_ad_distribution[itemp][level] /= sum[itemp];
  }
  Log ("Passive H2 data: loaded Cloudy/Launay91 associative-detachment weights and renormalized onto 70 levels\n");
  return 0;
}

void
passive_h2_ad_formation_distribution (double temperature, double *distribution)
{
  int lo, level;
  double upper_weight, sum = 0.0;
  if (temperature <= h2_ad_temperature[0])
  {
    lo = 0;
    upper_weight = 0.0;
  }
  else if (temperature >= h2_ad_temperature[6])
  {
    lo = 5;
    upper_weight = 1.0;
  }
  else
  {
    for (lo = 0; lo < 6 && temperature > h2_ad_temperature[lo + 1]; lo++);
    upper_weight = log (temperature / h2_ad_temperature[lo])
      / log (h2_ad_temperature[lo + 1] / h2_ad_temperature[lo]);
  }
  for (level = 0; level < H2_PASSIVE_LEVELS; level++)
  {
    distribution[level] = (1.0 - upper_weight) * h2_ad_distribution[lo][level]
      + upper_weight * h2_ad_distribution[lo + 1][level];
    sum += distribution[level];
  }
  if (sum > 0.0)
    for (level = 0; level < H2_PASSIVE_LEVELS; level++)
      distribution[level] /= sum;
}

double
passive_h2_lw_dissociation_rate (PlasmaPtr xplasma, int level)
{
  double rate = 0.0;
  size_t line_index;
  if (!h2_passive_data_ready || level < 0 || level >= H2_PASSIVE_LEVELS)
    return 0.0;
  for (line_index = 0; line_index < h2_uv_line_count; line_index++)
  {
    H2UVLineData *transition = &h2_uv_line[line_index];
    int bin;
    double jnu;
    if (transition->lower != level)
      continue;
    bin = (int) floor ((log10 (transition->frequency) - geo.cell_log_freq_min) / geo.cell_delta_lfreq);
    if (bin < 0 || bin >= NBINS_IN_CELL_SPEC)
      continue;
    jnu = xplasma->cell_spec_flux[bin];
    if (jnu <= 0.0)
      continue;
    rate += 4.0 * PI * jnu / (PLANCK * transition->frequency) * PI_E2_OVER_MC
      * transition->oscillator_strength * transition->dissociation_probability;
  }
  return rate;
}

double
passive_h2_shielding_factor (double column, double doppler_b)
{
  double x, b5, root;
  if (column <= 0.0)
    return 1.0;
  b5 = fmax (doppler_b, 1.0) / 1.e5;
  x = column / 5.e14;
  root = sqrt (1.0 + x);
  return 0.965 / pow (1.0 + x / b5, 2.0) + 0.035 / root * exp (-8.5e-4 * root);
}

static double
passive_h2_photoionization_sigma (double energy_ev)
{
  double energy_kev, x, root_x, x15, x2, cross_section;
  if (energy_ev < 15.4)
    return 0.0;
  x = energy_ev / 15.4;
  root_x = sqrt (x);
  x15 = x * root_x;
  x2 = x * x;
  energy_kev = energy_ev / 1000.0;
  if (energy_ev < 18.0)
    cross_section = 1.e7 * (1.0 - 197.448 / root_x + 438.823 / x - 260.481 / x15 + 17.915 / x2);
  else if (energy_ev <= 30.0)
    cross_section = (-145.528 + 351.394 * root_x - 274.294 * x + 74.320 * x15) / pow (energy_kev, 3.5);
  else if (energy_ev <= 85.0)
    cross_section = (65.304 - 91.762 * root_x + 51.778 * x - 9.364 * x15) / pow (energy_kev, 3.5);
  else
    cross_section = 45.57 * (1.0 - 2.003 / root_x - 4.806 / x + 50.577 / x15 - 171.044 / x2
                             + 231.608 / (root_x * x2) - 81.885 / (x * x2)) / pow (energy_kev, 3.5);
  return fmax (cross_section, 0.0) * 1.e-24;
}

static double
passive_h2plus_photodissociation_sigma (double energy_ev)
{
  static const double energy[] = { 6.75, 8.68, 10.54, 12.46, 14.28 };
  static const double cross_section_mb[] = { 0.24, 2.5, 7.1, 6.0, 2.7 };
  int index;
  if (energy_ev < energy[0] || energy_ev > energy[4])
    return 0.0;
  for (index = 0; index < 4 && energy_ev > energy[index + 1]; index++);
  return (cross_section_mb[index] + (cross_section_mb[index + 1] - cross_section_mb[index])
          * (energy_ev - energy[index]) / (energy[index + 1] - energy[index])) * 1.e-18;
}

static double
passive_continuum_photo_rate (PlasmaPtr xplasma, double (*cross_section) (double))
{
  double rate = 0.0;
  int bin;
  for (bin = 0; bin < NBINS_IN_CELL_SPEC; bin++)
  {
    double log_lo = geo.cell_log_freq_min + bin * geo.cell_delta_lfreq;
    double log_hi = log_lo + geo.cell_delta_lfreq;
    double freq_lo = pow (10.0, log_lo), freq_hi = pow (10.0, log_hi);
    double freq = sqrt (freq_lo * freq_hi);
    double sigma = cross_section (PLANCK * freq / EV2ERGS);
    if (sigma > 0.0 && xplasma->cell_spec_flux[bin] > 0.0)
      rate += 4.0 * PI * xplasma->cell_spec_flux[bin] * sigma / (PLANCK * freq) * (freq_hi - freq_lo);
  }
  return rate;
}

static double
passive_continuum_absorbed_power (PlasmaPtr xplasma, double (*cross_section) (double))
{
  double power = 0.0;
  int bin;
  for (bin = 0; bin < NBINS_IN_CELL_SPEC; bin++)
  {
    double log_lo = geo.cell_log_freq_min + bin * geo.cell_delta_lfreq;
    double log_hi = log_lo + geo.cell_delta_lfreq;
    double freq_lo = pow (10.0, log_lo), freq_hi = pow (10.0, log_hi);
    double freq = sqrt (freq_lo * freq_hi);
    double sigma = cross_section (PLANCK * freq / EV2ERGS);
    if (sigma > 0.0 && xplasma->cell_spec_flux[bin] > 0.0)
      power += 4.0 * PI * xplasma->cell_spec_flux[bin] * sigma * (freq_hi - freq_lo);
  }
  return power;
}

double
passive_h2_photoionization_rate (PlasmaPtr xplasma)
{
  return passive_continuum_photo_rate (xplasma, passive_h2_photoionization_sigma);
}

double
passive_h2plus_photodissociation_rate (PlasmaPtr xplasma)
{
  return passive_continuum_photo_rate (xplasma, passive_h2plus_photodissociation_sigma);
}

/* Cloudy C25 coarse gas-phase rates retained by the passive diagnostic. */
static double
h2plus_radiative_association (double temperature)
{
  /* Cloudy C25 standard (non-Leiden) radath closure.  The often quoted
   * 5.3e-19*(T/300)^1.85 expression belongs to its Leiden comparison branch. */
  return fmax (0.0, 2.325 * fmin (5000.0, fmax (temperature, 1.0)) - 1375.0) * 1.e-20;
}

static double
h2plus_dissociative_recombination (double temperature)
{
  return 1.6e-8 * pow (fmax (temperature, 1.0), -0.43);
}

static double
hminus_to_h2plus_forward (double temperature)
{
  temperature = fmax (temperature, 1.0);
  return temperature <= 7891.0 ? 6.9e-9 * pow (temperature, -0.35)
    : 9.6e-7 * pow (temperature, -0.9);
}

static double
hminus_to_h2plus_reverse (double temperature)
{
  const double reaction_energy_ev = 1.8965;
  double mass_factor, equilibrium_constant, exponent;
  temperature = fmax (temperature, 1.0);
  /* Effective scalar partition convention: all unresolved internal weights
   * are one.  Translational masses and the agreed chemical energy zero then
   * define n(H2+)ne/[n(H-)n(H+)]. */
  mass_factor = pow (2.0 * MELEC / MPROT, 1.5);
  exponent = fmin (reaction_energy_ev * EV2ERGS / (BOLTZMANN * temperature), 700.0);
  equilibrium_constant = mass_factor * exp (exponent);
  return hminus_to_h2plus_forward (temperature) / fmax (equilibrium_constant, 1.e-300);
}

static double
h2_level_collisional_dissociation (int level, double temperature, double collider_density)
{
  const double dissociation_temperature = 51997.0;
  double binding_temperature = fmax (dissociation_temperature - h2_level_data[level].energy_k, 0.0);
  return collider_density * 1.e-14 * exp (-binding_temperature / fmax (temperature, 1.0));
}

static int
h2_vj_level (int v, int j)
{
  int level = h2_find_level (v, j);
  return level >= 0 ? level : 0;
}

static int
h2_source_response (PlasmaPtr xplasma, double neutral_h_density, const double *source_fraction,
                    const double *destruction_rate, double *population, double *residual)
{
  double *matrix, *matrix_original, rhs[H2_PASSIVE_LEVELS] = { 0.0 };
  double max_residual = 0.0;
  int pair, level, row, error;

  memset (population, 0, H2_PASSIVE_LEVELS * sizeof (*population));
  *residual = 0.0;
  matrix = calloc (H2_PASSIVE_LEVELS * H2_PASSIVE_LEVELS, sizeof (*matrix));
  matrix_original = malloc (H2_PASSIVE_LEVELS * H2_PASSIVE_LEVELS * sizeof (*matrix_original));
  if (matrix == NULL || matrix_original == NULL)
  {
    free (matrix);
    free (matrix_original);
    return -1;
  }
  for (level = 0; level < H2_PASSIVE_LEVELS; level++)
    rhs[level] = -source_fraction[level];

  for (pair = 0; pair < H2_H_PAIRS; pair++)
  {
    H2HRateData *transition = &h2_h_rate_data[pair];
    int upper = transition->upper, lower = transition->lower;
    double down_collision = neutral_h_density * h2_h_downward_rate (transition, xplasma->t_e);
    double delta = h2_level_data[upper].energy_k - h2_level_data[lower].energy_k;
    double up = down_collision * h2_level_data[upper].g / h2_level_data[lower].g
      * exp (-delta / fmax (xplasma->t_e, 1.0));
    double down = down_collision + transition->spontaneous_rate;
    matrix[lower * H2_PASSIVE_LEVELS + upper] += down;
    matrix[upper * H2_PASSIVE_LEVELS + upper] -= down;
    matrix[upper * H2_PASSIVE_LEVELS + lower] += up;
    matrix[lower * H2_PASSIVE_LEVELS + lower] -= up;
  }
  /* Eliminate the short-lived electronic upper state.  Direct decays into
   * retained levels use the Cloudy A-values.  The remaining bound branch,
   * which lands in omitted X-state levels, is returned with a retained LTE
   * closure rather than being mistaken for molecular destruction. */
  {
    size_t line_index;
    double lte[H2_PASSIVE_LEVELS], lte_sum = 0.0;
    for (level = 0; level < H2_PASSIVE_LEVELS; level++)
    {
      lte[level] = h2_level_data[level].g
        * exp (-h2_level_data[level].energy_k / fmax (xplasma->t_e, 1.0));
      lte_sum += lte[level];
    }
    for (level = 0; level < H2_PASSIVE_LEVELS; level++)
      lte[level] /= lte_sum;
    for (line_index = 0; line_index < h2_uv_line_count; line_index++)
    {
      H2UVLineData *transition = &h2_uv_line[line_index];
      double total_a = h2_uv_upper_total[transition->upper_key], retained_probability = 0.0;
      double base_rate, missing_probability, jnu;
      int bin, target;
      bin = (int) floor ((log10 (transition->frequency) - geo.cell_log_freq_min) / geo.cell_delta_lfreq);
      if (bin < 0 || bin >= NBINS_IN_CELL_SPEC || total_a <= 0.0)
        continue;
      jnu = xplasma->cell_spec_flux[bin];
      if (jnu <= 0.0)
        continue;
      base_rate = 4.0 * PI * jnu / (PLANCK * transition->frequency) * PI_E2_OVER_MC
        * transition->oscillator_strength;
      for (target = 0; target < H2_PASSIVE_LEVELS; target++)
        retained_probability += h2_uv_upper_retained[transition->upper_key * H2_PASSIVE_LEVELS + target] / total_a;
      missing_probability = fmax (1.0 - transition->dissociation_probability - retained_probability, 0.0);
      for (target = 0; target < H2_PASSIVE_LEVELS; target++)
      {
        double branch = h2_uv_upper_retained[transition->upper_key * H2_PASSIVE_LEVELS + target] / total_a
          + missing_probability * lte[target];
        double rate = base_rate * branch;
        matrix[target * H2_PASSIVE_LEVELS + transition->lower] += rate;
        matrix[transition->lower * H2_PASSIVE_LEVELS + transition->lower] -= rate;
      }
    }
  }

  for (level = 0; level < H2_PASSIVE_LEVELS; level++)
  {
    matrix[level * H2_PASSIVE_LEVELS + level] -= destruction_rate[level];
  }
  memcpy (matrix_original, matrix, H2_PASSIVE_LEVELS * H2_PASSIVE_LEVELS * sizeof (*matrix));
  /* SIROCCO's matrix wrapper checks an absolute equation residual.  Molecular
   * rates span many decades, so scale each equation without changing its
   * solution and retain matrix_original for the physical relative residual. */
  for (row = 0; row < H2_PASSIVE_LEVELS; row++)
  {
    double row_scale = fabs (rhs[row]);
    for (level = 0; level < H2_PASSIVE_LEVELS; level++)
      row_scale = fmax (row_scale, fabs (matrix[row * H2_PASSIVE_LEVELS + level]));
    if (row_scale > 0.0)
    {
      rhs[row] /= row_scale;
      for (level = 0; level < H2_PASSIVE_LEVELS; level++)
        matrix[row * H2_PASSIVE_LEVELS + level] /= row_scale;
    }
  }
  error = solve_matrix (matrix, rhs, H2_PASSIVE_LEVELS, population, -1);
  if (error != 0)
  {
    free (matrix);
    free (matrix_original);
    return -1;
  }
  for (level = 0; level < H2_PASSIVE_LEVELS; level++)
  {
    if (!isfinite (population[level]) || population[level] < -1.e-12)
    {
      free (matrix);
      free (matrix_original);
      return -1;
    }
    population[level] = fmax (population[level], 0.0);
  }
  for (row = 0; row < H2_PASSIVE_LEVELS; row++)
  {
    double balance = source_fraction[row], scale = fabs (balance);
    for (level = 0; level < H2_PASSIVE_LEVELS; level++)
    {
      double term = matrix_original[row * H2_PASSIVE_LEVELS + level] * population[level];
      balance += term;
      scale += fabs (term);
    }
    if (scale > 0.0)
      max_residual = fmax (max_residual, fabs (balance) / scale);
  }
  free (matrix);
  free (matrix_original);
  *residual = max_residual;
  return 0;
}

typedef struct passive_radial_cell
{
  int nplasma;
  double radius;
} PassiveRadialCell;

static int
passive_radius_compare (const void *left, const void *right)
{
  const PassiveRadialCell *a = left, *b = right;
  return a->radius < b->radius ? -1 : a->radius > b->radius;
}

static double
passive_cell_path_length (PlasmaPtr xplasma)
{
  WindPtr cell = &wmain[xplasma->nwind];
  if (cell->dr > 0.0)
    return cell->dr;
  return cbrt (fmax (xplasma->vol, 0.0));
}

static void
passive_hydrogen_densities (PlasmaPtr xplasma, double *n_hi, double *n_hp)
{
  int element, first = -1;
  for (element = 0; element < nelements; element++)
    if (ele[element].z == 1)
    {
      first = ele[element].firstion;
      break;
    }
  *n_hi = first >= 0 && first < nions ? xplasma->density[first] : 0.0;
  *n_hp = first + 1 >= 0 && first + 1 < nions ? xplasma->density[first + 1] : 0.0;
}

static void
passive_hydrogen_level_densities (PlasmaPtr xplasma, double *n_h1, double *n_h2)
{
  int element, first = -1, config;
  *n_h1 = *n_h2 = 0.0;
  for (element = 0; element < nelements; element++)
    if (ele[element].z == 1)
    {
      first = ele[element].firstion;
      break;
    }
  if (first < 0 || first >= nions || ion[first].first_nlte_level < 0)
    return;
  for (config = ion[first].first_nlte_level;
       config < ion[first].first_nlte_level + ion[first].nlte; config++)
  {
    /* levden stores a fraction within the ion; den_config returns the
     * physical configuration density required by bimolecular rates. */
    double density = den_config (xplasma, config);
    if (xconfig[config].ilv == 1 || fabs (xconfig[config].q_num - 1.0) < 1.e-6)
      *n_h1 += density;
    if (xconfig[config].ilv == 2 || fabs (xconfig[config].q_num - 2.0) < 1.e-6)
      *n_h2 += density;
  }
}

static double
passive_h2_line_cooling (const double *population)
{
  double cooling = 0.0;
  int pair;
  for (pair = 0; pair < H2_H_PAIRS; pair++)
  {
    H2HRateData *transition = &h2_h_rate_data[pair];
    double delta_energy = BOLTZMANN * (h2_level_data[transition->upper].energy_k
                                      - h2_level_data[transition->lower].energy_k);
    cooling += population[transition->upper] * transition->spontaneous_rate * delta_energy;
  }
  return cooling;
}

static int
passive_solve_cell (PlasmaPtr xplasma, double shielding, int unshielded)
{
  PassiveMoleculeCell *molecule = &xplasma->molecules;
  double *population = unshielded ? molecule->h2_level_unshielded : molecule->h2_level;
  double response_ad[H2_PASSIVE_LEVELS], response_ct[H2_PASSIVE_LEVELS], response_lte[H2_PASSIVE_LEVELS];
  double response_threebody[H2_PASSIVE_LEVELS];
  double source_ad[H2_PASSIVE_LEVELS], source_ct[H2_PASSIVE_LEVELS] = { 0.0 };
  double source_lte[H2_PASSIVE_LEVELS], source_threebody[H2_PASSIVE_LEVELS], lte_sum = 0.0, threebody_sum = 0.0;
  double n_hi, n_hp, n_h1, n_h2, h2plus = 0.0, residual_ad = 0.0, residual_ct = 0.0;
  double residual_lte = 0.0, residual_threebody = 0.0, destruction[H2_PASSIVE_LEVELS];
  double lw_loss[H2_PASSIVE_LEVELS], collisional_loss[H2_PASSIVE_LEVELS];
  double formation_ad = xplasma->hminus.passive_ad_rate;
  double formation_lte, formation_threebody, h2_photoion, h2_total = 0.0, h2plus_source, h2plus_loss;
  double response_ad_total = 0.0, response_ct_total = 0.0, response_lte_total = 0.0;
  double response_threebody_total = 0.0, ct_coefficient, denominator, collider_density;
  double hminus_reverse_coefficient;
  int level, v4j0 = h2_vj_level (4, 0);

  passive_hydrogen_densities (xplasma, &n_hi, &n_hp);
  passive_hydrogen_level_densities (xplasma, &n_h1, &n_h2);
  passive_h2_ad_formation_distribution (xplasma->t_e, source_ad);
  source_ct[v4j0] = 1.0;
  for (level = 0; level < H2_PASSIVE_LEVELS; level++)
  {
    source_lte[level] = h2_level_data[level].g
      * exp (-h2_level_data[level].energy_k / fmax (xplasma->t_e, 1.0));
    lte_sum += source_lte[level];
    /* Cloudy big-H2 uses its level-specific dissociation coefficient and
     * the corresponding LTE ratio for the reverse three-body source.  With
     * its current 1e-14 exp[-D_bind/T] closure, the exponent cancels and the
     * retained-level source is proportional to g. */
    source_threebody[level] = h2_level_data[level].g;
    threebody_sum += source_threebody[level];
  }
  for (level = 0; level < H2_PASSIVE_LEVELS; level++)
  {
    source_lte[level] /= lte_sum;
    source_threebody[level] /= threebody_sum;
  }
  /* Cloudy C25 standard radasc is H(1s)+H(n=2), despite the compact
   * chemistry label H+H.  Its unresolved product distribution is assigned
   * to the retained LTE closure, as specified for residual formation. */
  formation_lte = 3.e-14 * n_h1 * n_h2;
  formation_threebody = n_hi * n_hi * n_hi * 1.e-14 * SAHA
    / pow (fmax (xplasma->t_e, 1.0), 1.5) * 0.25 * 3.634e-5 * threebody_sum;
  h2_photoion = passive_h2_photoionization_rate (xplasma);
  h2plus_source = xplasma->hminus.passive_h2plus_rate
    + h2plus_radiative_association (xplasma->t_e) * n_hi * n_hp
    + 3.27e-10 * pow (fmax (xplasma->t_e, 1.0) / 300.0, -0.35)
      * exp (-17829.0 / fmax (xplasma->t_e, 1.0)) * n_h1 * n_h2;
  hminus_reverse_coefficient = hminus_to_h2plus_reverse (xplasma->t_e);
  h2plus_loss = passive_h2plus_photodissociation_rate (xplasma)
    + 6.4e-10 * n_hi + h2plus_dissociative_recombination (xplasma->t_e) * xplasma->ne
    + hminus_reverse_coefficient * xplasma->ne;
  h2plus_loss = fmax (h2plus_loss, 1.e-40);
  collider_density = n_hi;
  for (level = 0; level < H2_PASSIVE_LEVELS; level++)
  {
    lw_loss[level] = shielding * passive_h2_lw_dissociation_rate (xplasma, level);
    collisional_loss[level] = h2_level_collisional_dissociation (level, xplasma->t_e, collider_density);
    destruction[level] = lw_loss[level] + h2_photoion + collisional_loss[level];
    destruction[level] = fmax (destruction[level], 1.e-40);
  }
  if (h2_source_response (xplasma, n_hi, source_ad, destruction, response_ad, &residual_ad) != 0
      || h2_source_response (xplasma, n_hi, source_ct, destruction, response_ct, &residual_ct) != 0
      || h2_source_response (xplasma, n_hi, source_lte, destruction, response_lte, &residual_lte) != 0
      || h2_source_response (xplasma, n_hi, source_threebody, destruction,
                             response_threebody, &residual_threebody) != 0)
    return -1;
  for (level = 0; level < H2_PASSIVE_LEVELS; level++)
  {
    response_ad_total += response_ad[level];
    response_ct_total += response_ct[level];
    response_lte_total += response_lte[level];
    response_threebody_total += response_threebody[level];
  }
  ct_coefficient = 6.4e-10 * n_hi;
  denominator = h2plus_loss - h2_photoion * ct_coefficient * response_ct_total;
  if (!(denominator > 0.0) || !isfinite (denominator))
  {
    molecule->flags |= 16u;
    denominator = h2plus_loss;
  }
  h2plus = (h2plus_source + h2_photoion
            * (formation_ad * response_ad_total + formation_lte * response_lte_total
               + formation_threebody * response_threebody_total)) / denominator;
  for (level = 0; level < H2_PASSIVE_LEVELS; level++)
  {
    population[level] = formation_ad * response_ad[level]
      + ct_coefficient * h2plus * response_ct[level] + formation_lte * response_lte[level]
      + formation_threebody * response_threebody[level];
    h2_total += population[level];
  }
  if (unshielded)
  {
    molecule->h2_total_unshielded = h2_total;
    return 0;
  }
  molecule->h2_total = h2_total;
  molecule->h2plus = h2plus;
  molecule->matrix_residual = fmax (residual_threebody, fmax (residual_lte, fmax (residual_ad, residual_ct)));
  molecule->shielding_factor = shielding;
  molecule->rate_formation = formation_ad + formation_lte + formation_threebody + 6.4e-10 * n_hi * h2plus;
  molecule->form_ad = formation_ad;
  molecule->form_charge_transfer = 6.4e-10 * n_hi * h2plus;
  molecule->form_excited_ra = formation_lte;
  molecule->form_threebody = formation_threebody;
  molecule->rate_destruction = molecule->destroy_lw = molecule->destroy_collisional = 0.0;
  molecule->destroy_photoionization = h2_photoion * h2_total;
  for (level = 0; level < H2_PASSIVE_LEVELS; level++)
  {
    molecule->rate_destruction += population[level] * destruction[level];
    molecule->destroy_lw += population[level] * lw_loss[level];
    molecule->destroy_collisional += population[level] * collisional_loss[level];
  }
  molecule->h2plus_form_hminus = xplasma->hminus.passive_h2plus_rate;
  molecule->h2plus_form_ra = h2plus_radiative_association (xplasma->t_e) * n_hi * n_hp;
  molecule->h2plus_form_assoc_ion = h2plus_source - molecule->h2plus_form_hminus - molecule->h2plus_form_ra;
  molecule->h2plus_form_photoionization = h2_photoion * h2_total;
  molecule->h2plus_destroy_pd = h2plus * passive_h2plus_photodissociation_rate (xplasma);
  molecule->h2plus_destroy_ct = h2plus * 6.4e-10 * n_hi;
  molecule->h2plus_destroy_dr = h2plus * h2plus_dissociative_recombination (xplasma->t_e) * xplasma->ne;
  molecule->h2plus_destroy_hminus_reverse = h2plus * hminus_reverse_coefficient * xplasma->ne;
  molecule->cooling_line = passive_h2_line_cooling (population) * xplasma->vol;
  molecule->lum_abs = xplasma->vol
    * (h2_total * passive_continuum_absorbed_power (xplasma, passive_h2_photoionization_sigma)
       + h2plus * passive_continuum_absorbed_power (xplasma, passive_h2plus_photodissociation_sigma));
  molecule->lum_net = molecule->lum_abs - molecule->lum_emit;
  {
    double forward = xplasma->hminus.passive_h2plus_rate;
    double reverse = h2plus * xplasma->ne * hminus_reverse_coefficient;
    double signed_power = (forward - reverse) * 1.8965 * EV2ERGS * xplasma->vol;
    molecule->heat_chemical = fmax (signed_power, 0.0);
    molecule->cool_chemical = fmax (-signed_power, 0.0);
  }
  molecule->trace_fraction = (2.0 * h2_total + 2.0 * h2plus) / fmax (n_hi + n_hp, 1.e-100);
  if (xplasma->cool_tot > 0.0)
    molecule->passive_cooling_fraction = molecule->cooling_line / xplasma->cool_tot;
  else
  {
    molecule->passive_cooling_fraction = 0.0;
    molecule->flags |= 32u;
  }
  if (molecule->trace_fraction > 1.e-2)
    molecule->flags |= 1u;
  if (molecule->matrix_residual > 1.e-6)
    molecule->flags |= 2u;
  if (h2_uv_uses_lte_backsubstitution)
    molecule->flags |= 128u;
  /* No frequency-resolved H2/H2+ continuum emission kernel is available.
   * Mark L_emit=0 as unresolved rather than physically zero. */
  molecule->flags |= 64u;
  return 0;
}

int
passive_molecules_update_all (void)
{
  PassiveRadialCell *order;
  double max_change = 0.0, total_h2_column = 0.0, total_h2plus_column = 0.0;
  int count = 0, index, iteration;

  if (!modes.passive_molecules_enabled || !h2_passive_data_ready)
    return 0;
  order = calloc (NPLASMA, sizeof (*order));
  if (order == NULL)
    return -1;
  for (index = 0; index < NPLASMA; index++)
    if (plasmamain[index].nwind >= 0 && plasmamain[index].nwind < NDIM2
        && wmain[plasmamain[index].nwind].inwind >= 0)
    {
      order[count].nplasma = index;
      order[count].radius = wmain[plasmamain[index].nwind].rcen;
      count++;
    }
  qsort (order, count, sizeof (*order), passive_radius_compare);

  for (index = 0; index < count; index++)
  {
    memset (&plasmamain[order[index].nplasma].molecules, 0, sizeof (PassiveMoleculeCell));
    passive_solve_cell (&plasmamain[order[index].nplasma], 1.0, TRUE);
  }
  for (iteration = 0; iteration < 30; iteration++)
  {
    double column = 0.0;
    max_change = 0.0;
    for (index = 0; index < count; index++)
    {
      PlasmaPtr xplasma = &plasmamain[order[index].nplasma];
      double path = passive_cell_path_length (xplasma);
      double old = xplasma->molecules.h2_total;
      double center_column = column + 0.5 * old * path;
      double thermal_b = sqrt (2.0 * BOLTZMANN * fmax (xplasma->t_e, 1.0) / (2.0 * MPROT));
      double shielding = passive_h2_shielding_factor (center_column, thermal_b);
      xplasma->molecules.h2_column_in = center_column;
      if (passive_solve_cell (xplasma, shielding, FALSE) < 0)
        xplasma->molecules.flags |= 4u;
      max_change = fmax (max_change, fabs (xplasma->molecules.h2_total - old)
                         / fmax (xplasma->molecules.h2_total, 1.e-100));
      column += 0.5 * (old + xplasma->molecules.h2_total) * path;
    }
    if (max_change < 1.e-5)
      break;
  }
  if (iteration >= 30)
    for (index = 0; index < count; index++)
      plasmamain[order[index].nplasma].molecules.flags |= 8u;
  for (index = 0; index < count; index++)
  {
    PlasmaPtr xplasma = &plasmamain[order[index].nplasma];
    double path = passive_cell_path_length (xplasma);
    total_h2_column += xplasma->molecules.h2_total * path;
    total_h2plus_column += xplasma->molecules.h2plus * path;
  }
  {
    double outward_h2 = total_h2_column, outward_h2plus = total_h2plus_column;
    for (index = 0; index < count; index++)
    {
      PlasmaPtr xplasma = &plasmamain[order[index].nplasma];
      double path = passive_cell_path_length (xplasma);
      double center_h2 = outward_h2 - 0.5 * xplasma->molecules.h2_total * path;
      double center_h2plus = outward_h2plus - 0.5 * xplasma->molecules.h2plus * path;
      /* Continuum maxima over the adopted fits: the H2+ table peaks at
       * 10.54 eV; H2 is sampled densely because its analytic fit is not
       * monotonic close to threshold. */
      double max_h2_sigma = 0.0;
      int sample;
      for (sample = 0; sample <= 400; sample++)
      {
        double energy = 15.4 * pow (1000.0 / 15.4, sample / 400.0);
        max_h2_sigma = fmax (max_h2_sigma, passive_h2_photoionization_sigma (energy));
      }
      xplasma->molecules.tau_h2_max = fmax (center_h2, 0.0) * max_h2_sigma;
      xplasma->molecules.tau_h2plus_max = fmax (center_h2plus, 0.0) * 7.1e-18;
      outward_h2 -= xplasma->molecules.h2_total * path;
      outward_h2plus -= xplasma->molecules.h2plus * path;
    }
  }
  Log ("Passive H2 update: cells=%d shielding_iterations=%d max_relative_change=%.3e\n",
       count, iteration + 1, max_change);
  free (order);
  return iteration < 30 ? 0 : 1;
}

static double
h2_h_downward_rate (const H2HRateData *transition, double temperature)
{
  const double *grid;
  int n, lo;
  double fraction, y0, y1;

  grid = transition->source == H2_H_LIQUE2015 ? lique_temperature : bossion_temperature;
  n = transition->ntemperature;
  if (temperature <= grid[0])
    return transition->rate[0];
  if (temperature >= grid[n - 1])
    return transition->rate[n - 1];
  for (lo = 0; lo + 1 < n && temperature > grid[lo + 1]; lo++);
  fraction = log (temperature / grid[lo]) / log (grid[lo + 1] / grid[lo]);
  y0 = transition->rate[lo];
  y1 = transition->rate[lo + 1];
  if (y0 > 0.0 && y1 > 0.0)
    return exp (log (y0) + fraction * log (y1 / y0));
  return y0 + fraction * (y1 - y0);
}

static int
h2_h_collision_selftest (void)
{
  static const double test_temperature[] = { 100.0, 500.0, 5000.0, 7000.0, 10000.0 };
  double worst = 0.0;
  int it, pair;

  for (it = 0; it < (int) (sizeof (test_temperature) / sizeof (test_temperature[0])); it++)
  {
    double temperature = test_temperature[it];
    for (pair = 0; pair < H2_H_PAIRS; pair++)
    {
      H2HRateData *transition = &h2_h_rate_data[pair];
      int upper = transition->upper, lower = transition->lower;
      double down = h2_h_downward_rate (transition, temperature);
      double boltzmann = exp (-(h2_level_data[upper].energy_k - h2_level_data[lower].energy_k) / temperature);
      double up = down * h2_level_data[upper].g / h2_level_data[lower].g * boltzmann;
      double pop_upper = h2_level_data[upper].g * exp (-h2_level_data[upper].energy_k / temperature);
      double pop_lower = h2_level_data[lower].g * exp (-h2_level_data[lower].energy_k / temperature);
      double scale = fmax (fabs (pop_upper * down), fabs (pop_lower * up));
      double residual = scale > 0.0 ? fabs (pop_upper * down - pop_lower * up) / scale : 0.0;
      if (!isfinite (down) || !isfinite (up) || down < 0.0 || up < 0.0)
        return -1;
      if (residual > worst)
        worst = residual;
    }
  }
  if (worst > 1.e-12)
  {
    Error ("Passive H2: H-collision detailed-balance residual %.3e exceeds tolerance\n", worst);
    return -1;
  }
  Log ("Passive H2 H-collision self-test: max LTE pair residual %.3e\n", worst);
  {
    double population[H2_PASSIVE_LEVELS], residual, lte[H2_PASSIVE_LEVELS], partition = 0.0, worst_lte = 0.0;
    int level;
    if (h2_normalized_populations (5000.0, 1.e10, population, &residual) != 0)
      return -1;
    for (level = 0; level < H2_PASSIVE_LEVELS; level++)
    {
      lte[level] = h2_level_data[level].g * exp (-h2_level_data[level].energy_k / 5000.0);
      partition += lte[level];
    }
    for (level = 0; level < H2_PASSIVE_LEVELS; level++)
    {
      double relative;
      lte[level] /= partition;
      relative = fabs (population[level] - lte[level]) / fmax (lte[level], 1.e-300);
      if (relative > worst_lte)
        worst_lte = relative;
    }
    if (residual > 1.e-10 || worst_lte > 1.e-5)
    {
      Error ("Passive H2: stationary matrix test failed residual=%.3e high-density LTE error=%.3e\n", residual, worst_lte);
      return -1;
    }
    Log ("Passive H2 matrix self-test: residual %.3e high-density LTE error %.3e\n", residual, worst_lte);
  }
  return 0;
}

static int
h2_normalized_populations (double temperature, double neutral_h_density, double *population, double *residual)
{
  double *matrix, source[H2_PASSIVE_LEVELS] = { 0.0 }, max_residual = 0.0;
  int pair, row, column, error;

  matrix = calloc (H2_PASSIVE_LEVELS * H2_PASSIVE_LEVELS, sizeof (double));
  if (matrix == NULL)
    return -1;
  for (pair = 0; pair < H2_H_PAIRS; pair++)
  {
    H2HRateData *transition = &h2_h_rate_data[pair];
    int upper = transition->upper, lower = transition->lower;
    double down_collision = neutral_h_density * h2_h_downward_rate (transition, temperature);
    double delta = h2_level_data[upper].energy_k - h2_level_data[lower].energy_k;
    double up = down_collision * h2_level_data[upper].g / h2_level_data[lower].g * exp (-delta / temperature);
    double down = down_collision + transition->spontaneous_rate;
    matrix[lower * H2_PASSIVE_LEVELS + upper] += down;
    matrix[upper * H2_PASSIVE_LEVELS + upper] -= down;
    matrix[upper * H2_PASSIVE_LEVELS + lower] += up;
    matrix[lower * H2_PASSIVE_LEVELS + lower] -= up;
  }
  for (column = 0; column < H2_PASSIVE_LEVELS; column++)
    matrix[(H2_PASSIVE_LEVELS - 1) * H2_PASSIVE_LEVELS + column] = 1.0;
  source[H2_PASSIVE_LEVELS - 1] = 1.0;
  error = solve_matrix (matrix, source, H2_PASSIVE_LEVELS, population, -1);
  free (matrix);
  if (error != 0)
    return -1;
  for (row = 0; row < H2_PASSIVE_LEVELS; row++)
  {
    double balance = 0.0, scale = 0.0;
    if (!isfinite (population[row]) || population[row] < -1.e-12)
      return -1;
    if (population[row] < 0.0)
      population[row] = 0.0;
    for (pair = 0; pair < H2_H_PAIRS; pair++)
    {
      H2HRateData *transition = &h2_h_rate_data[pair];
      int upper = transition->upper, lower = transition->lower;
      double down_collision = neutral_h_density * h2_h_downward_rate (transition, temperature);
      double up = down_collision * h2_level_data[upper].g / h2_level_data[lower].g
        * exp (-(h2_level_data[upper].energy_k - h2_level_data[lower].energy_k) / temperature);
      double down = down_collision + transition->spontaneous_rate;
      if (row == upper)
      {
        balance += population[lower] * up - population[upper] * down;
        scale += population[lower] * up + population[upper] * down;
      }
      else if (row == lower)
      {
        balance += population[upper] * down - population[lower] * up;
        scale += population[upper] * down + population[lower] * up;
      }
    }
    if (scale > 0.0 && fabs (balance) / scale > max_residual)
      max_residual = fabs (balance) / scale;
  }
  *residual = max_residual;
  return 0;
}
