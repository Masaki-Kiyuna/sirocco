/** @file hminus.c
 *
 * Runtime configuration and shared guards for the private H- implementation.
 * The physical network is connected incrementally; all options default off.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "atomic.h"
#include "sirocco.h"

#define HMINUS_THRESHOLD_EV 0.7542
#define HMINUS_RA_MAX_EV 10.95
#define KBOLTZ_EV 8.617333262e-5

typedef struct hminus_table
{
  size_t n;
  double *x;
  double *y;
} HminusTable;

static HminusTable pd_table = { 0, NULL, NULL };
static HminusTable ra_table = { 0, NULL, NULL };
static HminusTable ra_milne_table = { 0, NULL, NULL };

typedef struct hminus_cdf_bank
{
  size_t ntemp;
  double *temperature;
  double *quantile;
  double *knot_cdf;
  double *mean_energy;
} HminusCDFBank;

static HminusCDFBank ra_bank = { 0, NULL, NULL, NULL, NULL };

typedef struct neutral_h_ff_bank
{
  size_t ntemp, nquantile;
  unsigned int nsegment;
  int *segment;
  double *temperature, *power, *mean_energy, *quantile;
} NeutralHFFBank;

static NeutralHFFBank hff_bank = { 0, 0, 0, NULL, NULL, NULL, NULL, NULL };

typedef struct neutral_h_gaunt_table
{
  size_t ngam, nu;
  double gam_min, u_min, step;
  double *log_value;
} NeutralHGauntTable;

static NeutralHGauntTable hff_gaunt = { 0, 0, 0.0, 0.0, 0.0, NULL };

static double
hminus_lte_ratio (double temperature, double g_h_ground)
{
  double log_ratio;

  if (!(temperature > 0.0) || !(g_h_ground > 0.0))
    return 0.0;
  /* n(H-)/[n(H I,1s) ne], with g(H-)=1 and electron spin weight 2. */
  log_ratio = -log (2.0 * g_h_ground)
    + 1.5 * log (PLANCK * PLANCK / (2.0 * PI * MELEC * BOLTZMANN * temperature))
    + HMINUS_THRESHOLD_EV / (KBOLTZ_EV * temperature);
  if (log_ratio > 690.0)
    return exp (690.0);
  if (log_ratio < -690.0)
    return 0.0;
  return exp (log_ratio);
}
static const double ra_knots[HMINUS_RA_KNOTS] = { 1.0, 2.0, 3.4, HMINUS_RA_MAX_EV };

static int hminus_read_csv (const char *path, HminusTable *table);
static double hminus_log_interpolate (const HminusTable *table, double x);
static int hminus_ra_bank_init (const char *root);
static int neutral_h_ff_bank_init (const char *root);
static int neutral_h_ff_gaunt_init (const char *root);
static int neutral_h_ff_kirchhoff_selftest (void);
static void neutral_h_ff_temperature_bracket (double temperature, size_t *lo, size_t *hi, double *upper_weight);
static double neutral_h_ff_cdf_at_node (size_t itemp, double energy_ev);
static double hminus_ra_cdf_at_node (size_t itemp, double energy_ev);
static void hminus_ra_temperature_bracket (double temperature, size_t *lo, size_t *hi, double *upper_weight);

int
hminus_get_runtime_options (void)
{
  char answer[LINELENGTH];

  if (rdpar_keyword_exists ("@Hydrogen.hminus"))
  {
    strcpy (answer, "off");
    modes.hminus_enabled = rdchoice ("@Hydrogen.hminus(off,on)", "0,1", answer);
  }

  if (rdpar_keyword_exists ("@Hydrogen.hminus_ff_mode"))
  {
    strcpy (answer, "john88_guarded");
    modes.hminus_ff_mode =
      rdchoice ("@Hydrogen.hminus_ff_mode(john88_guarded,cloudy_c25_allt)", "0,1", answer);
  }

  if (rdpar_keyword_exists ("@Hydrogen.passive_molecules"))
  {
    strcpy (answer, "off");
    modes.passive_molecules_enabled = rdchoice ("@Hydrogen.passive_molecules(off,on)", "0,1", answer);
  }

  if (rdpar_keyword_exists ("@Hydrogen.passive_molecule_debug_spectra"))
  {
    strcpy (answer, "off");
    modes.passive_molecule_debug_spectra =
      rdchoice ("@Hydrogen.passive_molecule_debug_spectra(off,on)", "0,1", answer);
  }

  return hminus_validate_startup ();
}

int
hminus_validate_startup (void)
{
  if (modes.passive_molecule_debug_spectra && !modes.passive_molecules_enabled)
  {
    Error ("H- setup: passive molecule spectra require @Hydrogen.passive_molecules=on\n");
    return -1;
  }

  if (modes.hminus_enabled)
  {
    if (!(geo.ioniz_mode == IONMODE_MATRIX_BB
          || geo.ioniz_mode == IONMODE_MATRIX_SPECTRALMODEL
          || geo.ioniz_mode == IONMODE_MATRIX_MULTISHOT)
        || geo.rt_mode != RT_MODE_MACRO || geo.macro_simple != FALSE)
    {
      Error ("H- setup: active H- requires matrix_bb, matrix_pow, or matrix_multishot with full macro-atom line transfer\n");
      return -1;
    }
    Log ("Warning: private H- chemical macro-atom implementation enabled. Old windsave restarts are not supported.\n");
    Log ("H- neutral-H free-free mode: %s\n",
         modes.hminus_ff_mode == HMINUS_FF_CLOUDY_C25_ALLT ? "cloudy_c25_allt" : "john88_guarded");
    if (hminus_data_init () != 0)
    {
      return -1;
    }
  }

  if (modes.passive_molecules_enabled)
  {
    Log ("Warning: passive H2/H2+ diagnostics enabled; their populations do not alter the active H-nuclei ledger.\n");
    if (passive_molecules_init () != 0)
    {
      return -1;
    }
  }

  return 0;
}

static int
hminus_read_csv (const char *path, HminusTable *table)
{
  FILE *fp;
  char line[LINELENGTH];
  size_t capacity = 256;
  double x, y;

  fp = fopen (path, "r");
  if (fp == NULL)
  {
    Error ("H- data: unable to open %s\n", path);
    return -1;
  }

  table->x = calloc (capacity, sizeof (*table->x));
  table->y = calloc (capacity, sizeof (*table->y));
  if (table->x == NULL || table->y == NULL)
  {
    Error ("H- data: allocation failed for %s\n", path);
    fclose (fp);
    return -1;
  }

  while (fgets (line, sizeof (line), fp) != NULL)
  {
    if (sscanf (line, "%lf,%lf", &x, &y) != 2)
    {
      continue;
    }
    if (!isfinite (x) || !isfinite (y) || x <= 0.0 || y < 0.0 || (table->n > 0 && x <= table->x[table->n - 1]))
    {
      Error ("H- data: invalid or non-monotonic row %zu in %s\n", table->n + 2, path);
      fclose (fp);
      return -1;
    }
    if (table->n == capacity)
    {
      double *new_x, *new_y;
      capacity *= 2;
      new_x = realloc (table->x, capacity * sizeof (*table->x));
      new_y = realloc (table->y, capacity * sizeof (*table->y));
      if (new_x == NULL || new_y == NULL)
      {
        Error ("H- data: reallocation failed for %s\n", path);
        fclose (fp);
        return -1;
      }
      table->x = new_x;
      table->y = new_y;
    }
    table->x[table->n] = x;
    table->y[table->n] = y;
    table->n++;
  }
  fclose (fp);

  if (table->n < 2)
  {
    Error ("H- data: fewer than two usable rows in %s\n", path);
    return -1;
  }
  return 0;
}

int
hminus_data_init (void)
{
  char path[LINELENGTH];
  const char *root = getenv ("SIROCCO");

  if (pd_table.n > 0 && ra_table.n > 0 && ra_milne_table.n > 0)
  {
    return 0;
  }
  if (root == NULL || root[0] == '\0')
  {
    Error ("H- data: SIROCCO environment variable is not set\n");
    return -1;
  }

  snprintf (path, sizeof (path), "%s/data/hminus/mclaughlin_Hminus_single_detachment.csv", root);
  if (hminus_read_csv (path, &pd_table) != 0)
  {
    return -1;
  }
  snprintf (path, sizeof (path), "%s/data/hminus/mclaughlin_Hminus_attachment.csv", root);
  if (hminus_read_csv (path, &ra_table) != 0)
  {
    return -1;
  }
  snprintf (path, sizeof (path), "%s/data/hminus/Hminus_RA_milne_coefficient.csv", root);
  if (hminus_read_csv (path, &ra_milne_table) != 0)
  {
    return -1;
  }

  if (pd_table.x[0] <= HMINUS_THRESHOLD_EV || pd_table.x[pd_table.n - 1] < 1.e4 || ra_table.x[0] > 20.0
      || ra_table.x[ra_table.n - 1] < 1.e8 || ra_milne_table.n != 368
      || ra_milne_table.x[0] > 20.0 * (1.0 + 1.e-12)
      || ra_milne_table.x[ra_milne_table.n - 1] < 1.e5 * (1.0 - 1.e-12))
  {
    Error ("H- data: supplied tables do not cover the required validation domain: "
           "PD(n=%d, %.17g--%.17g eV), attachment(n=%d, %.17g--%.17g K), "
           "Milne(n=%d, %.17g--%.17g K)\n",
           pd_table.n, pd_table.x[0], pd_table.x[pd_table.n - 1],
           ra_table.n, ra_table.x[0], ra_table.x[ra_table.n - 1],
           ra_milne_table.n, ra_milne_table.x[0], ra_milne_table.x[ra_milne_table.n - 1]);
    return -1;
  }
  if (hminus_ra_bank_init (root) != 0)
  {
    return -1;
  }
  if (neutral_h_ff_bank_init (root) != 0)
  {
    return -1;
  }
  if (neutral_h_ff_gaunt_init (root) != 0)
  {
    return -1;
  }
  if (neutral_h_ff_kirchhoff_selftest () != 0)
  {
    return -1;
  }

  Log ("H- data: McLaughlin PD rows=%zu E=[%.7g,%.7g] eV; RA rows=%zu T=[%.7g,%.7g] K\n",
       pd_table.n, pd_table.x[0], pd_table.x[pd_table.n - 1], ra_table.n, ra_table.x[0], ra_table.x[ra_table.n - 1]);
  Log ("H- data check: sigma(1 eV)=%.8e cm2 sigma(10 eV)=%.8e cm2 "
       "k_RA_Milne(5000 K)=%.8e table/Milne=%.8f\n",
       hminus_pd_sigma (1.0), hminus_pd_sigma (10.0), hminus_ra_coefficient (5000.0),
       hminus_log_interpolate (&ra_table, 5000.0) / hminus_ra_coefficient (5000.0));
  Log ("H- RA CDF: temperatures=%zu quantiles=%d memory=%.3f MiB mean_E(5000K)=%.8e eV\n",
       ra_bank.ntemp, HMINUS_RA_QUANTILES,
       (ra_bank.ntemp * (HMINUS_RA_QUANTILES + HMINUS_RA_KNOTS + 2) * sizeof (double)) / (1024.0 * 1024.0),
       hminus_ra_mean_energy (5000.0));
  Log ("neutral-H ff CDF: mode=%s rows=%zu quantiles=%zu segments=%u memory=%.3f MiB\n",
       modes.hminus_ff_mode == HMINUS_FF_CLOUDY_C25_ALLT ? "cloudy_c25_allt" : "john88_guarded",
       hff_bank.ntemp, hff_bank.nquantile, hff_bank.nsegment,
       (hff_bank.ntemp * hff_bank.nquantile * sizeof (double)) / (1024.0 * 1024.0));
  Log ("neutral-H ff kernel check: pair(T=5000K,E=1eV)=%.12e cm5 "
       "pair(T=5000K,E=5eV)=%.12e cm5 pair(T=20000K,E=1eV)=%.12e cm5\n",
       neutral_h_ff_pair_opacity (5000.0, EV2ERGS / PLANCK),
       neutral_h_ff_pair_opacity (5000.0, 5.0 * EV2ERGS / PLANCK),
       neutral_h_ff_pair_opacity (20000.0, EV2ERGS / PLANCK));
  return 0;
}

void
hminus_data_free (void)
{
  free (pd_table.x);
  free (pd_table.y);
  free (ra_table.x);
  free (ra_table.y);
  free (ra_milne_table.x);
  free (ra_milne_table.y);
  free (ra_bank.temperature);
  free (ra_bank.quantile);
  free (ra_bank.knot_cdf);
  free (ra_bank.mean_energy);
  free (hff_bank.segment);
  free (hff_bank.temperature);
  free (hff_bank.power);
  free (hff_bank.mean_energy);
  free (hff_bank.quantile);
  free (hff_gaunt.log_value);
  pd_table.n = ra_table.n = ra_milne_table.n = 0;
  pd_table.x = pd_table.y = ra_table.x = ra_table.y = NULL;
  ra_milne_table.x = ra_milne_table.y = NULL;
  memset (&ra_bank, 0, sizeof (ra_bank));
  memset (&hff_bank, 0, sizeof (hff_bank));
  memset (&hff_gaunt, 0, sizeof (hff_gaunt));
}

static int
neutral_h_ff_gaunt_init (const char *root)
{
  char path[LINELENGTH], line[LINELENGTH];
  FILE *stream;
  int magic, atomic_number;
  size_t iu, ig;

  snprintf (path, sizeof (path), "%s/data/hminus/cloudy_gauntff_Z01.dat", root);
  stream = fopen (path, "r");
  if (stream == NULL)
  {
    Error ("neutral-H ff: unable to open %s\n", path);
    return -1;
  }
  if (fscanf (stream, "%d", &magic) != 1 || fgets (line, sizeof (line), stream) == NULL
      || fscanf (stream, "%zu %zu", &hff_gaunt.ngam, &hff_gaunt.nu) != 2 || fgets (line, sizeof (line), stream) == NULL
      || fscanf (stream, "%lf", &hff_gaunt.gam_min) != 1 || fgets (line, sizeof (line), stream) == NULL
      || fscanf (stream, "%lf", &hff_gaunt.u_min) != 1 || fgets (line, sizeof (line), stream) == NULL
      || fscanf (stream, "%lf", &hff_gaunt.step) != 1 || fgets (line, sizeof (line), stream) == NULL
      || fscanf (stream, "%d", &atomic_number) != 1 || fgets (line, sizeof (line), stream) == NULL)
  {
    Error ("neutral-H ff: malformed Cloudy gaunt header in %s\n", path);
    fclose (stream);
    return -1;
  }
  if (magic != 20140510 || atomic_number != 1 || hff_gaunt.ngam < 4 || hff_gaunt.nu < 4)
  {
    Error ("neutral-H ff: unsupported Cloudy gaunt table %s\n", path);
    fclose (stream);
    return -1;
  }
  /* Skip complete comment lines without reading a long numerical row into the
     fixed-size line buffer. Cloudy's comments may have leading whitespace. */
  for (;;)
  {
    int c;
    do
      c = fgetc (stream);
    while (c == ' ' || c == '\t' || c == '\r' || c == '\n');
    if (c == '#')
    {
      while (c != '\n' && c != EOF)
        c = fgetc (stream);
      continue;
    }
    if (c == EOF)
    {
      Error ("neutral-H ff: missing gaunt data in %s\n", path);
      fclose (stream);
      return -1;
    }
    ungetc (c, stream);
    break;
  }
  hff_gaunt.log_value = calloc (hff_gaunt.ngam * hff_gaunt.nu, sizeof (*hff_gaunt.log_value));
  if (hff_gaunt.log_value == NULL)
  {
    Error ("neutral-H ff: gaunt allocation failed\n");
    fclose (stream);
    return -1;
  }
  for (iu = 0; iu < hff_gaunt.nu; iu++)
  {
    for (ig = 0; ig < hff_gaunt.ngam; ig++)
    {
      double value;
      if (fscanf (stream, "%lf", &value) != 1 || !(value > 0.0))
      {
        Error ("neutral-H ff: truncated gaunt table at row %zu column %zu\n", iu, ig);
        fclose (stream);
        return -1;
      }
      hff_gaunt.log_value[iu * hff_gaunt.ngam + ig] = log (value);
    }
  }
  fclose (stream);
  return 0;
}

static void
neutral_h_ff_lagrange4 (double x, double x0, double step, double weight[4])
{
  int j, k;
  for (j = 0; j < 4; j++)
  {
    double node_j = x0 + j * step;
    weight[j] = 1.0;
    for (k = 0; k < 4; k++)
      if (k != j)
        weight[j] *= (x - (x0 + k * step)) / (node_j - (x0 + k * step));
  }
}

static double
neutral_h_ff_cloudy_gaunt (double temperature, double energy_ev)
{
  const double ryd_ev = 13.605693122994;
  const double te_1ryd = ryd_ev / KBOLTZ_EV;
  double lg, lu, wg[4], wu[4], along_u[4], log_gaunt;
  long ig, iu;
  int j, k;

  lg = log10 (te_1ryd / temperature);
  lu = log10 (energy_ev / (KBOLTZ_EV * temperature));
  ig = (long) floor ((lg - hff_gaunt.gam_min) / hff_gaunt.step);
  iu = (long) floor ((lu - hff_gaunt.u_min) / hff_gaunt.step);
  if (ig < 1)
    ig = 1;
  if (iu < 1)
    iu = 1;
  if (ig > (long) hff_gaunt.ngam - 3)
    ig = (long) hff_gaunt.ngam - 3;
  if (iu > (long) hff_gaunt.nu - 3)
    iu = (long) hff_gaunt.nu - 3;
  ig--;
  iu--;
  neutral_h_ff_lagrange4 (lg, hff_gaunt.gam_min + ig * hff_gaunt.step, hff_gaunt.step, wg);
  neutral_h_ff_lagrange4 (lu, hff_gaunt.u_min + iu * hff_gaunt.step, hff_gaunt.step, wu);
  for (k = 0; k < 4; k++)
  {
    along_u[k] = 0.0;
    for (j = 0; j < 4; j++)
      along_u[k] += hff_gaunt.log_value[(iu + k) * hff_gaunt.ngam + ig + j] * wg[j];
  }
  log_gaunt = 0.0;
  for (k = 0; k < 4; k++)
    log_gaunt += along_u[k] * wu[k];
  return exp (log_gaunt);
}

static double
neutral_h_ff_john_pair (double temperature, double wavelength_micron)
{
  static const double a[6][6] = {
    {0, 0, 0, 0, 0, 0},
    {2483.346, 285.827, -2054.291, 2827.776, -1341.537, 208.952},
    {-3449.889, -1158.382, 8746.523, -11485.632, 5303.609, -812.939},
    {2200.040, 2427.719, -13651.105, 16755.524, -7510.494, 1132.738},
    {-696.271, -1841.400, 8624.970, -10051.530, 4400.067, -655.020},
    {88.283, 444.517, -1863.864, 2095.288, -901.788, 132.985}
  };
  static const double b[6][6] = {
    {518.1021, -734.8666, 1021.1775, -479.0721, 93.1373, -6.4285},
    {473.2636, 1443.4137, -1977.3395, 922.3575, -178.9275, 12.3600},
    {-482.2089, -737.1616, 1096.8827, -521.1341, 101.7963, -7.0571},
    {115.5291, 169.6374, -245.6490, 114.2430, -21.9972, 1.5097},
    {0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0}
  };
  const double (*coeff)[6] = wavelength_micron >= 0.3645 ? a : b;
  double theta = 5040.0 / temperature;
  double result = 0.0;
  int i;

  for (i = 0; i < 6; i++)
  {
    double inverse = 1.0 / wavelength_micron;
    double polynomial = coeff[i][0] * wavelength_micron * wavelength_micron + coeff[i][1]
      + coeff[i][2] * inverse + coeff[i][3] * inverse * inverse
      + coeff[i][4] * inverse * inverse * inverse + coeff[i][5] * pow (inverse, 4.0);
    result += pow (theta, 0.5 * (i + 2)) * polynomial;
  }
  return fmax (0.0, result * 1.e-29 * BOLTZMANN * temperature);
}

static int
neutral_h_ff_kirchhoff_selftest (void)
{
  size_t representative[3], segment, row;
  double worst = 0.0;

  for (segment = 0; segment < hff_bank.nsegment; segment++)
  {
    representative[segment] = hff_bank.ntemp;
    for (row = 0; row < hff_bank.ntemp; row++)
      if ((unsigned int) hff_bank.segment[row] == segment)
      {
        representative[segment] = row;
        while (row + 1 < hff_bank.ntemp && (unsigned int) hff_bank.segment[row + 1] == segment)
          row++;
        representative[segment] = (representative[segment] + row) / 2;
        break;
      }
  }
  for (segment = 0; segment < hff_bank.nsegment; segment++)
  {
    const size_t itemp = representative[segment];
    const double temperature = hff_bank.temperature[itemp];
    const double *quantile = &hff_bank.quantile[itemp * hff_bank.nquantile];
    const double emin = quantile[0], emax = quantile[hff_bank.nquantile - 1];
    const int intervals = 32768;
    double integral = 0.0, previous = 0.0, dlog;
    int i;

    if (itemp >= hff_bank.ntemp || !(emax > emin) || !(emin > 0.0))
      return -1;
    dlog = log (emax / emin) / intervals;
    for (i = 0; i <= intervals; i++)
    {
      const double energy_ev = emin * exp (i * dlog);
      const double frequency = energy_ev * EV2ERGS / PLANCK;
      const double x = energy_ev / (KBOLTZ_EV * temperature);
      const double planck = x < 700.0 ? 2.0 * PLANCK * frequency * frequency * frequency
        / (VLIGHT * VLIGHT * expm1 (x)) : 0.0;
      const double integrand = 4.0 * PI * neutral_h_ff_pair_opacity (temperature, frequency)
        * planck * frequency;
      if (i > 0)
        integral += 0.5 * (previous + integrand) * dlog;
      previous = integrand;
    }
    {
      const double relative = fabs (integral / hff_bank.power[itemp] - 1.0);
      if (!isfinite (relative) || relative > worst)
        worst = relative;
      Log ("neutral-H ff Kirchhoff check: T=%.8e K direct=%.12e table=%.12e relative=%.3e\n",
           temperature, integral, hff_bank.power[itemp], relative);
    }
  }
  if (!isfinite (worst) || worst > 2.e-3)
  {
    Error ("neutral-H ff Kirchhoff self-test failed: worst relative error %.3e\n", worst);
    return -1;
  }
  return 0;
}

double
neutral_h_ff_pair_opacity (double temperature, double frequency)
{
  const double ryd_ev = 13.605693122994;
  const double ff_abs = 1.036997355477227e-38;
  double energy_ev, wavelength_micron, anu, x, gaunt;

  if (hff_gaunt.log_value == NULL || !(temperature > 0.0) || !(frequency > 0.0))
    return 0.0;
  energy_ev = PLANCK * frequency / EV2ERGS;
  wavelength_micron = 1.2398419843320026 / energy_ev;
  if (modes.hminus_ff_mode == HMINUS_FF_JOHN88_GUARDED && temperature >= 1400.0 && temperature <= 10080.0
      && wavelength_micron >= 0.1823)
    return neutral_h_ff_john_pair (temperature, wavelength_micron);
  anu = energy_ev / ryd_ev;
  x = energy_ev / (KBOLTZ_EV * temperature);
  gaunt = neutral_h_ff_cloudy_gaunt (temperature, energy_ev);
  return ff_abs / sqrt (temperature) * pow (anu, -3.0) * (-expm1 (-x)) * gaunt * sqrt (0.1175 * anu);
}

double
neutral_h_ff_opacity (PlasmaPtr xplasma, double frequency)
{
  int hion, hconfig, ndom;
  double n_h1;

  if (!modes.hminus_enabled || xplasma == NULL || !(xplasma->ne > 0.0))
    return 0.0;
  hion = 0;
  while (hion < nions && !(ion[hion].z == 1 && ion[hion].istate == 1))
    hion++;
  if (hion >= nions)
    return 0.0;
  hconfig = ion[hion].first_nlte_level;
  n_h1 = hconfig >= 0 ? den_config (xplasma, hconfig) : xplasma->density[hion];
  ndom = wmain[xplasma->nwind].ndom;
  return n_h1 * xplasma->ne * neutral_h_ff_pair_opacity (xplasma->t_e, frequency) * zdom[ndom].fill;
}

double
neutral_h_ff_luminosity (PlasmaPtr xplasma, double f1, double f2)
{
  int hion, hconfig;
  double n_h1, e1_ev, e2_ev, luminosity;

  if (!modes.hminus_enabled || xplasma == NULL || !(xplasma->t_e > 0.0) || !(f2 > f1))
    return 0.0;
  hion = 0;
  while (hion < nions && !(ion[hion].z == 1 && ion[hion].istate == 1))
    hion++;
  if (hion >= nions)
    return 0.0;
  hconfig = ion[hion].first_nlte_level;
  n_h1 = hconfig >= 0 ? den_config (xplasma, hconfig) : xplasma->density[hion];
  e1_ev = PLANCK * f1 / EV2ERGS;
  e2_ev = PLANCK * f2 / EV2ERGS;
  luminosity = n_h1 * xplasma->ne * neutral_h_ff_power_coefficient (xplasma->t_e)
    * neutral_h_ff_band_fraction (xplasma->t_e, e1_ev, e2_ev) * xplasma->vol;
  return isfinite (luminosity) && luminosity > 0.0 ? luminosity : 0.0;
}

double
neutral_h_ff_kpkt_cooling (PlasmaPtr xplasma)
{
  int hion, hconfig;
  double n_h1, coefficient;

  if (!modes.hminus_enabled || xplasma == NULL || !(xplasma->t_e > 0.0))
    return 0.0;
  hion = 0;
  while (hion < nions && !(ion[hion].z == 1 && ion[hion].istate == 1))
    hion++;
  if (hion >= nions)
    return 0.0;
  hconfig = ion[hion].first_nlte_level;
  n_h1 = hconfig >= 0 ? den_config (xplasma, hconfig) : xplasma->density[hion];
  coefficient = n_h1 * neutral_h_ff_power_coefficient (xplasma->t_e);
  return isfinite (coefficient) && coefficient > 0.0 ? coefficient : 0.0;
}

static int
neutral_h_ff_bank_init (const char *root)
{
  char path[LINELENGTH], magic[8];
  FILE *stream;
  uint32_t version, nrow, nq, nsegment;
  size_t i;

  snprintf (path, sizeof (path), "%s/data/hminus/neutral_H_ff_%s.bin", root,
            modes.hminus_ff_mode == HMINUS_FF_CLOUDY_C25_ALLT ? "cloudy_c25_allT" : "john88_guarded");
  stream = fopen (path, "rb");
  if (stream == NULL)
  {
    Error ("neutral-H ff: unable to open %s\n", path);
    return -1;
  }
  if (fread (magic, 1, sizeof (magic), stream) != sizeof (magic)
      || memcmp (magic, "HNFFCDF1", sizeof (magic)) != 0
      || fread (&version, sizeof (version), 1, stream) != 1
      || fread (&nrow, sizeof (nrow), 1, stream) != 1
      || fread (&nq, sizeof (nq), 1, stream) != 1
      || fread (&nsegment, sizeof (nsegment), 1, stream) != 1)
  {
    Error ("neutral-H ff: invalid header in %s\n", path);
    fclose (stream);
    return -1;
  }
  if (version != 1 || nrow < 2 || nq < 2 || nsegment < 1 || nsegment > 3)
  {
    Error ("neutral-H ff: unsupported table dimensions in %s\n", path);
    fclose (stream);
    return -1;
  }
  hff_bank.ntemp = nrow;
  hff_bank.nquantile = nq;
  hff_bank.nsegment = nsegment;
  hff_bank.segment = calloc (nrow, sizeof (*hff_bank.segment));
  hff_bank.temperature = calloc (nrow, sizeof (*hff_bank.temperature));
  hff_bank.power = calloc (nrow, sizeof (*hff_bank.power));
  hff_bank.mean_energy = calloc (nrow, sizeof (*hff_bank.mean_energy));
  hff_bank.quantile = calloc ((size_t) nrow * nq, sizeof (*hff_bank.quantile));
  if (hff_bank.segment == NULL || hff_bank.temperature == NULL || hff_bank.power == NULL
      || hff_bank.mean_energy == NULL || hff_bank.quantile == NULL)
  {
    Error ("neutral-H ff: allocation failed for %s\n", path);
    fclose (stream);
    return -1;
  }
  for (i = 0; i < nrow; i++)
  {
    int32_t segment;
    if (fread (&segment, sizeof (segment), 1, stream) != 1
        || fread (&hff_bank.temperature[i], sizeof (double), 1, stream) != 1
        || fread (&hff_bank.power[i], sizeof (double), 1, stream) != 1
        || fread (&hff_bank.mean_energy[i], sizeof (double), 1, stream) != 1
        || fread (&hff_bank.quantile[i * nq], sizeof (double), nq, stream) != nq)
    {
      Error ("neutral-H ff: truncated row %zu in %s\n", i, path);
      fclose (stream);
      return -1;
    }
    hff_bank.segment[i] = segment;
    if (segment < 0 || (unsigned int) segment >= nsegment || !(hff_bank.temperature[i] > 0.0)
        || !(hff_bank.power[i] >= 0.0) || !(hff_bank.mean_energy[i] > 0.0))
    {
      Error ("neutral-H ff: invalid row %zu in %s\n", i, path);
      fclose (stream);
      return -1;
    }
  }
  fclose (stream);
  return 0;
}

static void
neutral_h_ff_temperature_bracket (double temperature, size_t *lo, size_t *hi, double *upper_weight)
{
  int segment;
  size_t first, last, mid;

  if (hff_bank.nsegment == 1)
    segment = 0;
  else if (temperature < 1400.0)
    segment = 0;
  else if (temperature <= 10080.0)
    segment = 1;
  else
    segment = 2;
  first = 0;
  while (first < hff_bank.ntemp && hff_bank.segment[first] != segment)
    first++;
  last = first;
  while (last + 1 < hff_bank.ntemp && hff_bank.segment[last + 1] == segment)
    last++;
  if (temperature <= hff_bank.temperature[first])
  {
    *lo = *hi = first;
    *upper_weight = 0.0;
    return;
  }
  if (temperature >= hff_bank.temperature[last])
  {
    *lo = *hi = last;
    *upper_weight = 0.0;
    return;
  }
  while (last - first > 1)
  {
    mid = (first + last) / 2;
    if (temperature < hff_bank.temperature[mid])
      last = mid;
    else
      first = mid;
  }
  *lo = first;
  *hi = last;
  *upper_weight = (log (temperature) - log (hff_bank.temperature[first]))
    / (log (hff_bank.temperature[last]) - log (hff_bank.temperature[first]));
}

static double
neutral_h_ff_cdf_at_node (size_t itemp, double energy_ev)
{
  const double *q = &hff_bank.quantile[itemp * hff_bank.nquantile];
  size_t lo, hi, mid;

  if (energy_ev <= q[0])
    return 0.0;
  if (energy_ev >= q[hff_bank.nquantile - 1])
    return 1.0;
  lo = 0;
  hi = hff_bank.nquantile - 1;
  while (hi - lo > 1)
  {
    mid = (lo + hi) / 2;
    if (energy_ev < q[mid])
      hi = mid;
    else
      lo = mid;
  }
  return (lo + (energy_ev - q[lo]) / (q[hi] - q[lo])) / (hff_bank.nquantile - 1.0);
}

double
neutral_h_ff_power_coefficient (double temperature)
{
  size_t lo, hi;
  double upper_weight;

  if (hff_bank.ntemp == 0 || !(temperature > 0.0))
    return 0.0;
  neutral_h_ff_temperature_bracket (temperature, &lo, &hi, &upper_weight);
  if (lo == hi)
    return hff_bank.power[lo];
  return exp ((1.0 - upper_weight) * log (hff_bank.power[lo]) + upper_weight * log (hff_bank.power[hi]));
}

double
neutral_h_ff_band_fraction (double temperature, double e1_ev, double e2_ev)
{
  size_t lo, hi;
  double upper_weight, flo, fhi;

  if (hff_bank.ntemp == 0 || !(temperature > 0.0) || !(e2_ev > e1_ev))
    return 0.0;
  neutral_h_ff_temperature_bracket (temperature, &lo, &hi, &upper_weight);
  flo = neutral_h_ff_cdf_at_node (lo, e2_ev) - neutral_h_ff_cdf_at_node (lo, e1_ev);
  if (lo == hi)
    return fmax (0.0, fmin (1.0, flo));
  fhi = neutral_h_ff_cdf_at_node (hi, e2_ev) - neutral_h_ff_cdf_at_node (hi, e1_ev);
  return fmax (0.0, fmin (1.0, (1.0 - upper_weight) * flo + upper_weight * fhi));
}

double
neutral_h_ff_sample_energy (double temperature, double e1_ev, double e2_ev)
{
  size_t lo, hi, itemp, iq;
  double upper_weight, u1, u2, u, x, frac;
  const double *q;

  if (hff_bank.ntemp == 0 || !(e2_ev > e1_ev))
    return -1.0;
  neutral_h_ff_temperature_bracket (temperature, &lo, &hi, &upper_weight);
  itemp = (lo != hi && random_number (0.0, 1.0) < upper_weight) ? hi : lo;
  u1 = neutral_h_ff_cdf_at_node (itemp, e1_ev);
  u2 = neutral_h_ff_cdf_at_node (itemp, e2_ev);
  if (!(u2 > u1))
    return -1.0;
  u = random_number (u1, u2);
  x = u * (hff_bank.nquantile - 1.0);
  iq = (size_t) floor (x);
  if (iq >= hff_bank.nquantile - 1)
    return hff_bank.quantile[itemp * hff_bank.nquantile + hff_bank.nquantile - 1];
  frac = x - iq;
  q = &hff_bank.quantile[itemp * hff_bank.nquantile];
  return (1.0 - frac) * q[iq] + frac * q[iq + 1];
}

static int
hminus_ra_bank_init (const char *root)
{
  char path[LINELENGTH], line[LINELENGTH];
  FILE *fp;
  size_t itemp, ngrid, i, q, k;
  double *energy, *weight, *cumulative;
  double max_milne_relative = 0.0, worst_milne_temperature = 0.0;

  snprintf (path, sizeof (path), "%s/data/hminus/Hminus_RA_adaptive_temperature_nodes.csv", root);
  fp = fopen (path, "r");
  if (fp == NULL)
  {
    Error ("H- RA CDF: unable to open %s\n", path);
    return -1;
  }
  ra_bank.temperature = calloc (HMINUS_RA_TNODES, sizeof (double));
  ra_bank.quantile = calloc ((size_t) HMINUS_RA_TNODES * HMINUS_RA_QUANTILES, sizeof (double));
  ra_bank.knot_cdf = calloc ((size_t) HMINUS_RA_TNODES * HMINUS_RA_KNOTS, sizeof (double));
  ra_bank.mean_energy = calloc (HMINUS_RA_TNODES, sizeof (double));
  if (ra_bank.temperature == NULL || ra_bank.quantile == NULL || ra_bank.knot_cdf == NULL || ra_bank.mean_energy == NULL)
  {
    Error ("H- RA CDF: bank allocation failed\n");
    fclose (fp);
    return -1;
  }
  while (fgets (line, sizeof (line), fp) != NULL && ra_bank.ntemp < HMINUS_RA_TNODES)
  {
    unsigned int index;
    double temperature;
    if (sscanf (line, "%u,%lf", &index, &temperature) == 2)
    {
      if (index != ra_bank.ntemp || !isfinite (temperature) || temperature <= 0.0 ||
          (ra_bank.ntemp > 0 && temperature <= ra_bank.temperature[ra_bank.ntemp - 1]))
      {
        Error ("H- RA CDF: invalid temperature row %u in %s\n", index, path);
        fclose (fp);
        return -1;
      }
      ra_bank.temperature[ra_bank.ntemp++] = temperature;
    }
  }
  fclose (fp);
  if (ra_bank.ntemp != HMINUS_RA_TNODES)
  {
    Error ("H- RA CDF: expected %d temperature nodes but read %zu\n", HMINUS_RA_TNODES, ra_bank.ntemp);
    return -1;
  }

  ngrid = 0;
  while (ngrid < pd_table.n && pd_table.x[ngrid] < HMINUS_RA_MAX_EV)
    ngrid++;
  /* Add the threshold and the exact agreed upper endpoint.  The McLaughlin
   * energy mesh does not in general contain 10.95 eV exactly. */
  ngrid += 2;
  energy = calloc (ngrid, sizeof (double));
  weight = calloc (ngrid, sizeof (double));
  cumulative = calloc (ngrid, sizeof (double));
  if (energy == NULL || weight == NULL || cumulative == NULL)
  {
    Error ("H- RA CDF: integration workspace allocation failed\n");
    return -1;
  }
  energy[0] = HMINUS_THRESHOLD_EV;
  for (i = 1; i < ngrid - 1; i++)
    energy[i] = pd_table.x[i - 1];
  energy[ngrid - 1] = HMINUS_RA_MAX_EV;

  for (itemp = 0; itemp < ra_bank.ntemp; itemp++)
  {
    double total, emission_total, number_total, temperature = ra_bank.temperature[itemp];
    size_t cursor;
    weight[0] = cumulative[0] = 0.0;
    emission_total = number_total = 0.0;
    for (i = 1; i < ngrid; i++)
    {
      double delta = energy[i] - energy[i - 1];
      double midpoint = 0.5 * (energy[i] + energy[i - 1]);
      double quarter1 = 0.5 * (energy[i - 1] + midpoint);
      double quarter3 = 0.5 * (midpoint + energy[i]);
      double middle_number = midpoint * midpoint * hminus_pd_sigma (midpoint) *
        exp (-(midpoint - HMINUS_THRESHOLD_EV) / (KBOLTZ_EV * temperature));
      double quarter1_number = quarter1 * quarter1 * hminus_pd_sigma (quarter1) *
        exp (-(quarter1 - HMINUS_THRESHOLD_EV) / (KBOLTZ_EV * temperature));
      double quarter3_number = quarter3 * quarter3 * hminus_pd_sigma (quarter3) *
        exp (-(quarter3 - HMINUS_THRESHOLD_EV) / (KBOLTZ_EV * temperature));
      weight[i] = pow (energy[i], 3.0) * hminus_pd_sigma (energy[i]) *
        exp (-(energy[i] - HMINUS_THRESHOLD_EV) / (KBOLTZ_EV * temperature));
      cumulative[i] = cumulative[i - 1] + 0.5 * (weight[i - 1] + weight[i]) * delta;
      if (energy[i - 1] > 0.0 && energy[i] > 0.0)
      {
        number_total += delta / 12.0 *
          (weight[i - 1] / energy[i - 1] + 4.0 * quarter1_number + 2.0 * middle_number +
           4.0 * quarter3_number + weight[i] / energy[i]);
        emission_total += delta / 12.0 *
          (weight[i - 1] + 4.0 * quarter1 * quarter1_number + 2.0 * midpoint * middle_number +
           4.0 * quarter3 * quarter3_number + weight[i]);
      }
    }
    total = cumulative[ngrid - 1];
    if (!isfinite (total) || total <= 0.0)
    {
      Error ("H- RA CDF: non-positive normalization at T=%.8e K\n", temperature);
      free (energy);
      free (weight);
      free (cumulative);
      return -1;
    }
    if (!isfinite (number_total) || number_total <= 0.0)
    {
      Error ("H- RA CDF: non-positive photon-number normalization at T=%.8e K\n", temperature);
      free (energy);
      free (weight);
      free (cumulative);
      return -1;
    }
    {
      const double saha0 = 0.25 * pow (PLANCK * PLANCK /
        (2.0 * M_PI * MELEC * BOLTZMANN * temperature), 1.5);
      const double cdf_milne = 8.0 * M_PI * saha0 /
        (VLIGHT * VLIGHT * PLANCK * PLANCK * PLANCK) *
        number_total * EV2ERGS * EV2ERGS * EV2ERGS;
      const double tabulated_milne = hminus_ra_coefficient (temperature);
      const double relative = fabs (cdf_milne - tabulated_milne) /
        fmax (tabulated_milne, 1.e-300);
      if (relative > max_milne_relative)
      {
        max_milne_relative = relative;
        worst_milne_temperature = temperature;
      }
    }
    /* The CDF samples indivisible energy packets.  Cooling needs energy per
     * reaction, i.e. total emissivity divided by photon-number emissivity. */
    ra_bank.mean_energy[itemp] = emission_total / number_total;
    cursor = 1;
    for (q = 0; q < HMINUS_RA_QUANTILES; q++)
    {
      double target = total * q / (HMINUS_RA_QUANTILES - 1.0);
      while (cursor < ngrid - 1 && cumulative[cursor] < target)
        cursor++;
      if (q == 0)
        ra_bank.quantile[itemp * HMINUS_RA_QUANTILES + q] = HMINUS_THRESHOLD_EV;
      else
      {
        double denom = cumulative[cursor] - cumulative[cursor - 1];
        double fraction = denom > 0.0 ? (target - cumulative[cursor - 1]) / denom : 0.0;
        ra_bank.quantile[itemp * HMINUS_RA_QUANTILES + q] = energy[cursor - 1] + fraction * (energy[cursor] - energy[cursor - 1]);
      }
    }
    for (k = 0; k < HMINUS_RA_KNOTS; k++)
    {
      cursor = 1;
      while (cursor < ngrid - 1 && energy[cursor] < ra_knots[k])
        cursor++;
      if (ra_knots[k] >= HMINUS_RA_MAX_EV)
        ra_bank.knot_cdf[itemp * HMINUS_RA_KNOTS + k] = 1.0;
      else
      {
        double fraction = (ra_knots[k] - energy[cursor - 1]) / (energy[cursor] - energy[cursor - 1]);
        double value = cumulative[cursor - 1] + fraction * (cumulative[cursor] - cumulative[cursor - 1]);
        ra_bank.knot_cdf[itemp * HMINUS_RA_KNOTS + k] = value / total;
      }
    }
  }
  free (energy);
  free (weight);
  free (cumulative);
  Log ("H- RA Milne/CDF check: worst_relative=%.6e at T=%.8e K\n",
       max_milne_relative, worst_milne_temperature);
  if (max_milne_relative > 2.e-3)
  {
    Error ("H- RA CDF: photon-number integral is inconsistent with the Milne coefficient table\n");
    return -1;
  }
  return 0;
}

static void
hminus_ra_temperature_bracket (double temperature, size_t *lo, size_t *hi, double *upper_weight)
{
  if (temperature <= ra_bank.temperature[0])
  {
    *lo = *hi = 0;
    *upper_weight = 0.0;
    return;
  }
  if (temperature >= ra_bank.temperature[ra_bank.ntemp - 1])
  {
    *lo = *hi = ra_bank.ntemp - 1;
    *upper_weight = 0.0;
    return;
  }
  *lo = 0;
  *hi = ra_bank.ntemp - 1;
  while (*hi - *lo > 1)
  {
    size_t mid = *lo + (*hi - *lo) / 2;
    if (ra_bank.temperature[mid] <= temperature)
      *lo = mid;
    else
      *hi = mid;
  }
  *upper_weight = log (temperature / ra_bank.temperature[*lo]) /
    log (ra_bank.temperature[*hi] / ra_bank.temperature[*lo]);
}

static double
hminus_ra_cdf_at_node (size_t itemp, double energy_ev)
{
  size_t lo, hi, k;
  const double *quantile = &ra_bank.quantile[itemp * HMINUS_RA_QUANTILES];
  if (energy_ev <= HMINUS_THRESHOLD_EV)
    return 0.0;
  if (energy_ev >= HMINUS_RA_MAX_EV)
    return 1.0;
  for (k = 0; k < HMINUS_RA_KNOTS; k++)
    if (fabs (energy_ev - ra_knots[k]) < 1.e-12 * ra_knots[k])
      return ra_bank.knot_cdf[itemp * HMINUS_RA_KNOTS + k];
  lo = 0;
  hi = HMINUS_RA_QUANTILES - 1;
  while (hi - lo > 1)
  {
    size_t mid = lo + (hi - lo) / 2;
    if (quantile[mid] <= energy_ev)
      lo = mid;
    else
      hi = mid;
  }
  if (quantile[hi] <= quantile[lo])
    return lo / (HMINUS_RA_QUANTILES - 1.0);
  return (lo + (energy_ev - quantile[lo]) / (quantile[hi] - quantile[lo])) / (HMINUS_RA_QUANTILES - 1.0);
}

double
hminus_ra_band_fraction (double temperature, double e1_ev, double e2_ev)
{
  size_t lo, hi;
  double upper_weight, flo, fhi;
  if (ra_bank.ntemp == 0 || e2_ev <= e1_ev)
    return 0.0;
  hminus_ra_temperature_bracket (temperature, &lo, &hi, &upper_weight);
  flo = hminus_ra_cdf_at_node (lo, e2_ev) - hminus_ra_cdf_at_node (lo, e1_ev);
  if (lo == hi)
    return flo;
  fhi = hminus_ra_cdf_at_node (hi, e2_ev) - hminus_ra_cdf_at_node (hi, e1_ev);
  return (1.0 - upper_weight) * flo + upper_weight * fhi;
}

double
hminus_ra_sample_energy (double temperature, double e1_ev, double e2_ev)
{
  size_t lo, hi, itemp, qlo, qhi;
  double upper_weight, cdf1, cdf2, target, qposition, fraction;
  const double *quantile;
  if (ra_bank.ntemp == 0 || e2_ev <= e1_ev)
    return -1.0;
  hminus_ra_temperature_bracket (temperature, &lo, &hi, &upper_weight);
  itemp = (lo != hi && random_number (0.0, 1.0) < upper_weight) ? hi : lo;
  cdf1 = hminus_ra_cdf_at_node (itemp, e1_ev);
  cdf2 = hminus_ra_cdf_at_node (itemp, e2_ev);
  if (cdf2 <= cdf1)
    return -1.0;
  target = cdf1 + random_number (0.0, 1.0) * (cdf2 - cdf1);
  qposition = target * (HMINUS_RA_QUANTILES - 1.0);
  qlo = (size_t) floor (qposition);
  qhi = qlo < HMINUS_RA_QUANTILES - 1 ? qlo + 1 : qlo;
  fraction = qposition - qlo;
  quantile = &ra_bank.quantile[itemp * HMINUS_RA_QUANTILES];
  return quantile[qlo] + fraction * (quantile[qhi] - quantile[qlo]);
}

double
hminus_ra_mean_energy (double temperature)
{
  size_t lo, hi;
  double upper_weight;
  if (ra_bank.ntemp == 0)
    return 0.0;
  hminus_ra_temperature_bracket (temperature, &lo, &hi, &upper_weight);
  if (lo == hi)
    return ra_bank.mean_energy[lo];
  return (1.0 - upper_weight) * ra_bank.mean_energy[lo] + upper_weight * ra_bank.mean_energy[hi];
}

double
hminus_ra_luminosity (PlasmaPtr xplasma, double f1, double f2)
{
  int h1_ion, h1_config, n;
  double e1_ev, e2_ev, event_rate, luminosity;

  if (!modes.hminus_enabled || xplasma == NULL || !(f2 > f1) || !(xplasma->t_e > 0.0))
    return 0.0;
  h1_ion = -1;
  for (n = 0; n < nions; n++)
  {
    if (ion[n].z == 1 && ion[n].istate == 1)
    {
      h1_ion = n;
      break;
    }
  }
  if (h1_ion < 0)
    return 0.0;
  h1_config = ion[h1_ion].first_nlte_level;
  e1_ev = PLANCK * f1 / EV2ERGS;
  e2_ev = PLANCK * f2 / EV2ERGS;
  event_rate = hminus_ra_coefficient (xplasma->t_e)
    * den_config (xplasma, h1_config) * xplasma->ne;
  luminosity = event_rate * hminus_ra_mean_energy (xplasma->t_e) * EV2ERGS
    * xplasma->vol * hminus_ra_band_fraction (xplasma->t_e, e1_ev, e2_ev);
  return isfinite (luminosity) && luminosity > 0.0 ? luminosity : 0.0;
}

double
hminus_ra_kpkt_cooling (PlasmaPtr xplasma)
{
  int n, h1_ion, h1_config;
  double n_h1, excess_energy;

  if (!modes.hminus_enabled || xplasma == NULL || !(xplasma->t_e > 0.0))
    return 0.0;
  h1_ion = -1;
  for (n = 0; n < nions; n++)
  {
    if (ion[n].z == 1 && ion[n].istate == 1)
    {
      h1_ion = n;
      break;
    }
  }
  if (h1_ion < 0)
    return 0.0;
  h1_config = ion[h1_ion].first_nlte_level;
  n_h1 = den_config (xplasma, h1_config);
  excess_energy = fmax (hminus_ra_mean_energy (xplasma->t_e) - HMINUS_THRESHOLD_EV, 0.0) * EV2ERGS;
  return n_h1 * hminus_ra_coefficient (xplasma->t_e) * excess_energy;
}

/**
 * Process a gross H- photodetachment absorption as an indivisible energy
 * packet. The photoelectron excess energy becomes a k-packet with probability
 * (h nu-chi)/h nu. The affinity-energy branch activates the H(1s)+e chemical
 * continuum and terminates through spontaneous RA or an inverse collisional
 * attachment, using the same coefficients as the population matrix.
 *
 * Stimulated RA is included in population and thermal balance, but is not
 * subtracted here: packet transport uses the agreed non-negative gross PD
 * opacity rather than a frequency-dependent net opacity.
 */
int
hminus_pd_event (PhotPtr p, PlasmaPtr xplasma, int *nres, int *which_out)
{
  int n, h1_ion, h1_config;
  double energy_ev, p_thermal, n_hi, n_he1, ne, lte_ratio;
  double radiative_rate, collisional_rate, total_rate, energy_ev_out;

  if (!modes.hminus_enabled || p == NULL || xplasma == NULL || nres == NULL || which_out == NULL)
    return -1;

  xplasma->hminus.pd_events++;
  if (geo.matom_radiation == 1)
  {
    p->w = 0.0;
    *which_out = MATOM;
    return 0;
  }

  energy_ev = PLANCK * p->freq / EV2ERGS;
  p_thermal = energy_ev > 0.0 ? 1.0 - HMINUS_THRESHOLD_EV / energy_ev : 0.0;
  if (p_thermal < 0.0)
    p_thermal = 0.0;
  if (p_thermal > 1.0)
    p_thermal = 1.0;
  if (random_number (0.0, 1.0) < p_thermal)
  {
    macro_gov (p, nres, KPKT, which_out);
    return 0;
  }

  h1_ion = -1;
  n_he1 = 0.0;
  for (n = 0; n < nions; n++)
  {
    if (ion[n].z == 1 && ion[n].istate == 1)
      h1_ion = n;
    if (ion[n].z == 2 && ion[n].istate == 1)
      n_he1 = xplasma->density[n];
  }
  if (h1_ion < 0)
  {
    Error ("hminus_pd_event: H I ground state is unavailable in cell %d\n", xplasma->nplasma);
    p->istat = P_ERROR_MATOM;
    return -1;
  }
  h1_config = ion[h1_ion].first_nlte_level;
  n_hi = xplasma->density[h1_ion];
  ne = xplasma->ne;
  lte_ratio = hminus_lte_ratio (xplasma->t_e, xconfig[h1_config].g);

  radiative_rate = hminus_ra_coefficient (xplasma->t_e) * ne;
  collisional_rate = hminus_detachment_e (xplasma->t_e) * lte_ratio * ne * ne
    + hminus_detachment_h (xplasma->t_e) * lte_ratio * n_hi * ne
    + hminus_detachment_he (xplasma->t_e) * lte_ratio * n_he1 * ne;
  total_rate = radiative_rate + collisional_rate;
  if (!(total_rate > 0.0) || !isfinite (total_rate))
  {
    Error ("hminus_pd_event: no finite chemical-macro exit in cell %d; using RA closure\n", xplasma->nplasma);
    radiative_rate = total_rate = 1.0;
  }

  /* This one-transient-state absorbing matrix has the same normalized final
   * probabilities as an explicit MC jump. Both Matom transition modes thus
   * consume exactly the same chemical edge weights. */
  if (random_number (0.0, 1.0) * total_rate < radiative_rate)
  {
    energy_ev_out = hminus_ra_sample_energy (xplasma->t_e, HMINUS_THRESHOLD_EV, HMINUS_RA_MAX_EV);
    if (!(energy_ev_out > 0.0) || !isfinite (energy_ev_out))
      energy_ev_out = HMINUS_THRESHOLD_EV;
    p->freq = energy_ev_out * EV2ERGS / PLANCK;
    p->emission_process = EMISSION_PROCESS_HMINUS_RA;
    *nres = NRES_HMINUS_PD;
    *which_out = MATOM;
    if (p->origin < 10)
      p->origin += 10;
  }
  else
  {
    macro_gov (p, nres, KPKT, which_out);
  }
  return 0;
}

static double
hminus_log_interpolate (const HminusTable *table, double x)
{
  size_t lo = 0, hi = table->n - 1;
  double fraction;

  while (hi - lo > 1)
  {
    const size_t mid = lo + (hi - lo) / 2;
    if (table->x[mid] <= x)
      lo = mid;
    else
      hi = mid;
  }
  fraction = log (x / table->x[lo]) / log (table->x[hi] / table->x[lo]);
  if (table->y[lo] > 0.0 && table->y[hi] > 0.0)
  {
    return exp (log (table->y[lo]) + fraction * log (table->y[hi] / table->y[lo]));
  }
  return table->y[lo] + fraction * (table->y[hi] - table->y[lo]);
}

double
hminus_pd_sigma (double energy_ev)
{
  if (energy_ev <= HMINUS_THRESHOLD_EV || pd_table.n == 0)
  {
    return 0.0;
  }
  if (energy_ev < pd_table.x[0])
  {
    const double fraction = (energy_ev - HMINUS_THRESHOLD_EV) / (pd_table.x[0] - HMINUS_THRESHOLD_EV);
    return pd_table.y[0] * pow (fraction, 1.5);
  }
  if (energy_ev > pd_table.x[pd_table.n - 1])
  {
    return 0.0;
  }
  if (energy_ev == pd_table.x[pd_table.n - 1])
  {
    return pd_table.y[pd_table.n - 1];
  }
  return hminus_log_interpolate (&pd_table, energy_ev);
}

double
hminus_pd_opacity (PlasmaPtr xplasma, double frequency)
{
  int ndom;
  double energy_ev;

  if (!modes.hminus_enabled || xplasma == NULL || !(frequency > 0.0) || !(xplasma->hminus.density > 0.0))
    return 0.0;
  ndom = wmain[xplasma->nwind].ndom;
  energy_ev = PLANCK * frequency / EV2ERGS;
  return xplasma->hminus.density * hminus_pd_sigma (energy_ev) * zdom[ndom].fill;
}

double
hminus_ra_coefficient (double temperature)
{
  if (ra_milne_table.n == 0 || !isfinite (temperature) || temperature <= 0.0)
  {
    return 0.0;
  }
  if (temperature <= ra_milne_table.x[0])
  {
    return ra_milne_table.y[0];
  }
  if (temperature >= ra_milne_table.x[ra_milne_table.n - 1])
  {
    return ra_milne_table.y[ra_milne_table.n - 1];
  }
  return hminus_log_interpolate (&ra_milne_table, temperature);
}

double
hminus_detachment_e (double temperature)
{
  if (temperature <= 0.0)
    return 0.0;
  return 9.28372e-10 * sqrt (temperature / 300.0) * exp (-21787.8 / temperature);
}

double
hminus_detachment_h (double temperature)
{
  static const double a[10] =
    { -20.37260896, 1.13944933, -0.14210135, 0.00846446, -0.00143276,
    0.00020123, 0.00008664, -0.00002585, 2.4555e-6, -8.0684e-8
  };
  double te, log_te, log_rate, power;
  int i;

  if (temperature <= 0.0)
    return 0.0;
  te = temperature / 11604.51812;
  if (te <= 0.1)
    return 2.5634e-9 * pow (te, 1.78186);
  log_te = log (te);
  log_rate = 0.0;
  power = 1.0;
  for (i = 0; i < 10; i++)
  {
    log_rate += a[i] * power;
    power *= log_te;
  }
  return exp (log_rate);
}

double
hminus_detachment_he (double temperature)
{
  if (temperature <= 0.0)
    return 0.0;
  return 4.1e-17 * temperature * temperature * exp (-19870.0 / temperature);
}

double
hminus_mutual_neutralization (double temperature)
{
  double bounded_temperature;
  if (temperature <= 0.0)
    return 0.0;
  bounded_temperature = temperature;
  if (bounded_temperature < 10.0)
    bounded_temperature = 10.0;
  if (bounded_temperature > 1.e4)
    bounded_temperature = 1.e4;
  return 1.89e-7 * pow (bounded_temperature / 300.0, -0.442);
}

double
hminus_associative_detachment (double temperature)
{
  double numerator, denominator;
  if (temperature <= 0.0)
    return 0.0;
  /* Kreckel/Miller fit used by the agreed passive Launay91 formation
     closure.  This is diagnostic only and is not an active H- sink. */
  numerator = pow (temperature, 0.098493) + 0.32852 * pow (temperature, 0.55610)
    + 2.771e-7 * pow (temperature, 2.1826);
  denominator = 1.0 + 6.191e-3 * pow (temperature, 1.0461)
    + 8.9712e-11 * pow (temperature, 3.0424) + 3.2576e-14 * pow (temperature, 3.7741);
  return 1.35e-9 * numerator / denominator;
}

double
hminus_proton_to_h2plus (double temperature)
{
  if (temperature <= 0.0)
    return 0.0;
  /* GP98_CloudyC25 scalar-ground H2+ diagnostic rate [cm3 s-1]. */
  if (temperature < 7891.0)
    return 6.9e-9 * pow (temperature, -0.35);
  return 9.6e-7 * pow (temperature, -0.9);
}

double
hminus_thermal_cooling (PlasmaPtr xplasma, double temperature)
{
  int n, h1_ion, hp_ion, h1_config, h3_config;
  double n_h1, n_h3, n_hi, n_hp, n_he1, n_hm, ne, volume;
  double lte_ratio, event_ra, mean_ra_energy;
  double det_events, attach_events, mn_events, mn_inverse_events;
  double mn_inverse_coefficient, q_mn, degeneracy_ratio;
  double mn_heat, mn_cool;

  if (!modes.hminus_enabled || !(temperature > 0.0))
    return 0.0;

  h1_ion = -1;
  n_he1 = 0.0;
  for (n = 0; n < nions; n++)
  {
    if (ion[n].z == 1 && ion[n].istate == 1)
      h1_ion = n;
    if (ion[n].z == 2 && ion[n].istate == 1)
      n_he1 = xplasma->density[n];
  }
  if (h1_ion < 0 || h1_ion + 1 >= nions)
    return 0.0;
  hp_ion = h1_ion + 1;
  h1_config = ion[h1_ion].first_nlte_level;
  h3_config = -1;
  for (n = h1_config; n < h1_config + ion[h1_ion].nlte; n++)
  {
    if (xconfig[n].ilv == 3 || fabs (xconfig[n].q_num - 3.0) < 1.0e-6)
    {
      h3_config = n;
      break;
    }
  }
  if (h3_config < 0)
    return 0.0;

  n_h1 = den_config (xplasma, h1_config);
  n_h3 = den_config (xplasma, h3_config);
  n_hi = xplasma->density[h1_ion];
  n_hp = xplasma->density[hp_ion];
  n_hm = xplasma->hminus.density;
  ne = xplasma->ne;
  volume = xplasma->vol;

  event_ra = hminus_ra_coefficient (temperature) * n_h1 * ne;
  mean_ra_energy = hminus_ra_mean_energy (temperature);
  xplasma->hminus.ra_cool = event_ra * fmax (mean_ra_energy - HMINUS_THRESHOLD_EV, 0.0) * EV2ERGS * volume;
  xplasma->hminus.ra_cool += n_h1 * ne
    * fmax (xplasma->hminus.ra_stim_e_rate_old - xplasma->hminus.ra_stim_rate_old, 0.0)
    * HMINUS_THRESHOLD_EV * EV2ERGS * volume;
  xplasma->hminus.ra_lum = event_ra * mean_ra_energy * EV2ERGS * volume;

  det_events = n_hm * (hminus_detachment_e (temperature) * ne
                       + hminus_detachment_h (temperature) * n_hi
                       + hminus_detachment_he (temperature) * n_he1);
  lte_ratio = hminus_lte_ratio (temperature, xconfig[h1_config].g);
  attach_events = n_h1 * (hminus_detachment_e (temperature) * lte_ratio * ne * ne
                          + hminus_detachment_h (temperature) * lte_ratio * n_hi * ne
                          + hminus_detachment_he (temperature) * lte_ratio * n_he1 * ne);
  xplasma->hminus.coll_det_cool = det_events * HMINUS_THRESHOLD_EV * EV2ERGS * volume;
  xplasma->hminus.coll_attach_heat = attach_events * HMINUS_THRESHOLD_EV * EV2ERGS * volume;

  q_mn = ion[h1_ion].ip - HMINUS_THRESHOLD_EV * EV2ERGS - xconfig[h3_config].ex;
  degeneracy_ratio = ion[hp_ion].g / (xconfig[h1_config].g * xconfig[h3_config].g);
  mn_inverse_coefficient = hminus_mutual_neutralization (temperature) * degeneracy_ratio
    * exp (-q_mn / (BOLTZMANN * temperature));
  mn_events = hminus_mutual_neutralization (temperature) * n_hm * n_hp;
  mn_inverse_events = mn_inverse_coefficient * n_h1 * n_h3;
  mn_heat = mn_events * q_mn * volume;
  mn_cool = mn_inverse_events * q_mn * volume;
  xplasma->hminus.mn_heat = mn_heat - mn_cool;

  return xplasma->hminus.ra_cool + xplasma->hminus.coll_det_cool + mn_cool
    - xplasma->hminus.coll_attach_heat - mn_heat;
}

int
hminus_extend_rate_matrix (PlasmaPtr xplasma, double xne, int index_element,
                           double rate_matrix[NLEVELS_MACRO][NLEVELS_MACRO],
                           int conf_to_matrix[NLEVELS_MACRO], int n_macro_lvl)
{
  int h1_ion, hp_ion, h1_config, h3_config, hp_config;
  int h1_matrix, h3_matrix, hp_matrix, hm_matrix, n;
  double n_hi, n_h1, n_h3, n_hp, n_hm_old, n_he1, rate_up, rate_down;
  double lte_ratio, inv_e_coefficient, inv_h_coefficient, inv_he_coefficient, inverse_rate;
  double mn_coefficient, inverse_coefficient, q_mn, degeneracy_ratio;

  if (!modes.hminus_enabled || ion[ele[index_element].firstion].z != 1)
    return n_macro_lvl;
  if (n_macro_lvl >= NLEVELS_MACRO)
  {
    Error ("H- matrix: no room for chemical state in NLEVELS_MACRO=%d\n", NLEVELS_MACRO);
    return -1;
  }

  h1_ion = ele[index_element].firstion;
  hp_ion = h1_ion + 1;
  h1_config = ion[h1_ion].first_nlte_level;
  hp_config = ion[hp_ion].first_nlte_level;
  h1_matrix = conf_to_matrix[h1_config];
  hp_matrix = conf_to_matrix[hp_config];
  h3_config = -1;
  for (n = ion[h1_ion].first_nlte_level; n < ion[h1_ion].first_nlte_level + ion[h1_ion].nlte; n++)
  {
    if (xconfig[n].ilv == 3 || fabs (xconfig[n].q_num - 3.0) < 1.0e-6)
    {
      h3_config = n;
      break;
    }
  }
  if (h1_matrix < 0 || hp_matrix < 0 || h3_config < 0 || conf_to_matrix[h3_config] < 0)
  {
    Error ("H- matrix missing state in cell %d\n", xplasma->nplasma);
    Error ("H- map H1: config %d matrix %d\n", h1_config, h1_matrix);
    Error ("H- map H3: config %d matrix %d\n", h3_config,
           h3_config >= 0 ? conf_to_matrix[h3_config] : -1);
    Error ("H- map Hp: config %d matrix %d; nmacro %d\n", hp_config, hp_matrix, n_macro_lvl);
    return -1;
  }
  h3_matrix = conf_to_matrix[h3_config];
  hm_matrix = n_macro_lvl;
  n_hi = xplasma->density[h1_ion];
  n_h1 = den_config (xplasma, h1_config);
  n_h3 = den_config (xplasma, h3_config);
  n_hp = xplasma->density[hp_ion];
  n_hm_old = xplasma->hminus.density_old > 0.0 ? xplasma->hminus.density_old : xplasma->hminus.density;
  n_he1 = 0.0;
  for (n = 0; n < nions; n++)
  {
    if (ion[n].z == 2 && ion[n].istate == 1)
      n_he1 = xplasma->density[n];
  }

  /* The chemical equilibrium uses spontaneous + induced RA.  Transport
     remains a separate gross-PD process and never uses this as negative
     opacity. */
  rate_up = (hminus_ra_coefficient (xplasma->t_e) + xplasma->hminus.ra_stim_rate_old) * xne;
  rate_down = xplasma->hminus.pd_rate_old;
  xplasma->hminus.det_e_rate = hminus_detachment_e (xplasma->t_e) * xne;
  xplasma->hminus.det_h_rate = hminus_detachment_h (xplasma->t_e) * n_hi;
  xplasma->hminus.det_he_rate = hminus_detachment_he (xplasma->t_e) * n_he1;
  rate_down += xplasma->hminus.det_e_rate + xplasma->hminus.det_h_rate + xplasma->hminus.det_he_rate;
  lte_ratio = hminus_lte_ratio (xplasma->t_e, xconfig[h1_config].g);
  inv_e_coefficient = hminus_detachment_e (xplasma->t_e) * lte_ratio;
  inv_h_coefficient = hminus_detachment_h (xplasma->t_e) * lte_ratio;
  inv_he_coefficient = hminus_detachment_he (xplasma->t_e) * lte_ratio;
  xplasma->hminus.inv_e_rate = inv_e_coefficient * xne * xne;
  xplasma->hminus.inv_h_rate = inv_h_coefficient * n_hi * xne;
  xplasma->hminus.inv_he_rate = inv_he_coefficient * n_he1 * xne;
  inverse_rate = xplasma->hminus.inv_e_rate + xplasma->hminus.inv_h_rate + xplasma->hminus.inv_he_rate;
  xplasma->hminus.ra_rate = rate_up;
  xplasma->hminus.active_sink_rate = rate_down;

  rate_matrix[h1_matrix][h1_matrix] -= rate_up;
  rate_matrix[hm_matrix][h1_matrix] += rate_up;
  rate_matrix[h1_matrix][h1_matrix] -= inverse_rate;
  rate_matrix[hm_matrix][h1_matrix] += inverse_rate;
  rate_matrix[hm_matrix][hm_matrix] -= rate_down;
  rate_matrix[h1_matrix][hm_matrix] += rate_down;

  /*
   * Picard linearisation of H- + H+ <-> H(1s) + H(n=3).  Assigning one
   * product to each reactant column preserves one H nucleus in every matrix
   * column; at a converged Picard fixed point the two column fluxes coincide.
   */
  mn_coefficient = hminus_mutual_neutralization (xplasma->t_e);
  q_mn = ion[h1_ion].ip - HMINUS_THRESHOLD_EV * EV2ERGS - xconfig[h3_config].ex;
  if (!(q_mn > 0.0) || !(ion[hp_ion].g > 0.0) || !(xconfig[h1_config].g > 0.0) || !(xconfig[h3_config].g > 0.0))
  {
    Error ("H- matrix: invalid mutual-neutralization closure Q=%.8e erg g(H+)=%.8e g(H1)=%.8e g(H3)=%.8e\n",
           q_mn, ion[hp_ion].g, xconfig[h1_config].g, xconfig[h3_config].g);
    return -1;
  }
  degeneracy_ratio = ion[hp_ion].g /
    (xconfig[h1_config].g * xconfig[h3_config].g);
  inverse_coefficient = mn_coefficient * degeneracy_ratio * exp (-q_mn / (BOLTZMANN * xplasma->t_e));
  if (!isfinite (inverse_coefficient) || inverse_coefficient < 0.0)
    inverse_coefficient = 0.0;

  rate_matrix[hm_matrix][hm_matrix] -= mn_coefficient * n_hp;
  rate_matrix[h1_matrix][hm_matrix] += mn_coefficient * n_hp;
  rate_matrix[hp_matrix][hp_matrix] -= mn_coefficient * n_hm_old;
  rate_matrix[h3_matrix][hp_matrix] += mn_coefficient * n_hm_old;
  xplasma->hminus.active_sink_rate = rate_down + mn_coefficient * n_hp;

  rate_matrix[h1_matrix][h1_matrix] -= inverse_coefficient * n_h3;
  rate_matrix[hm_matrix][h1_matrix] += inverse_coefficient * n_h3;
  rate_matrix[h3_matrix][h3_matrix] -= inverse_coefficient * n_h1;
  rate_matrix[hp_matrix][h3_matrix] += inverse_coefficient * n_h1;
  xplasma->hminus.mn_rate = mn_coefficient * n_hp;
  xplasma->hminus.mn_inverse_rate = inverse_coefficient * n_h3;

  return n_macro_lvl + 1;
}

int
hminus_copy_population (PlasmaPtr xplasma, int index_element, const double *populations, int hminus_matrix_index)
{
  double total_h;

  if (!modes.hminus_enabled || ion[ele[index_element].firstion].z != 1 || hminus_matrix_index < 0)
    return 0;
  total_h = ele[index_element].abun * xplasma->rho * rho2nh;
  xplasma->hminus.density_old = xplasma->hminus.density;
  xplasma->hminus.density = populations[hminus_matrix_index] * total_h;
  if (!isfinite (xplasma->hminus.density) || xplasma->hminus.density < 0.0)
  {
    Error ("H- matrix: invalid density %.8e in plasma cell %d\n", xplasma->hminus.density, xplasma->nplasma);
    return -1;
  }
  return 0;
}

void
hminus_update_conservation (PlasmaPtr xplasma)
{
  double total_h, accounted_h, net_charge, n_hi = 0.0, n_hp = 0.0;
  double passive_sink, active_sink;
  int h_element, n;

  if (!modes.hminus_enabled)
    return;
  h_element = -1;
  for (n = 0; n < nelements; n++)
  {
    if (ele[n].z == 1)
    {
      h_element = n;
      break;
    }
  }
  if (h_element < 0)
    return;

  total_h = ele[h_element].abun * xplasma->rho * rho2nh;
  for (n = ele[h_element].firstion; n < ele[h_element].firstion + ele[h_element].nions; n++)
  {
    if (ion[n].istate == 1)
      n_hi = xplasma->density[n];
    else if (ion[n].istate == 2)
      n_hp = xplasma->density[n];
  }
  xplasma->hminus.passive_ad_rate = xplasma->hminus.density * n_hi * hminus_associative_detachment (xplasma->t_e);
  xplasma->hminus.passive_h2plus_rate = xplasma->hminus.density * n_hp * hminus_proton_to_h2plus (xplasma->t_e);
  passive_sink = xplasma->hminus.passive_ad_rate + xplasma->hminus.passive_h2plus_rate;
  active_sink = xplasma->hminus.density * xplasma->hminus.active_sink_rate;
  xplasma->hminus.passive_sink_ratio = active_sink > 0.0 ? passive_sink / active_sink : 0.0;
  accounted_h = xplasma->hminus.density;
  net_charge = -xplasma->hminus.density;
  for (n = 0; n < nions; n++)
  {
    net_charge += xplasma->density[n] * (ion[n].istate - 1);
    if (ion[n].z == 1)
      accounted_h += xplasma->density[n];
  }
  xplasma->hminus.h_nuclei_residual = total_h > 0.0 ? (accounted_h - total_h) / total_h : 0.0;
  xplasma->hminus.charge_residual = net_charge - xplasma->ne;
  if (!isfinite (xplasma->hminus.h_nuclei_residual) || !isfinite (xplasma->hminus.charge_residual))
    Error ("H- matrix: non-finite conservation residual in plasma cell %d\n", xplasma->nplasma);
}

void
hminus_reduce_estimators (void)
{
#ifdef MPI_ON
  double *local_double, *global_double;
  unsigned long *local_count, *global_count;
  int n;

  if (!modes.hminus_enabled)
    return;

  local_double = calloc ((size_t) NPLASMA * 8, sizeof (double));
  global_double = calloc ((size_t) NPLASMA * 8, sizeof (double));
  local_count = calloc ((size_t) NPLASMA * 6, sizeof (unsigned long));
  global_count = calloc ((size_t) NPLASMA * 6, sizeof (unsigned long));
  if (local_double == NULL || global_double == NULL || local_count == NULL || global_count == NULL)
  {
    Error ("H- estimator reduction: allocation failed\n");
    Exit (EXIT_FAILURE);
  }

  for (n = 0; n < NPLASMA; n++)
  {
    local_double[8 * n] = plasmamain[n].hminus.pd_rate_est / np_mpi_global;
    local_double[8 * n + 1] = plasmamain[n].hminus.pd_energy_est / np_mpi_global;
    local_double[8 * n + 2] = plasmamain[n].hminus.pd_heat_est / np_mpi_global;
    local_double[8 * n + 3] = plasmamain[n].hminus.ra_stim_rate_est / np_mpi_global;
    local_double[8 * n + 4] = plasmamain[n].hminus.ra_stim_e_rate_est / np_mpi_global;
    local_double[8 * n + 5] = plasmamain[n].hminus.hff_abs_est / np_mpi_global;
    local_double[8 * n + 6] = plasmamain[n].hminus.unresolved_hard_rate / np_mpi_global;
    local_double[8 * n + 7] = plasmamain[n].hminus.unresolved_hard_power / np_mpi_global;
    local_count[6 * n] = plasmamain[n].hminus.pd_packet_paths;
    local_count[6 * n + 1] = plasmamain[n].hminus.pd_events;
    local_count[6 * n + 2] = plasmamain[n].hminus.hff_john_cells;
    local_count[6 * n + 3] = plasmamain[n].hminus.hff_cloudy_cells;
    local_count[6 * n + 4] = plasmamain[n].hminus.cdf_endpoint_uses;
    local_count[6 * n + 5] = plasmamain[n].hminus.negative_net_opacity_count;
  }

  MPI_Allreduce (local_double, global_double, NPLASMA * 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (local_count, global_count, NPLASMA * 6, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  for (n = 0; n < NPLASMA; n++)
  {
    plasmamain[n].hminus.pd_rate_est = global_double[8 * n];
    plasmamain[n].hminus.pd_energy_est = global_double[8 * n + 1];
    plasmamain[n].hminus.pd_heat_est = global_double[8 * n + 2];
    plasmamain[n].hminus.ra_stim_rate_est = global_double[8 * n + 3];
    plasmamain[n].hminus.ra_stim_e_rate_est = global_double[8 * n + 4];
    plasmamain[n].hminus.hff_abs_est = global_double[8 * n + 5];
    /* The neutral-H ff absorption is pure thermal deposition. Keep one MPI
       reduction as the authoritative energy ledger for both diagnostics. */
    plasmamain[n].hminus.hff_heat = global_double[8 * n + 5];
    plasmamain[n].hminus.unresolved_hard_rate = global_double[8 * n + 6];
    plasmamain[n].hminus.unresolved_hard_power = global_double[8 * n + 7];
    plasmamain[n].hminus.pd_packet_paths = global_count[6 * n];
    plasmamain[n].hminus.pd_events = global_count[6 * n + 1];
    plasmamain[n].hminus.hff_john_cells = global_count[6 * n + 2];
    plasmamain[n].hminus.hff_cloudy_cells = global_count[6 * n + 3];
    plasmamain[n].hminus.cdf_endpoint_uses = global_count[6 * n + 4];
    plasmamain[n].hminus.negative_net_opacity_count = global_count[6 * n + 5];
  }

  free (local_double);
  free (global_double);
  free (local_count);
  free (global_count);
#endif
}

void
hminus_reset_cycle_diagnostics (void)
{
  int n;
  if (!modes.hminus_enabled)
    return;
  for (n = 0; n < NPLASMA; n++)
  {
    plasmamain[n].hminus.pd_events = 0;
    plasmamain[n].hminus.cdf_endpoint_uses = 0;
    plasmamain[n].hminus.negative_net_opacity_count = 0;
    plasmamain[n].hminus.unresolved_hard_rate = 0.0;
    plasmamain[n].hminus.unresolved_hard_power = 0.0;
    plasmamain[n].hminus.hff_abs_est = 0.0;
    plasmamain[n].hminus.hff_heat = 0.0;
    plasmamain[n].hminus.picard_iterations = 0;
    plasmamain[n].hminus.picard_iterations_max = 0;
    plasmamain[n].hminus.picard_calls = 0;
    plasmamain[n].hminus.picard_iterations_total = 0;
  }
}

void
hminus_cycle_summary (void)
{
  double max_fraction, max_density, max_pd_rate, max_h_residual, max_charge_residual;
  double sum_pd_heat, sum_ra_cool, sum_coll_det_cool, sum_coll_attach_heat, sum_mn_heat;
  double max_hard_rate, sum_hard_power;
  long long sum_pd_paths, sum_pd_events;
  unsigned int max_picard_iterations;
  unsigned long total_picard_calls, total_picard_iterations;
  int max_cell, n;

  if (!modes.hminus_enabled || rank_global != 0)
    return;

  max_fraction = max_density = max_pd_rate = 0.0;
  max_h_residual = max_charge_residual = 0.0;
  sum_pd_heat = sum_ra_cool = sum_coll_det_cool = sum_coll_attach_heat = sum_mn_heat = 0.0;
  max_hard_rate = sum_hard_power = 0.0;
  sum_pd_paths = sum_pd_events = 0;
  max_picard_iterations = 0;
  total_picard_calls = total_picard_iterations = 0;
  max_cell = -1;
  for (n = 0; n < NPLASMA; n++)
  {
    double total_h, fraction;
    total_h = 0.0;
    if (plasmamain[n].nwind >= 0)
    {
      int hion = 0;
      while (hion < nions && ion[hion].z != 1)
        hion++;
      if (hion < nions)
        total_h = ele[ion[hion].nelem].abun * plasmamain[n].rho * rho2nh;
    }
    fraction = total_h > 0.0 ? plasmamain[n].hminus.density / total_h : 0.0;
    if (fraction > max_fraction)
    {
      max_fraction = fraction;
      max_density = plasmamain[n].hminus.density;
      max_cell = n;
    }
    if (plasmamain[n].hminus.pd_rate_old > max_pd_rate)
      max_pd_rate = plasmamain[n].hminus.pd_rate_old;
    if (fabs (plasmamain[n].hminus.h_nuclei_residual) > max_h_residual)
      max_h_residual = fabs (plasmamain[n].hminus.h_nuclei_residual);
    if (fabs (plasmamain[n].hminus.charge_residual) > max_charge_residual)
      max_charge_residual = fabs (plasmamain[n].hminus.charge_residual);
    sum_pd_heat += plasmamain[n].hminus.pd_heat;
    sum_ra_cool += plasmamain[n].hminus.ra_cool;
    sum_coll_det_cool += plasmamain[n].hminus.coll_det_cool;
    sum_coll_attach_heat += plasmamain[n].hminus.coll_attach_heat;
    sum_mn_heat += plasmamain[n].hminus.mn_heat;
    if (plasmamain[n].hminus.unresolved_hard_rate > max_hard_rate)
      max_hard_rate = plasmamain[n].hminus.unresolved_hard_rate;
    sum_hard_power += plasmamain[n].hminus.unresolved_hard_power;
    if (plasmamain[n].hminus.picard_iterations_max > max_picard_iterations)
      max_picard_iterations = plasmamain[n].hminus.picard_iterations_max;
    total_picard_calls += plasmamain[n].hminus.picard_calls;
    total_picard_iterations += plasmamain[n].hminus.picard_iterations_total;
    sum_pd_paths += plasmamain[n].hminus.pd_packet_paths;
    sum_pd_events += plasmamain[n].hminus.pd_events;
  }
  Log ("H- cycle summary: max_fraction=%.8e density=%.8e cm-3 cell=%d max_PD_rate=%.8e s-1\n",
       max_fraction, max_density, max_cell, max_pd_rate);
  Log ("H- conservation: max_abs_H_residual=%.8e max_abs_charge_residual=%.8e cm-3\n",
       max_h_residual, max_charge_residual);
  Log ("H- thermal ledger [erg/s]: PD_heat=%.8e RA_cool=%.8e coll_det_cool=%.8e coll_attach_heat=%.8e MN_net_heat=%.8e\n",
       sum_pd_heat, sum_ra_cool, sum_coll_det_cool, sum_coll_attach_heat, sum_mn_heat);
  Log ("H- packet diagnostics: PD_path_samples=%lld PD_absorption_events=%lld\n",
       sum_pd_paths, sum_pd_events);
  Log ("H- solver diagnostics: calls=%lu total_iterations=%lu mean_iterations=%.6g max_iterations=%u\n",
       total_picard_calls, total_picard_iterations,
       total_picard_calls > 0 ? (double) total_picard_iterations / total_picard_calls : 0.0,
       max_picard_iterations);
  Log ("H- hard diagnostics: max_exposure_rate=%.8e s-1 per H- hard_exposure_power=%.8e erg/s\n",
       max_hard_rate, sum_hard_power);
}
