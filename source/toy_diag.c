/** Diagnostic-only toy-model ledgers. No checkpoint structure or RNG changes. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef MPI_ON
#include <mpi.h>
#endif
#include "atomic.h"
#include "sirocco.h"
#include "toy_diag.h"

#define TB_MAX 280
#define TB_FIELDS 49
#define TB_PATH 0
#define TB_NUMBER 1
#define TB_RADIAL 2
#define TB_FF 3
#define TB_HFF 4
#define TB_ES 5
#define TB_ESPOWER 6
#define TB_PATHCOUNT 7
#define TB_ESWEIGHT2 8
#define TB_GAMMA 9
#define TB_GAMMAE 29

static int enabled = -1, nb = 0, conf[21], ready = 0;
static double edges[TB_MAX + 1], *bank = NULL, *jraw = NULL, *jstim = NULL;
static double *matrix_beta = NULL;
static int *matrix_cycle = NULL;
static FILE *fcell, *fbb, *fbf, *fspec, *fhmb, *fheat, *fcont, *fpop;

int
toy_diag_enabled (void)
{
  if (enabled < 0)
    enabled = getenv ("SIROCCO_TOY_DIAG") && atoi (getenv ("SIROCCO_TOY_DIAG")) != 0;
  return enabled;
}

int
toy_diag_options (void)
{
  char answer[32];
  toy_diag_enabled ();
  if (rdpar_keyword_exists ("@Diag.toy_model"))
  {
    strcpy (answer, enabled ? "on" : "off");
    enabled = rdchoice ("@Diag.toy_model(off,on)", "0,1", answer);
  }
  if (enabled)
    Log ("Toy diagnostic: n=1--20 gross rates, raw J, frequency-resolved path moments; no physics changes.\n");
  return 0;
}

static int
cmp_double (const void *a, const void *b)
{
  double x = *(const double *) a, y = *(const double *) b;
  return x < y ? -1 : x > y;
}

static void
init (void)
{
  int i, n, count = 0;
  if (ready || !toy_diag_enabled () || NPLASMA <= 0)
    return;
  for (i = 0; i <= 240; i++)
    edges[count++] = pow (10., -2. + 6. * i / 240.) * EV2ERGS / PLANCK;
  for (n = 0; n <= 20; n++)
    conf[n] = -1;
  for (i = 0; i < nlevels_macro; i++)
    if (xconfig[i].z == 1 && xconfig[i].istate == 1)
    {
      n = xconfig[i].ilv;
      if (n >= 1 && n <= 20)
        conf[n] = i;
    }
  for (i = 0; i < nphot_total; i++)
    if (phot_top[i].macro_info && phot_top[i].nlev >= 0 &&
        xconfig[phot_top[i].nlev].z == 1 && xconfig[phot_top[i].nlev].istate == 1 &&
        phot_top[i].freq[0] > edges[0] && phot_top[i].freq[0] < edges[240])
      edges[count++] = phot_top[i].freq[0];
  qsort (edges, count, sizeof (double), cmp_double);
  nb = 0;
  for (i = 1; i < count; i++)
    if (edges[i] > edges[nb] * (1. + 1.e-12))
      edges[++nb] = edges[i];
  bank = calloc ((size_t) NPLASMA * nb * TB_FIELDS, sizeof (double));
  jraw = calloc ((size_t) NPLASMA * size_Jbar_est, sizeof (double));
  jstim = calloc ((size_t) NPLASMA * size_Jbar_est, sizeof (double));
  matrix_beta = calloc ((size_t) NPLASMA * nlines, sizeof (double));
  matrix_cycle = calloc ((size_t) NPLASMA * nlines, sizeof (int));
  if (!bank || !jraw || !jstim || !matrix_beta || !matrix_cycle || nb >= TB_MAX)
  {
    Error ("Toy diagnostic allocation or bin count failed\n");
    Exit (EXIT_FAILURE);
  }
  ready = 1;
  Log ("Toy diagnostic memory %.3f MiB/rank; %d frequency bins\n",
       ((double) NPLASMA * (nb * TB_FIELDS + 2 * size_Jbar_est) * sizeof (double)) / 1048576., nb);
}

double
toy_diag_escape_used (struct lines *line_ptr, PlasmaPtr xp, double beta)
{
  size_t idx;
  if (!toy_diag_enabled ()) return beta;
  init ();
  idx = (size_t) xp->nplasma*nlines + (line_ptr-line);
  matrix_beta[idx] = beta;
  matrix_cycle[idx] = geo.wcycle+1;
  return beta;
}

static int
bin_of (double freq)
{
  int lo = 0, hi = nb, mid;
  if (freq < edges[0] || freq >= edges[nb])
    return -1;
  while (hi - lo > 1)
  {
    mid = (hi + lo) / 2;
    if (freq >= edges[mid]) lo = mid;
    else hi = mid;
  }
  return lo;
}

static double *
at (int cell, int bin)
{
  return bank + ((size_t) cell * nb + bin) * TB_FIELDS;
}

void
toy_diag_path (PlasmaPtr xp, PhotPtr p, double ds)
{
  int b;
  double *v, r, mu, wd;
  if (!toy_diag_enabled () || geo.ioniz_or_extract != CYCLE_IONIZ)
    return;
  init ();
  b = bin_of (p->freq);
  if (b < 0) return;
  v = at (xp->nplasma, b);
  wd = p->w * ds;
  r = length (p->x);
  mu = r > 0 ? dot (p->x, p->lmn) / r : 0;
  v[TB_PATH] += wd;
  v[TB_NUMBER] += wd / (PLANCK * p->freq);
  v[TB_RADIAL] += wd * mu / (PLANCK * p->freq);
  v[TB_FF] += wd * kappa_ff (xp, p->freq);
  if (modes.hminus_enabled)
    v[TB_HFF] += wd * neutral_h_ff_opacity (xp, p->freq);
  v[TB_PATHCOUNT] += 1;
}

void
toy_diag_bf (PlasmaPtr xp, int lower, double freq, double ft, double y)
{
  int n, b;
  double *v;
  if (!toy_diag_enabled () || geo.ioniz_or_extract != CYCLE_IONIZ)
    return;
  if (xconfig[lower].z != 1 || xconfig[lower].istate != 1) return;
  init ();
  n = xconfig[lower].ilv;
  b = bin_of (freq);
  if (n < 1 || n > 20 || b < 0) return;
  v = at (xp->nplasma, b);
  v[TB_GAMMA+n-1] += y / freq;
  v[TB_GAMMAE+n-1] += y / ft;
}

void
toy_diag_es (PhotPtr p, int grid)
{
  struct photon local;
  int b, cell;
  double *v;
  if (!toy_diag_enabled () || geo.ioniz_or_extract != CYCLE_IONIZ) return;
  init ();
  cell = wmain[grid].nplasma;
  if (cell < 0 || cell >= NPLASMA) return;
  if (p->frame == F_OBSERVER)
    observer_to_local_frame (p, &local);
  else
    stuff_phot (p, &local);
  b = bin_of (local.freq);
  if (b < 0) return;
  v = at (cell, b);
  v[TB_ES] += 1;
  v[TB_ESPOWER] += local.w;
  v[TB_ESWEIGHT2] += local.w * local.w;
}

void
toy_diag_jbar (PlasmaPtr xp, int index, double raw, double stimfac)
{
  if (!toy_diag_enabled ()) return;
  init ();
  jraw[(size_t) xp->nplasma * size_Jbar_est + index] = raw;
  jstim[(size_t) xp->nplasma * size_Jbar_est + index] = stimfac;
}

static FILE *
open_table (const char *kind, const char *columns)
{
  char name[1024];
  FILE *fp;
  snprintf (name, sizeof name, "diag_%s/%s.toy.latest.rank%04d.%s.tsv.tmp", files.root, files.root,
            rank_global, kind);
  fp = fopen (name, "w");
  if (!fp) { Error ("Toy diagnostic cannot open %s\n", name); Exit (EXIT_FAILURE); }
  fprintf (fp, "# phase=0 transport state; phase=1 post-update state. Radiation estimators always from the just-finished transport.\n");
  fprintf (fp, "# Rates cm^-3 s^-1, powers erg cm^-3 s^-1 unless named cell; Jnu cgs. phase-lag terms are explicitly retained.\n");
  fprintf (fp, "cycle rank phase cell r_cm depth_cm Te_K ne_cm3 volume_cm3 %s\n", columns);
  return fp;
}

void
toy_diag_begin (void)
{
  size_t count;
  int cell, bin, field;
  if (!toy_diag_enabled ()) return;
  init ();
  count = (size_t) NPLASMA * nb * TB_FIELDS;
#ifdef MPI_ON
  MPI_Allreduce (MPI_IN_PLACE, bank, (int) count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  /* SIROCCO launches a full-luminosity packet ensemble on each rank and
     averages energy/rate estimators, while actual packet counters are summed. */
  for(cell=0;cell<NPLASMA;cell++) for(bin=0;bin<nb;bin++) for(field=0;field<TB_FIELDS;field++)
    if(field!=TB_ES && field!=TB_PATHCOUNT)
      at(cell,bin)[field] /= field==TB_ESWEIGHT2 ? (double)np_mpi_global*np_mpi_global : np_mpi_global;
#endif
  fcell = open_table ("cells", "nHI nHII nHminus nH_total rho dvds_ave rinner_cm router_cm depth_inner_cm depth_outer_cm");
fbb = open_table ("bb", "lower_n upper_n line_index energy_erg nlower nupper A_s B_lu B_ul tau_sob beta_used beta_uncached J_raw J_used stimfac_transport abs_MC stim_current stim_applied spontaneous_bare spontaneous_escape trapped_eliminated rad_net_solver coll_e_up coll_e_down coll_HI_up coll_HI_down q_e_up q_e_down q_HI_up q_HI_down electron_cool electron_heat Hcoll_cool Hcoll_heat rad_abs_power_MC rad_stim_power_applied rad_emit_power_bare rad_emit_power_escape rawJ_available J_local_closure abs_local_closure stim_local_closure closure_valid beta_last_matrix matrix_cycle_valid");
  fbf = open_table ("bf", "level_n continuum_index threshold_Hz nlower nHII gamma gamma_e alpha_sp alpha_sp_e alpha_st alpha_st_e ion_photo rec_sp rec_st ion_coll rec_3body photo_heat photo_abs_power rr_thermal_cool rr_emit_power stim_thermal_cool coll_ion_cool threebody_heat bf_net_into_HI alpha_sp_current alpha_sp_e_current");
  fspec = open_table ("radiation", "bin elo_eV ehi_eV nlevel photon_fluence_4piJ_over_hnu radial_photon_flux Jnu_binmean gamma_bin gamma_e_bin effective_sigma photo_rate photo_heat photo_abs_power ff_abs_cell hff_abs_cell es_events es_weight es_weight2 paths_number rr_emit_cell rr_thermal_cell");
  fhmb = open_table ("chemistry", "process from_n to_n event_rate electron_heat electron_cool radiation_emit radiation_absorb active");
  fheat = open_table ("thermal", "heat_total_cell cool_total_cell photo_heat_cell ff_heat_cell ff_lum_cell rr_cool_cell rr_lum_cell bb_heat_cell bb_cool_cell bf_cool_cell comp_heat_cell comp_cool_cell induced_comp_cell adiabatic_cell shock_cell hminus_pd_heat_cell hminus_ra_cool_cell hminus_det_cool_cell hminus_attach_heat_cell hminus_mn_net_heat_cell hff_heat_cell hff_cool_cell passive_H2_formation passive_H2_cooling_cell");
  fcont = open_table ("continuum", "bin elo_eV ehi_eV mean_Jnu radial_photon_flux es_events es_power_weight es_weight2 ff_abs_cell ff_emit_Kirchhoff_cell hff_abs_cell hff_emit_Kirchhoff_cell two_photon_emit_cell tau_es_local rho_rad_E flux_E");
  fpop = open_table ("populations", "level_n density fraction_H gross_in gross_out net_into residual_over_max gross_rad_in gross_rad_out gross_coll_in gross_coll_out gross_bf_in gross_bf_out gross_chem_in gross_chem_out");
}

static void
prefix (FILE *fp, PlasmaPtr xp, int phase)
{
  WindPtr w = &wmain[xp->nwind];
  double rin = wmain[zdom[w->ndom].nstart].r;
  fprintf (fp, "%d %d %d %d %.17e %.17e %.17e %.17e %.17e ", geo.wcycle+1,
           rank_global, phase, xp->nwind, w->rcen, w->rcen-rin, xp->t_e, xp->ne, xp->vol);
}

static void
chemical_row (PlasmaPtr xp, int phase, const char *name, int from, int to, double rate,
              double heat, double cool, double emit, double absorb)
{
  prefix (fhmb, xp, phase);
  fprintf (fhmb, "%s %d %d %.17e %.17e %.17e %.17e %.17e %d\n", name, from, to,
           rate, heat, cool, emit, absorb, modes.hminus_enabled);
}

void
toy_diag_cell (PlasmaPtr xp, int phase)
{
  int i, j, k, l, u, n, idx, iph, bfd;
  double nl, nu, b, tau, bu, bd, jr, js, sf, de, ce, ch, eu, ed, hu, hd, up;
  double jlocal, denom, emitted, planck, dnu, two, dr, rrbin, rrheat, fnu, df, alpha;
  double h0, hp, hm, h1, h3, g, ge, rs, rse, st, ste, ph, rr, sr, ci, tb, energy;
  double inv, a0, ae, q, event, vol, *v, fluence, ga, gae, ft, norm;
  LinePtr line_ptr;
  TopPhotPtr cp;
  MacroPtr mp = &macromain[xp->nplasma];
  WindPtr w = &wmain[xp->nwind];
  double pop[23]={0}, rin[23]={0}, rout[23]={0}, cin[23]={0}, cout[23]={0};
  double fin[23]={0}, fout[23]={0}, hin[23]={0}, hout[23]={0}, tin, tout;
  if (!toy_diag_enabled () || w->inwind != W_ALL_INWIND) return;
  h0 = xp->density[0]; hp = xp->density[1]; hm = xp->hminus.density;
  pop[21]=hp; pop[22]=hm;
  for(n=1;n<=20;n++) if(conf[n]>=0) pop[n]=den_config(xp,conf[n]);
  vol = xp->vol;
  prefix (fcell, xp, phase);
  fprintf (fcell, "%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
           h0, hp, hm, h0+hp+hm, xp->rho, w->dvds_ave, w->r, wmain[xp->nwind+1].r,
           w->r-wmain[zdom[w->ndom].nstart].r, wmain[xp->nwind+1].r-wmain[zdom[w->ndom].nstart].r);
  for (i = 0; i < nlines; i++)
  {
    line_ptr = &line[i];
    if (!line_ptr->macro_info || line_ptr->z != 1 || line_ptr->istate != 1) continue;
    l = line_ptr->levl; u = line_ptr->levu;
    if (l < 1 || u > 20 || l >= u) continue;
    nl = den_config (xp, line_ptr->nconfigl); nu = den_config (xp, line_ptr->nconfigu);
    idx = -1;
    for (j = 0; j < xconfig[line_ptr->nconfigl].n_bbu_jump; j++)
      if (xconfig[line_ptr->nconfigl].bbu_jump[j] == i)
        idx = xconfig[line_ptr->nconfigl].bbu_indx_first+j;
    if (idx < 0) continue;
    jr = jraw[(size_t) xp->nplasma * size_Jbar_est+idx];
    sf = jstim[(size_t) xp->nplasma * size_Jbar_est+idx];
    js = mp->jbar_old[idx];
    tau = sobolev (w, w->x, xp->density[line_ptr->nion], line_ptr, w->dvds_ave);
    b = line_escape_factor (line_ptr, xp);
    bu = b12 (line_ptr); bd = bu * line_ptr->gl / line_ptr->gu;
    energy = PLANCK * line_ptr->freq;
    de = macro_bb_coll_deexc_rate (line_ptr, xp);
    up = macro_bb_coll_exc_rate (line_ptr, xp);
    ch = 0;
    if (h1_omukai2001_line_applies (line_ptr))
    {
      ce = h1_omukai2001_deexc_coeff_e (line_ptr, xp->t_e) * xp->ne;
      ch = h1_omukai2001_deexc_coeff_h (line_ptr, xp->t_e) * h0;
    }
    else ce = xp->ne * q21 (line_ptr, xp->t_e);
    ed = nu * ce; hd = nu * ch;
    eu = nl * ce * line_ptr->gu / line_ptr->gl * exp (-energy / BOLTZMANN / xp->t_e);
    hu = nl * ch * line_ptr->gu / line_ptr->gl * exp (-energy / BOLTZMANN / xp->t_e);
    /* Match the standard upward rate exactly, including atomic collision data choices. */
    if (!h1_omukai2001_line_applies (line_ptr)) eu = nl * up;
    rin[u]+=nl*bu*jr; rout[l]+=nl*bu*jr;
    rin[l]+=nl*bu*jr*(1-sf)+nu*a21(line_ptr)*b;
    rout[u]+=nl*bu*jr*(1-sf)+nu*a21(line_ptr)*b;
    cin[u]+=eu+hu; cout[l]+=eu+hu; cin[l]+=ed+hd; cout[u]+=ed+hd;
    prefix (fbb, xp, phase);
    fprintf (fbb, "%d %d %d", l, u, i);
#define D(x) fprintf (fbb, " %.17e", (double) (x))
    D(energy); D(nl); D(nu); D(a21(line_ptr)); D(bu); D(bd); D(tau); D(b); D(p_escape_from_tau(tau));
    D(jr); D(js); D(sf); D(nl*bu*jr); D(nu*bd*jr); D(nl*bu*jr*(1-sf));
    D(nu*a21(line_ptr)); D(nu*a21(line_ptr)*b); D(nu*a21(line_ptr)*(1-b));
    D(nl*bu*js-nu*a21(line_ptr)*b); D(eu); D(ed); D(hu); D(hd);
    D(xp->ne>0 ? (nl>0 ? eu/nl/xp->ne : 0):0); D(xp->ne>0 ? ce/xp->ne:0);
    D(h0>0 ? (nl>0 ? hu/nl/h0:0):0); D(h0>0 ? ch/h0:0);
    D(eu*energy); D(ed*energy); D(hu*energy); D(hd*energy);
    D(nl*bu*jr*energy); D(nl*bu*jr*(1-sf)*energy); D(nu*a21(line_ptr)*energy); D(nu*a21(line_ptr)*b*energy); D(1);
    denom=nl*bu-nu*bd;
    jlocal=denom>0 ? jr+(1-b)*nu*a21(line_ptr)/denom : 0;
    D(jlocal);D(nl*bu*jlocal);D(nu*bd*jlocal);D(denom>0);
    D(matrix_beta[(size_t)xp->nplasma*nlines+i]);
    D(matrix_cycle[(size_t)xp->nplasma*nlines+i]==geo.wcycle+1);
    fprintf (fbb, "\n");
#undef D
    (void) de;
  }
  for (n = 1; n <= 20; n++)
  {
    if (conf[n] < 0) continue;
    nl = den_config (xp, conf[n]);
    for (j = 0; j < xconfig[conf[n]].n_bfu_jump; j++)
    {
      iph = xconfig[conf[n]].bfu_jump[j]; cp = &phot_top[iph];
      idx = xconfig[conf[n]].bfu_indx_first+j;
      g=mp->gamma_old[idx]; ge=mp->gamma_e_old[idx]; st=mp->alpha_st_old[idx]; ste=mp->alpha_st_e_old[idx];
      bfd = -1;
      for (k=0;k<xconfig[cp->uplev].n_bfd_jump;k++)
        if (xconfig[cp->uplev].bfd_jump[k]==iph) bfd=xconfig[cp->uplev].bfd_indx_first+k;
      rs = bfd>=0 ? mp->recomb_sp[bfd]:0;
      rse = bfd>=0 ? mp->recomb_sp_e[bfd]:0;
      energy=PLANCK*cp->freq[0];
      ph=nl*g; rr=hp*xp->ne*rs; sr=hp*xp->ne*st;
      ci=nl*xp->ne*q_ioniz(cp,xp->t_e); tb=hp*xp->ne*xp->ne*q_recomb(cp,xp->t_e);
      fin[n]+=rr+sr+tb; fout[21]+=rr+sr+tb;
      fout[n]+=ph+ci; fin[21]+=ph+ci;
      a0=alpha_sp(cp,xp,0); ae=alpha_sp(cp,xp,2);
      prefix (fbf,xp,phase);
      fprintf(fbf,"%d %d",n,iph);
#define D(x) fprintf(fbf," %.17e",(double)(x))
      D(cp->freq[0]);D(nl);D(hp);D(g);D(ge);D(rs);D(rse);D(st);D(ste);
      D(ph);D(rr);D(sr);D(ci);D(tb);D(nl*(ge-g)*energy);D(nl*ge*energy);
      D(hp*xp->ne*rse*energy);D(hp*xp->ne*(rs+rse)*energy);
      D(hp*xp->ne*(ste-st)*energy);D(ci*energy);D(tb*energy);D(rr+sr+tb-ph-ci);D(a0);D(ae);
      fprintf(fbf,"\n");
#undef D
    }
  }
  if (conf[2]>=0)
  {
    i=two_photon_process_for_upper(conf[2]);
    if (i>=0)
    {
      event=den_config(xp,conf[2])*two_photon_rate(i);
      rin[1]+=event; rout[2]+=event;
      chemical_row(xp,phase,"two_photon",2,1,event,0,0,event*(xconfig[conf[2]].ex-xconfig[conf[1]].ex),0);
    }
  }
  if (modes.hminus_enabled && conf[1]>=0 && conf[3]>=0)
  {
    HminusCell *h=&xp->hminus;
    double det[3]={h->det_e_rate,h->det_h_rate,h->det_he_rate};
    double att[3]={h->inv_e_rate,h->inv_h_rate,h->inv_he_rate};
    const char *dn[3]={"detach_e","detach_HI","detach_He"};
    const char *an[3]={"attach_e","attach_HI","attach_He"};
    h1=den_config(xp,conf[1]); h3=den_config(xp,conf[3]);
    chemical_row(xp,phase,"photodetach_used",-1,1,hm*h->pd_rate_old,h->pd_heat/vol,0,0,h->pd_abs_power/vol);
    hin[1]+=hm*h->pd_rate_old; hout[22]+=hm*h->pd_rate_old;
    chemical_row(xp,phase,"photodetach_lagged",-1,1,hm*h->pd_rate,0,0,0,0);
    event=h1*xp->ne*hminus_ra_coefficient(xp->t_e);
    hout[1]+=event;hin[22]+=event;
    q=hminus_ra_mean_energy(xp->t_e)*EV2ERGS;
    chemical_row(xp,phase,"RA_spontaneous",1,-1,event,0,event*(q-0.7542*EV2ERGS),event*q,0);
    event=h1*xp->ne*h->ra_stim_rate_old;
    hout[1]+=event;hin[22]+=event;
    chemical_row(xp,phase,"RA_stimulated",1,-1,event,0,h1*xp->ne*(h->ra_stim_e_rate_old-h->ra_stim_rate_old)*0.7542*EV2ERGS,0,0);
    for (i=0;i<3;i++)
    {
      chemical_row(xp,phase,dn[i],-1,1,hm*det[i],0,hm*det[i]*0.7542*EV2ERGS,0,0);
      chemical_row(xp,phase,an[i],1,-1,h1*att[i],h1*att[i]*0.7542*EV2ERGS,0,0,0);
      hin[1]+=hm*det[i];hout[22]+=hm*det[i];hout[1]+=h1*att[i];hin[22]+=h1*att[i];
    }
    q=ion[0].ip-0.7542*EV2ERGS-xconfig[conf[3]].ex;
    event=hm*h->mn_rate;
    hin[1]+=event;hin[3]+=event;hout[21]+=event;hout[22]+=event;
    chemical_row(xp,phase,"MN_Hminus_HII_to_H1_H3",-1,3,event,event*q,0,0,0);
    event=h1*h->mn_inverse_rate;
    hout[1]+=event;hout[3]+=event;hin[21]+=event;hin[22]+=event;
    chemical_row(xp,phase,"MN_inverse_H1_H3_to_Hminus_HII",3,-1,event,0,event*q,0,0);
    (void) h3;
  }
  for(n=1;n<=22;n++)
  {
    tin=rin[n]+cin[n]+fin[n]+hin[n];tout=rout[n]+cout[n]+fout[n]+hout[n];
    prefix(fpop,xp,phase);
    fprintf(fpop,"%d %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
      n==22?-1:n,pop[n],pop[n]/(h0+hp+hm),tin,tout,tin-tout,(tin-tout)/fmax(fmax(tin,tout),1.e-300),
      rin[n],rout[n],cin[n],cout[n],fin[n],fout[n],hin[n],hout[n]);
  }
  norm=w->vol/w->xgamma_cen;
  for (i=0;i<nb;i++)
  {
    v=at(xp->nplasma,i); ft=sqrt(edges[i]*edges[i+1]);
    fluence=v[TB_NUMBER]/norm;
    if (phase==0)
    {
      dnu=edges[i+1]-edges[i];
      planck=H_OVER_K*ft/xp->t_e<700 ? 2*PLANCK*ft*ft*ft/VLIGHT/VLIGHT/expm1(H_OVER_K*ft/xp->t_e) : 0;
      emitted=4*PI*kappa_ff(xp,ft)*planck*dnu*vol;
      two=0;
      if(conf[2]>=0 && (j=two_photon_process_for_upper(conf[2]))>=0)
        two=den_config(xp,conf[2])*two_photon_rate(j)*(xconfig[conf[2]].ex-xconfig[conf[1]].ex)*vol*two_photon_band_fraction(j,edges[i],edges[i+1]);
      dr=wmain[xp->nwind+1].r-w->r;
      prefix(fcont,xp,phase);
      fprintf(fcont,"%d",i);
#define D(x) fprintf(fcont," %.17e",(double)(x))
      D(edges[i]*PLANCK/EV2ERGS);D(edges[i+1]*PLANCK/EV2ERGS);D(v[TB_PATH]/norm/4/PI/dnu);
      D(v[TB_RADIAL]/norm);D(v[TB_ES]);D(v[TB_ESPOWER]);D(v[TB_ESWEIGHT2]);D(v[TB_FF]);D(emitted);
      D(v[TB_HFF]);D(modes.hminus_enabled?4*PI*neutral_h_ff_opacity(xp,ft)*planck*dnu*vol:0);D(two);
      D(xp->ne*THOMPSON*dr);D(v[TB_PATH]/norm/VLIGHT);D(v[TB_RADIAL]/norm*PLANCK*ft);
      fprintf(fcont,"\n");
#undef D
    }
    if(phase!=0)continue;
    for(n=1;n<=20;n++)
    {
      if(conf[n]<0)continue;
      ga=v[TB_GAMMA+n-1]/PLANCK/norm; gae=v[TB_GAMMAE+n-1]/PLANCK/norm;
      nl=den_config(xp,conf[n]);
      energy=xconfig[conf[n]].n_bfu_jump>0 ? PLANCK*phot_top[xconfig[conf[n]].bfu_jump[0]].freq[0]:0;
      rrbin=rrheat=0;
      for(j=0;j<xconfig[conf[n]].n_bfu_jump;j++)
      {
        cp=&phot_top[xconfig[conf[n]].bfu_jump[j]];
        if(edges[i+1]<=cp->freq[0])continue;
        df=(edges[i+1]-fmax(edges[i],cp->freq[0]))/8;
        for(k=0;k<8;k++)
        {
          fnu=fmax(edges[i],cp->freq[0])+(k+0.5)*df;
          alpha=5.79618e-36*xconfig[conf[n]].g/xconfig[cp->uplev].g*pow(xp->t_e,-1.5)
            *sigma_phot(cp,fnu)*fnu*fnu*exp(H_OVER_K*(cp->freq[0]-fnu)/xp->t_e)*df;
          rrbin+=hp*xp->ne*alpha*PLANCK*fnu*vol;
          rrheat+=hp*xp->ne*alpha*PLANCK*(fnu-cp->freq[0])*vol;
        }
      }
      prefix(fspec,xp,phase);
      fprintf(fspec,"%d %.17e %.17e %d",i,edges[i]*PLANCK/EV2ERGS,edges[i+1]*PLANCK/EV2ERGS,n);
#define D(x) fprintf(fspec," %.17e",(double)(x))
      D(fluence);D(v[TB_RADIAL]/norm);D(v[TB_PATH]/norm/4/PI/(edges[i+1]-edges[i]));
      D(ga);D(gae);D(fluence>0?ga/fluence:0);D(nl*ga);D(nl*(gae-ga)*energy);D(nl*gae*energy);
      D(v[TB_FF]);D(v[TB_HFF]);D(v[TB_ES]);D(v[TB_ESPOWER]);D(v[TB_ESWEIGHT2]);D(v[TB_PATHCOUNT]);
      D(rrbin);D(rrheat);
      fprintf(fspec,"\n");
#undef D
    }
    (void) ft;
  }
  prefix(fheat,xp,phase);
  fprintf(fheat,"%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
          xp->heat_tot,xp->cool_tot,xp->heat_photo,xp->heat_ff,xp->lum_ff,xp->cool_rr,xp->lum_rr,
          xp->heat_lines_macro,xp->cool_lines_macro,xp->cool_bf_macro,xp->heat_comp,xp->cool_comp,
          xp->heat_ind_comp,xp->cool_adiabatic,xp->heat_shock,xp->hminus.pd_heat,xp->hminus.ra_cool,
          xp->hminus.coll_det_cool,xp->hminus.coll_attach_heat,xp->hminus.mn_heat,
          xp->hminus.hff_heat,xp->hminus.hff_cool,xp->molecules.rate_formation,xp->molecules.cooling_line);
  (void) inv;
}

void
toy_diag_end (void)
{
  if (!toy_diag_enabled ()) return;
  fclose(fcell);fclose(fbb);fclose(fbf);fclose(fspec);fclose(fhmb);fclose(fheat);
  fclose(fcont);
  fclose(fpop);
  memset(bank,0,(size_t)NPLASMA*nb*TB_FIELDS*sizeof(double));
  memset(jraw,0,(size_t)NPLASMA*size_Jbar_est*sizeof(double));
  memset(jstim,0,(size_t)NPLASMA*size_Jbar_est*sizeof(double));
}

void
toy_diag_commit (void)
{
  const char *kinds[]={"cells","bb","bf","radiation","chemistry","thermal","continuum","populations"};
  char pending[2048], complete[1024];
  FILE *fp;
  int i;
  if (!toy_diag_enabled ()) return;
#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  for(i=0;i<8;i++)
  {
    snprintf(complete,sizeof complete,"diag_%s/%s.toy.latest.rank%04d.%s.tsv",files.root,files.root,rank_global,kinds[i]);
    snprintf(pending,sizeof pending,"%s.tmp",complete);
    if(rename(pending,complete)!=0)
    {Error("Toy diagnostic commit failed: %s\n",complete);Exit(EXIT_FAILURE);}
  }
#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if(rank_global==0)
  {
    snprintf(complete,sizeof complete,"diag_%s/%s.toy.latest.complete",files.root,files.root);
    snprintf(pending,sizeof pending,"%s.tmp",complete);
    fp=fopen(pending,"w");
    if(!fp){Error("Toy diagnostic commit marker open failed\n");Exit(EXIT_FAILURE);}
    fprintf(fp,"saved_cycle %d\nranks %d\n",geo.wcycle,np_mpi_global);
    fclose(fp);
    if(rename(pending,complete)!=0){Error("Toy diagnostic marker commit failed\n");Exit(EXIT_FAILURE);}
  }
}
