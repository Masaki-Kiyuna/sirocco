#ifndef SIROCCO_TOY_DIAG_H
#define SIROCCO_TOY_DIAG_H
int toy_diag_options (void);
int toy_diag_enabled (void);
double h1_omukai2001_deexc_coeff_e (struct lines *line_ptr, double t);
double h1_omukai2001_deexc_coeff_h (struct lines *line_ptr, double t);
void toy_diag_path (PlasmaPtr xp, PhotPtr p, double ds);
void toy_diag_bf (PlasmaPtr xp, int lower, double freq, double ft, double y);
void toy_diag_jbar (PlasmaPtr xp, int index, double raw, double stimfac);
void toy_diag_es (PhotPtr p, int grid);
void toy_diag_begin (void);
void toy_diag_cell (PlasmaPtr xp, int phase);
void toy_diag_end (void);
void toy_diag_commit (void);
double toy_diag_escape_used (struct lines *line_ptr, PlasmaPtr xp, double beta);
#endif
