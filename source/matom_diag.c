
/***********************************************************/
/** @file  matom_diag.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  diagnostic macro-atom printouts.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


/**********************************************************/
/** 
 * @brief      a routine which reports the matom level and kpkt emissivities summed
 * 	over all cells. It is called in define_phot() in the spectral cycles, after get_matom_f()
 *
 * @return     simply logs the level emissivites and kpkt emissivites in macro atom mode.
 *
 * @details
 * prints out the level emissivities in macro atoms, summed over each cell.
 *     	Should only be used in macro atom mode in the spectral cycles
 *     	(when geo.matom_radiation = 1).
 *
 * ### Notes ###
 *
 **********************************************************/

int
matom_emiss_report ()
{

  int n, m;
  double emiss_sum, abs_sum, two_photon_emiss_sum, two_photon_fraction;

  two_photon_emiss_sum = 0.0;

  /* Cycle over macro atom levels and log emissivities */

  for (m = 0; m < nlevels_macro; m++)
  {
    emiss_sum = 0.0;
    abs_sum = 0.0;
    for (n = 0; n < NPLASMA; n++)
    {
      emiss_sum += macromain[n].matom_emiss[m];
      abs_sum += macromain[n].matom_abs[m];
      if (two_photon_process_for_upper (m) >= 0)
      {
        two_photon_fraction = f_matom_two_photon_fraction (&plasmamain[n], m, geo.sfmin, geo.sfmax);
        two_photon_emiss_sum += macromain[n].matom_emiss[m] * two_photon_fraction;
      }
    }

    Log ("Macro Atom level emissivities (summed): z %2d i %2d macro %2d n %2d matom_abs %8.4e matom_emiss %8.4e\n",
         xconfig[m].z, xconfig[m].istate, xconfig[m].macro_info, m, abs_sum, emiss_sum);
  }

  Log ("TwoPhotonMatomEmissivity band_min=%.8e band_max=%.8e luminosity=%.8e\n",
       geo.sfmin, geo.sfmax, two_photon_emiss_sum);



  /* Log kpkt emissivities as well */
  emiss_sum = 0.0;
  abs_sum = 0.0;

  for (n = 0; n < NPLASMA; n++)
  {
    emiss_sum += plasmamain[n].kpkt_emiss;
    abs_sum += plasmamain[n].kpkt_abs;

  }

  Log ("Kpkt emissivities (summed over cells): kpkt_abs %8.4e kpkt_emiss %8.4e\n", abs_sum, emiss_sum);

  /* Log totals */
  Log ("Totals: f_matom %le f_kpkt %le\n", geo.f_matom, geo.f_kpkt);

  return (0);
}
