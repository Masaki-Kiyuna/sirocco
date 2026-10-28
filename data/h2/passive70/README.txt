Private passive H2/H2+ diagnostic data
======================================

These tables are used only when @Hydrogen.passive_molecules=on.  The resulting
H2/H2+ populations do not feed back on active chemistry, temperature, opacity,
or packet transport.

Source preparation and validation records are retained under:
  /work/kiyunams/sirocco_run/ssh_handoff_hminus_h2_h2plus_20261027/
  /work/kiyunams/sirocco_run/ssh_update_h2_h_public_20261027/

Files and SHA256
----------------
4bb6d681ab392a533e35aed9d1e20e501184a05a9bccc47c2b7c768608bf3b1d  cloudy_C25_H2_UV_dissociation.csv
153f35bf716db9a52bafab6f8f60de5d6d4d0c586c433154df33dc11f9174dfd  cloudy_C25_H2_UV_lines.csv
4dd3c321584987a9bccc3a6e8fba432ad15922518d4bc6e58244bf6500423ab1  downward_rates70.csv
6ae0ef2f7c3f54bbb34b353a6044cef0c3c43ee7ba72498c84d735fd1dd882b5  hminus_deposit.dat
cc54ca7907d29404b6c0e0eebfe11924bf88c13abf2442e924f9f4a4bc127d59  levels70.csv
a99d3a2ccff93c128a2b9cb84194ab1d65340b91a42b5f5d3fc7667cc00de820  roueff_H2_IR_transitions.csv

source/passive_molecules.c validates row counts, level identities, pair
coverage, non-negative rates, endpoint interpolation, detailed balance, and a
70-level matrix test before a run starts.
