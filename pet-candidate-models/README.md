# Candidate PET models — University Square, Paper 1 (UNVERIFIED)

Candidate implementations of the two-node MEMI / PET model used for
Paper 1 (University Square, Melbourne, 6/18/19 January 2018).

**Status: none of these has yet been shown to reproduce the stored
production results** (`Production_PET_Stations1-8_2018Jan.csv`). Do not
cite any of them as "the code behind the paper" until one passes the
test below.

| Folder | What it is | Reproduces stored Tc/Ts? |
|---|---|---|
| `A_thermoregulated_excerpts/` | Gagge-type thermoregulated MEMI (vasomotion + sweating; core set point 36.6 °C), verbatim excerpts from the ChatGPT development threads. Not yet runnable code. | Not yet tested — most likely candidate |
| `B_free_core_two_node/` | Tuned free-core two-node model (M=70, k=8, alpha_sw=0.50), Broyden solver. Basis of the paper's tuning demonstration. | No (median miss ~11 °C) |
| `C_archived_broyden_production/` | Broyden script archived as the production run (M=80, alpha_sw=0.60 defaults; SVF-attenuated shortwave). | No (median miss ~9.8 °C) |

Key clue: stored core temperature lies between 36.29 and 36.55 °C in
every row, which a free-core model (B, C) cannot produce.

## Acceptance test

A candidate is the production model if, fed the stored inputs
(`Ta_C, RH_pct, WS_ms, Kg_Wm2` or `Kg_eff_Wm2, Tg_used_C`), it reproduces
stored `Tc_C`, `Ts_C` and `PET_C` to within 0.05 °C on at least 95% of
rows (excluding `variant == 'recon_A'`).
