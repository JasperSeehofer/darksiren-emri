# M192 — the paper's wording: two rulings for the author

Context: M193 withdrew the "M_z halves the H0 error" headline. M192 item 1 cleared four
thesis-era remnants and restated the study's aim in `appendix_parameters.tex` as ruling 3's
question. Two items are the author's to rule ([RULE], C-12).

## Q1 — the title (`paper/main.tex:26`)

Current: *"Constraints on the Hubble constant from EMRI dark sirens with LISA using the massive black hole mass"*

The title promises an H0 constraint from the mass, which ruling 1 withdrew and ruling 3 reframed
("what LISA EMRIs teach us", with H0 as one output).

- **Suggested:** *"What LISA EMRI dark sirens measure: the H0 information budget once population, selection and the black-hole–galaxy relation are marginalised"*
- Shorter: *"The H0 information in LISA EMRI dark sirens, conditional and marginalised"*
- Keep the current title (not recommended: it restates the withdrawn framing).

`paper/PAPER-CONFIG.json` follows whatever title is chosen.

## Q2 — P12, the Reines–Volonteri scatter wording

- `paper/sections/discussion.tex:115–116`: *"The Reines & Volonteri relation has an intrinsic scatter of approximately 0.55 dex, with best-fit parameters α = 7.45 ± 0.08 and β = 1.05 ± 0.11."*
- `discussion.tex:175–176` (same wording, one ruling covers both): *"Moreover, the intrinsic scatter of the relation (~0.55 dex) is the dominant source of uncertainty in the galaxy-side mass estimate."*

The issue: 0.55 dex is the scatter RV15 report for their AGN sample. Calling it "intrinsic" overstates
what is known; the published low-mass relations differ in normalisation by more than 1 dex, and the
Rank-1 budget treats both normalisation (prior 0.5 dex) and scatter (prior 0.2 dex) as nuisances.

- **Suggested:** say "a scatter of about 0.55 dex (RV15, AGN sample)" and add that the normalisation
  and scatter are uncertain at low mass and are marginalised in the Rank-1 budget.
- Keep as is.
