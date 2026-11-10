# M193 R1 — REPORT

## Budget
- Wall time: first `date` 14:55:30 CET, report written about 15:00 CET (2026-11-10). Tool calls: about 14. No subagents, no Workflow.
- Least-confident rows (read these first):
  1. **P5** (nearest grammatical form). The survey replaces only the tail of the sentence "The comparison between the GW-measured $\Mz/(1+z)$ and the catalog-predicted $\Mbh$ … not exploited in Ref.~[Laghi]". I replaced the whole sentence with the survey sentence ("Weighting EMRI host candidates by such a relation …"), which now follows "…empirical scaling relations~[Reines]" directly, so "such a relation" has its antecedent. Keeping the old sentence's opening would have needed new wording about the mass comparison's power, which the page forbids.
  2. **P3**. The page's `\pending{…}` has no closing full stop; written as given, so it renders "[RESULT PENDING: …] We discuss the impact …".
  3. **P11**. The `\cite{Reines:2015moa}` and `[Eq.~\eqref{eq:Lcat_mass}]` inside the removed sentence go with it; both stay cited/labelled elsewhere (integrity check below).
  4. **P8** keeps the line-30–32 sentence ("Second, we introduce the MBH mass channel, which exploits the precise measurement of $\Mz$ … as an additional discriminant for host galaxy identification.") unchanged, as the page's row starts at line 32.
- Nothing touched method.tex, a constant, a model choice or an equation other than the two thesis-result equations (P14/P15). No STOP condition hit.

## Rows (all edits in commit `a41bf60c`)
| # | file:line (main 8ed76e73) | After | Commit |
|---|---|---|---|
| P1 | main.tex:26, PAPER-CONFIG.json | unchanged (title on M192's card) | — |
| P2 | main.tex:45–48 | "LISA measures $\Mz$ to a fractional precision far below the scatter of any black-hole--galaxy relation, so the mass channel can at best re-weight candidate hosts within that scatter." | a41bf60c |
| P3 | main.tex:49–54 | Gray sentence stands alone ("We use …"); + ruling-1 sentence + `\pending{Rank-1: every $\HH$ uncertainty conditional on and marginalised …}` | a41bf60c |
| P4 | introduction.tex:66–74 | Laghi corrected: Fisher errors from Babak+17, mock complete above 10^10 Msun to z≲0.5; 1.6--3.6% (90%), 2.5% M1, 1.1--2.6% for 10 yr | a41bf60c |
| P5 | introduction.tex:85–88 | "Weighting … proposed as future work by Laghi et al. \cite[Sec.~6.4]{Laghi:2021pqk}, and Zhu et al. \cite{Zhu:2024qpp} … using it at least to estimate catalogue completeness; neither quantifies …" | a41bf60c |
| P6 | introduction.tex:96 | "We compare two likelihood configurations" | a41bf60c |
| P7 | introduction.tex:101–104 | ruling-1 sentence + `\pending{Rank-1: what the mass channel adds …}` | a41bf60c |
| P8 | discussion.tex:32–34 | "The idea is not new: Laghi … Sec.~6.4 … Zhu \cite{Zhu:2024qpp} … (0.4 dex, window [-3σ,+2σ]), used at least for completeness and without an ablation … What this work adds is that ablation … on a real catalogue (GLADE+) \pending{…}" | a41bf60c |
| P9 | discussion.tex:72 | \subsection{What the MBH mass can and cannot add} (label kept) | a41bf60c |
| P10 | discussion.tex:76–79 | "…limited by the scatter of the black-hole--galaxy relation, not by the precision of $\Mz$. The thesis-era factor ~2 is withdrawn. \pending{…}" | a41bf60c |
| P11 | discussion.tex:94–103 | keeps the 10^-6 clause and its \ref; "the width of the mass filter is set by the scatter of the Mbh--Mstellar relation (see below), not by the GW measurement. \pending{Rank-1 census …}" | a41bf60c |
| P12 | discussion.tex:115–116 | unchanged (M192 card) | — |
| P13 | conclusions.tex:8–17 | "The thesis-era answer, a factor ~2, is withdrawn: no surviving result reproduces it. … \pending{Rank-1 summary}." | a41bf60c |
| P14 | results.tex:130–139 | equation removed; `\pending{Rank-1 / production value for h^(no mass)}` + withdrawal sentence | a41bf60c |
| P15 | results.tex:141–150 | as P14 for h^(mass) | a41bf60c |
| P16 | results.tex:152–168 | \paragraph{Mass-channel ablation.} \pending{Rank-1: σ_h with and without the mass term …} | a41bf60c |
| P17 | results.tex:179–182 | caption: last sentence and both Eq. refs removed | a41bf60c |
| P18 | results.tex:327–329 | \pending{Rank-1 census of candidates per event with and without the mass window, GLADE+ today}. | a41bf60c |
| P19 | results.tex:369–379 | thesis sub-header, its three rows and the \hline before it removed | a41bf60c |
| P20 | discussion.tex:5–9 | "We now place the analysis in the context of prior work, discuss what the MBH mass can and cannot add, and assess the limitations of the analysis." | a41bf60c |
| P21 | discussion.tex:21–24 | "First, Laghi et al. took Fisher-matrix errors from … \cite{Babak:2017tow}; we compute the full Cramér--Rao bound … for every simulated event (Sec.~\ref{sec:fisher}), which captures …" | a41bf60c |
| P22 | discussion.tex:105–111 | "Any gain from the mass channel should be largest …"; "… the mass information adds less because …" | a41bf60c |
| P23 | results.tex:186–193, 199–200, 205–206 | paragraph → \pending{Rank-1: single-event likelihoods … (Fig.~\ref{fig:single_event}).}; todo "with MBH mass."; caption "additionally using $\Mz$." | a41bf60c |
| P24 | results.tex:295–298 | \pending{Show the effect of the completeness correction on the $\HH$ posterior bias (production run, $h_\mathrm{true} = 0.73$).} | a41bf60c |
| P25 | results.tex:215–218, 229–230 | text → \pending{Rank-1: posterior width against the number of events, per likelihood channel.}; caption ends "…expected $N^{-1/2}$ scaling.}" | a41bf60c |

Also: CHANGELOG `[Unreleased]` "### Changed (paper, 2026-11-10)"; `paper/references.bib` + `Zhu:2024qpp`. `appendix_parameters.tex` unchanged.

## Done-criterion grep (`/usr/bin/grep`, after merging main)
0 lines (main: 45). Per file: all seven paper files 0.

## Integrity check (one-off stdlib script in a temp dir, not committed)
```
labels=46 refs=51 cites=89 bibkeys=26
unresolved refs: none
missing cite keys: none
main.tex: begin=3 end=3 braces 80/80 OK
sections/appendix_parameters.tex: begin=2 end=2 braces 24/24 OK
sections/conclusions.tex: begin=0 end=0 braces 7/7 OK
sections/discussion.tex: begin=0 end=0 braces 89/89 OK
sections/introduction.tex: begin=0 end=0 braces 46/46 OK
sections/method.tex: begin=15 end=15 braces 306/306 OK
sections/results.tex: begin=10 end=10 braces 185/185 OK
INTEGRITY: PASS
```

## Diff names (`git diff --name-only $(git merge-base HEAD main) HEAD`, before the REPORT commit)
```
.planning/missions/M193-darksiren-paper-headline-withdrawal/inspire-2403.04950.bib
CHANGELOG.md
paper/main.tex
paper/references.bib
paper/sections/conclusions.tex
paper/sections/discussion.tex
paper/sections/introduction.tex
paper/sections/results.tex
```
The REPORT commit adds only `REPORT-R1.md` and `publish-round-r1/**` under the mission folder. Mission files are git-ignored (`.planning/*`) and added with `git add -f`, as M162's were.

## Tree hash
`git rev-parse HEAD:darksiren_emri` = `71acb98b0fd8a0f9d3aad9595eb9fc084cae0b6d` (before and after the merge of main; main is still `8ed76e73`, merge "Already up to date").

## INSPIRE
GET `https://inspirehep.net/api/arxiv/2403.04950?format=bibtex` succeeded: texkey `Zhu:2024qpp`, Zhu, Fan, Chen, Hu & Zhang, ApJS 273(2), 24, doi 10.3847/1538-4365/ad5446. Raw response saved as `inspire-2403.04950.bib` and appended verbatim. A second fetch matched it after whitespace normalisation. Nothing is left open.

## Digit check on added paper lines
Every digit token in the added lines comes from the literature quotations (Laghi: 10^10, 0.5, 1.6–3.6, 90, 2.5, 1.1–2.6, Sec. 6.4; Zhu: 0.4, -3/+2 σ), the existing context (Gray PRD 101 122001 2020, 14×14, 10^-6, h_true 0.73, N^-1/2), the withdrawn "factor ~2" itself, or the bib record.

## Left for M192 (not edited, on purpose)
- `results.tex:50` "(thesis baseline: $N_\mathrm{det}^\mathrm{post} = 223$)"
- `results.tex` systematics "$\Delta h \sim 0.02$" (main :314–315)
- `results.tex` "were obtained with the detection probability $\pdet$ set to unity" (main :239–240) and `discussion.tex:230` "use" the $\pdet = 1$ baseline that P14, P15 and P19 withdraw
- P12 scatter wording (`discussion.tex` 0.55 dex, and :175–176); title (P1); appendix aim (`appendix_parameters.tex:23–26`)
- Literal "tighter mass filter" in the discussion paragraph after P12 (main :120). It is not in the grep and not in scope, and it describes a calibration improvement, not the headline.
