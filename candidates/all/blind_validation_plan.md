# Multi-season blind validation plan

Written 2026-09-20, before any label exists for the decks listed below. The point of writing it first is that the scope hypothesis in section 3 was tuned on labels, so it can only be tested on labels it has never seen.

## 1. Decks

Built by `scripts/make_blind_deck.py`. Each deck is 30 episodes: 18 drawn from tier 1 and 12 from tier 2, round-robin across days so the largest days cannot fill it, then shuffled. Cutouts are rendered without the model probability into `analysis/blind_cutouts/`. Tier and probability per episode are in `blind_key_<name>_SEALED.csv` and stay closed until every row of that deck is labeled.

| deck | season | seed | direction | note |
|---|---|---|---|---|
| 2024 | tail2024 | 2024 | parallel | solar maximum, largest catalog |
| 2012b | tail2012 | 2012 | anti-parallel | the 15 episodes of the earlier non-blind deck are excluded |
| 2013 | tail2013 | 2013 | parallel | coincidence check untestable, labels are the only evidence |
| 2018 | tail2018 | 2018 | anti-parallel | only 14 tier-1 episodes exist, all are drawn, so this deck has 26 episodes; coincidence untestable |
| 2019 | tail2019 | 2019 | anti-parallel | coincidence untestable |
| 2026 | tail2026 | 2026 | parallel | solar maximum |

Label with the vocabulary and rules in `candidates/LABELING.md`: `positive`, `negative`, or `borderline`, plus `label_confidence` and `label_note`. Suggested order: 2024, 2012b, 2013, 2018, 2019, 2026.

## 2. Labels that already exist and are not blind

`tier1_deck_2011_labeled.csv` and `tier1_deck_2012_labeled.csv` hold 15 tier-1 episodes each, labeled 2026-09-19 with the model probability visible in every cutout title. They are a tuning set, not a validation set.

| season | positive | negative | borderline | contamination among decided |
|---|---|---|---|---|
| 2011 | 7 | 5 | 3 | 5/12 = 42% |
| 2012 | 11 | 3 | 1 | 3/14 = 21% |
| both | 18 | 8 | 4 | 8/26 = 31% |

Against 1/14 = 7% in the 2015 blind deck. The label notes attribute nearly every negative to warm electrons or a plasma-sheet edge, not to morphology.

## 3. Hypotheses, fixed before labeling

- **H1, the current catalog.** Tier-1 contamination at the 0.82 cut, per deck and pooled, by episode and by independent structure. A deck is a concern if its point estimate exceeds 20%.
- **H2, a tighter scope.** On the 26 decided tuning labels, `te_med < 125 eV` with `beta_med < 0.01` keeps 17 of 18 positives and 1 of 8 negatives. The prediction is that on the new decks this scope brings pooled tier-1 contamination to 10% or below while keeping at least 85% of labeled positives. Known cost: in 2015 it removes 13 tier-1 episodes, 9 labeled positive and none negative, so it may be the wrong shape for seasons that already validate. If H2 fails either half of the prediction, the scope is not adopted and context goes back to being a labeling consideration only.
- **H3, separation.** AUC of model probability against verdict over decided episodes of both tiers, per deck and pooled, against 0.975 for the 2015 blind deck under the same model.

No other scope values are to be tried on these decks. If a different cut looks better after unsealing, it needs a fresh deck.

## 4. Scoring

`python scripts/score_blind_decks.py` refuses to open the key of any deck that still has unlabeled rows. `--scope 125 0.01` adds the H2 numbers. Structures are same-probe episodes within 480 s of each other, the same tau `link_events.py` uses, and a structure counts as positive if any member is. Intervals are Clopper-Pearson.

An intra-rater repeat deck, two episodes from each deck above, is to be built after the six are labeled and labeled last.

## 5. Extension to 2011 through 2016

Written 2026-09-21, before any label exists for the four decks below. The comparison set is now the six seasons 2011 to 2016. January to February 2016 was run on 2026-09-21 with the frozen model as season `tail2016`: 114 intervals, 59 tier 1, 338 tier 2.

| deck | season | seed | content | note |
|---|---|---|---|---|
| 2016 | tail2016 | 2016 | 18 tier 1 + 12 tier 2 | never seen by the model, same tail pass as 2015 |
| 2014census | tail2014 | 2014 | all 108 scope-gated episodes: 15 tier 1, 59 tier 2, 34 untiered | drawn from `beam_catalog_2014_census.csv`, a census so there is no sampling error and recall inside the gates is measurable |
| 2015held | tail2015 | 2015 | the 20 tier-1 episodes with no label + 12 unlabeled tier 2 | none of these were in the fit, labeled ids are excluded through `labeled_ids_2015.csv` |
| 2011b | tail2011 | 2011 | the 5 tier-1 episodes not in the earlier non-blind deck + 12 tier 2 | 2011 tier 1 is otherwise used up, the 15 non-blind labels stay a tuning set |

Decks 2012b and 2013 from section 1 complete the set. Their current labels are model-rater labels, so they count toward this comparison only once the user has labeled or spot-checked them.

H1, H2 and H3 apply unchanged, with the same thresholds and the same single scope cut. Two additions:

- **H4, recall inside the gates.** On the 2014 census, the share of labeled positives that are tier 1, and tier 1 or 2. Reference is the literature recall under the same model, 80% and 96% over 117 in-scope positives. The scorer prints this for every deck but it only means something for a census.
- **H5, held-out 2015.** Contamination among the 20 held-out tier-1 episodes against 1/14 in the original 2015 blind deck.

Label-free result already in hand for 2016 and not part of the blind test: of 37 in-scope literature positives, 33 match a tier-1 episode and 4 a tier-2 episode, none are missed.

## 6. Validation of the Te and beta refit

Written 2026-09-23, before the decks below exist. The model refit on 2026-09-23 adds log Te and log beta as inputs and is fit on every user label from 2011 to 2016, so the decks in sections 1 and 5 are now training data and cannot validate it. Held-out season AUC at fit time was 0.970 pooled, 0.905 for 2013. The cut stays at 0.82.

| deck | season | seed | direction | note |
|---|---|---|---|---|
| 2017 | tail2017 | 2017 | parallel | never labeled |
| 2023 | tail2023 | 2023 | mostly anti-parallel | never labeled |
| 2019v2 | tail2019 | 2019 | anti-parallel | tier 1 went from 33 to 96 under the refit, episodes in the 2026-09-20 deck excluded |

- **H6, contamination.** Tier-1 contamination at 0.82, per deck and pooled, by episode and by structure. The refit fails if any deck's point estimate exceeds 20%, or if the pooled estimate exceeds the 13% of the old model on the 2011 to 2016 decks.
- **H7, separation.** Pooled AUC of 0.96 or better.
- **H8, promoted episodes.** Tier-1 episodes that were tier 2 or untiered under the old model are scored as a separate group. If their contamination exceeds 20%, the Te and beta terms promote weak lines and need revisiting.

No threshold or scope is to be tuned on these decks.

## 7. Validation of the shape-feature model

Written 2026-09-24, after the two decks below were built and before any of their rows was labeled. The model refit on 2026-09-24 adds four inputs to the section 6 model: isolation (beam-flagged snapshots within 600 s of the episode on the same interval), duty cycle, line stability (share of members within 0.12 dex of the median beam energy) and beam-energy drift. It is fit on the same 2011 to 2016 user labels. The feature choice was made after the section 6 decks were scored, so those three decks are a design set, not a validation set.

Fit-time held-out-season numbers at the 0.82 cut: pooled AUC 0.974 (section 6 model 0.970), tier 1 holding 8 negatives of 166 and 158 of 193 positives (section 6 model 9 of 162 and 153). Interaction terms of prominence with Te and beta were tested and not adopted, they changed nothing.

| deck | season | seed | note |
|---|---|---|---|
| 2022 | tail2022 | 2022 | never had a deck, beams in both directions |
| 2024v2 | tail2024 | 2024 | episodes of the 2026-09-20 2024 deck excluded, no overlap |

- **H9, contamination.** Tier-1 contamination at 0.82, per deck and pooled, by episode and by structure. The model fails if either deck exceeds 20%.
- **H10, separation.** Pooled AUC of 0.96 or better.
- **H11, promotion.** Tier-1 draws that were not tier 1 under the section 6 model (catalog backups `beam_catalog_<year>_tebeta20260923.csv`) are scored separately, as in H8. The section 6 decks showed 4 of 21 there, this model should do better.

Who labels these decks is recorded with the result. No threshold, scope or feature is to be tuned on them.

