# GermSynth-CR — LATEST HANDOFF

Release-state parent commit: `9e3ebe640f373ba91cd60eaef4369daf5cd179f8`. The checked-out `main` SHA is authoritative because a commit cannot contain its own SHA.

| Status | Value |
|---|---|
| BASELINE_REPRODUCTION | **PASS** |
| SOURCE_FORENSICS_PASS | **FAIL** |
| PARSER_CONVENTIONS_EXHAUSTED | **PASS** |
| PARENT_3736_EXACT_VALIDATION | **FAIL** |
| RESIDUAL_LOCALIZED | **PASS** |
| TRUNCATION_8_CONFIRMED | **FAIL** |
| RESIDUAL_COMPLETION_RANK | **UNKNOWN** |
| RANK_3744_RECONSTRUCTED | **FAIL** |
| RAW_3825_RESOLVED | **PASS** |
| CXX_INDEPENDENT_VERIFICATION | **FAIL** |
| CONTACT_COMPILER_PASS | **FAIL** |
| KNOWN_CONTACT_PASS | **FAIL** |
| NOVEL_CONTACT_PASS | **FAIL** |
| NEW_EXPLICIT_CERTIFICATE | **FAIL** |
| NEW_VERIFIED_RANK | **FAIL** |
| NEW_PARAMETERIZED_CONSTRUCTOR | **FAIL** |
| NEW_EXPONENT | **FAIL** |
| TURING_PATH_PASS | **FAIL** |
| GITHUB_PUSH | **FAIL** |
| PUBLIC_REPO_ACCESS | **FAIL** |

## Required conclusions

- The rank-3736 failure still holds: Tensor and LRP are factor-identical and have 28,098 exact residual coordinates.
- 3736 is not a proven truncation of 3744. Its 21-term local component has flattening lower bound 9, so eight missing rank-one terms cannot repair it.
- Residual completion rank: UNKNOWN (machine bounds 9 <= kappa <= 51 for the tested local component).
- A valid rank-3744 scheme was not reconstructed: the trusted rank-250 base has six shared-U pairs and zero shared-V pairs, incompatible with the pinned <4,3,3>:29 constructor.
- raw 3825 is resolved: it is a separate valid rank-3825 non-commutative Q algorithm.
- Contact compiler and the three-benchmark aggregate do not pass because benchmark A is unavailable.
- Known rank-48 and Pan pair controls pass; rank-48 contact is global/unclassified.
- No new algorithm, rank, parameterized constructor, or exponent was produced.

## Reproduce

```bash
./reproduce_contact_repair.sh
```

## GitHub delivery

Release-state parent: `9e3ebe640f373ba91cd60eaef4369daf5cd179f8`. Work branch and main contain this corrective state after synchronization. Annotated tag `germsynth-cr-v1` points to the preceding release commit `9e3ebe640f373ba91cd60eaef4369daf5cd179f8` and is not force-moved. SSH branch/main/tag pushes passed, but the combined GITHUB gate is FAIL because anonymous HTTPS access fails; repository visibility was not changed.
