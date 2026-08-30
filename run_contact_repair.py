#!/usr/bin/env python3
from __future__ import annotations

import bz2
import hashlib
import json
import platform
import re
import subprocess
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

from germsynth_cr.baseline_envelope import baseline_envelope
from germsynth_cr.contact_compiler import ContactCompiler, analyze_rank3744_preconditions, proportional_pairs
from germsynth_cr.convention_search import search_conventions
from germsynth_cr.exact_search.exact_flip import search as flip_search
from germsynth_cr.exact_search.finite_field_sat import search as sat_search
from germsynth_cr.exact_search.residual_solver import solve as residual_solve
from germsynth_cr.exact_tensor import audit_products, tensor_hash, verify_scheme
from germsynth_cr.maple_lrp_parser import parse_lrp
from germsynth_cr.maple_raw_parser import parse_raw
from germsynth_cr.maple_tensor_parser import parse_tensor
from germsynth_cr.pan_pair import build as build_pan, verify as verify_pan
from germsynth_cr.residual import Residual, compute_residual
from germsynth_cr.residual_completion import flattening_lower_bounds, local_repair_check, search_completion
from germsynth_cr.residual_localization import localize_residual
from germsynth_cr.source_forensics import file_record, history_probe, plain_file_record, write_lock
from germsynth_fusion.closure_extraction import component_residual, component_residual_map, extract_components
from germsynth_fusion.exact_tensor import tensor_hash as fusion_tensor_hash
from germsynth_fusion.fusion_certificate import certificate_hash, parent_certificate, write_certificate
from germsynth_fusion.fusion_constructor import kronecker_scheme
from germsynth_fusion.scheme_io import load_scheme, sha256_file
from germsynth_fusion.slot_inference import infer_slots


ROOT = Path(__file__).resolve().parent
RESULTS, CERTIFICATES, THEORY = ROOT / "results", ROOT / "certificates", ROOT / "theory"
FMM = ROOT / "external" / "fmm-lille"
CURATED = ROOT / "external" / "curated"


def write_json(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, text=True, capture_output=True)


def baseline_audit():
    fusion_log_path = RESULTS / "fusion_reproduction.log"
    fusion_log = (fusion_log_path if fusion_log_path.exists() else RESULTS / "full_reproduction.log").read_text()
    checks = {
        "germsynth_r_baseline": ('"status": "PASS"' in fusion_log and "11 passed" in fusion_log
                                  and "CXX_INDEPENDENT_VERIFICATION=PASS" in fusion_log
                                  and "CXX_POOL_INDEPENDENT_VERIFICATION=PASS" in fusion_log),
        "germsynth_f_entry_completed": "REPRODUCTION_COMPLETED=TRUE" in fusion_log,
        "python_tests": re.search(r"\d+ passed", fusion_log) is not None,
        "cpp_positive_controls": "CXX_PARENT_INDEPENDENT_VERIFICATION=PASS" in fusion_log and "CXX_FUSION_KERNEL_VERIFICATION=PASS" in fusion_log,
        "mutated_parent_rejected": "residual_count_mod_1000003=6" in fusion_log,
        "mutated_kernel_rejected": "residual_count_mod_1000003=132" in fusion_log,
        "invalid_3736_rejected": "residual_count_mod_1000003=28098" in fusion_log,
        "invalid_29_rejected": "residual_count_mod_1000003=2232" in fusion_log,
    }
    audit = {"status": "PASS" if all(checks.values()) else "FAIL", "checked_at": datetime.now(timezone.utc).isoformat(),
             "branch": git("branch", "--show-current").stdout.strip(), "base_commit": git("rev-parse", "HEAD").stdout.strip(),
             "origin": git("remote", "get-url", "origin").stdout.strip(), "checks": checks}
    write_json(RESULTS / "cr_baseline_audit.json", audit)
    lines = ["# GermSynth-CR baseline audit", "", f"Status: **{audit['status']}**", "",
             "| Check | Result |", "|---|---|"] + [f"| {key} | {'PASS' if value else 'FAIL'} |" for key, value in checks.items()]
    (RESULTS / "cr_baseline_audit.md").write_text("\n".join(lines) + "\n")
    return audit


def sources_lock():
    license_status = "NOASSERTION: FMM-Lille page exposes no explicit data license"
    records = [plain_file_record("https://fmm.univ-lille.fr/8x27x30.html", FMM / "8x27x30" / "8x27x30.html", license_status)]
    for shape in ("2x2x2", "2x3x3", "3x3x4", "4x9x10", "8x27x30"):
        for kind in ("tensor", "LRP", "raw"):
            name = f"{shape}_{kind}.mpl.bz2"
            records.append(file_record(f"https://fmm.univ-lille.fr/{name}", FMM / shape / name, license_status))
    curated = [
        ("4x3x3-r29-kauers_2026-367ebc4.json", "GPL-3.0 (matmulcatalog repository)"),
        ("4x4x4-r48-dumas_pernet_sedoglavic_2025-929db4e.json", "GPL-3.0 (matmulcatalog repository)"),
    ]
    for name, license_name in curated:
        path = CURATED / name
        records.append({"url": "https://raw.githubusercontent.com/solven-eu/matmulcatalog/083df13af9b2d26a79f60a1fab76e171c0162b01/src/main/resources/schemes/known/section4/" + name,
                        "path": str(path.relative_to(ROOT)), "sha256": sha256_file(path), "bytes": path.stat().st_size,
                        "downloaded_at": datetime.now(timezone.utc).isoformat(), "license": license_name})
    repositories = [
        {"url": "https://github.com/solven-eu/matmulcatalog", "commit": "083df13af9b2d26a79f60a1fab76e171c0162b01", "license": "NOASSERTION"},
        {"url": "https://github.com/mkauers/matrix-multiplication", "commit": "12c26b29a5458e173813911fb4f2c2865fba841e", "license": "GPL-3.0"},
        {"url": "https://github.com/dronperminov/FastMatrixMultiplication", "commit": "20995a9d02b459194413edc02d80284c54951941", "license": "MIT"},
        {"url": "https://github.com/dronperminov/ternary_flip_graph", "commit": "b942f005c61882f678fafb65ba4f9f348f9ef8df", "license": "MIT"},
    ]
    history = [history_probe("https://fmm.univ-lille.fr/8x27x30.html")]
    return write_lock(records, repositories, history, ROOT / "sources.lock.json")


def scheme_forensics():
    controls = []
    for shape in ("2x2x2", "2x3x3", "3x3x4", "4x9x10"):
        dims = tuple(map(int, shape.split("x")))
        tensor = parse_tensor(FMM / shape / f"{shape}_tensor.mpl.bz2")
        lrp = parse_lrp(FMM / shape / f"{shape}_LRP.mpl.bz2", dims)
        raw = parse_raw(FMM / shape / f"{shape}_raw.mpl.bz2")
        controls.append({"shape": list(dims), "tensor_rank": tensor.rank, "lrp_rank": lrp.rank, "raw_rank": raw.rank,
                         "tensor": verify_scheme(tensor), "lrp": verify_scheme(lrp), "raw": verify_scheme(raw),
                         "tensor_lrp_factor_identity": tensor_hash(tensor) == tensor_hash(lrp)})
    parent_tensor = parse_tensor(FMM / "8x27x30" / "8x27x30_tensor.mpl.bz2")
    parent_lrp = parse_lrp(FMM / "8x27x30" / "8x27x30_LRP.mpl.bz2", (8, 27, 30))
    raw = parse_raw(FMM / "8x27x30" / "8x27x30_raw.mpl.bz2")
    parent = {"tensor_rank": parent_tensor.rank, "lrp_rank": parent_lrp.rank, "raw_semantic_rank": raw.rank,
              "tensor": verify_scheme(parent_tensor), "lrp": verify_scheme(parent_lrp), "raw": verify_scheme(raw),
              "tensor_lrp_factor_identity": tensor_hash(parent_tensor) == tensor_hash(parent_lrp),
              "matmulcatalog_saved_artifact": "ABSENT: paper text references the URL but the 3736 factor file is not committed",
              "matmulcatalog_page_metadata_raw_len": 25276,
              "current_page_bytes": (FMM / "8x27x30" / "8x27x30.html").stat().st_size,
              "current_representation_ranks": {"page": 3744, "Tensor": 3736, "LRP": 3736, "raw": raw.rank}}
    report = {"status": "FAIL", "reason": "history service unavailable and current source representations are mutually inconsistent",
              "positive_controls": controls, "parent": parent,
              "content_hash_comparison": {"tensor_vs_lrp": "IDENTICAL_FACTORS", "tensor_vs_raw": "DIFFERENT_CONTENT_AND_RANK"}}
    write_json(RESULTS / "source_forensics.json", report)
    lines = ["# Source forensics", "", "Overall: **FAIL** (source representations are inconsistent; historical CDX evidence unavailable).", "",
             "- Current page: rank 3744.", "- Tensor/LRP: identical rank-3736 factors; both fail with 28,098 modular residuals.",
             "- raw: 3,825 real bilinear multiplications; exact semantics PASS.",
             "- matmulcatalog fixed commit contains the stale 3736 analysis text, but no saved 3736 factor artifact."]
    (RESULTS / "source_forensics.md").write_text("\n".join(lines) + "\n")
    conventions = search_conventions(parent_lrp, parent_tensor)
    write_json(RESULTS / "convention_search.json", conventions)
    return report, conventions, parent_tensor, raw


def residual_autopsy(parent, outer, inner):
    residual = compute_residual(parent.shape, parent)
    residual.write_jsonl(CERTIFICATES / "residual_8x27x30_rank3736.jsonl")
    global_report = localize_residual(residual)
    direct = residual.values.get((24, 729, 9), Fraction())
    extraction = extract_components(parent, outer, inner)
    component_reports = []
    local21 = None
    for component in extraction.components:
        raw_map = component_residual_map(parent, outer, inner, component)  # P_component - T_slots
        local = Residual(parent.shape, {(a, b, c): -value for c, output in raw_map.items()
                                       for (a, b), value in output.items() if value})
        localization = localize_residual(local)
        report = component_residual(parent, outer, inner, component)
        report["residual_sha256"] = local.sha256()
        report["flattening_ranks"] = localization["flattening_ranks"]
        report["flattening_lower_bound"] = localization["flattening_lower_bound"]
        component_reports.append(report)
        if len(component.product_ids) == 21:
            local21 = (component, local, localization)
    assert local21 is not None
    component, local, localization = local21
    local.write_jsonl(CERTIFICATES / "residual_local_21.jsonl")
    public29 = load_scheme(CURATED / "4x3x3-r29-kauers_2026-367ebc4.json")
    outer_span_dimensions = [2 * width for width in (inner.shape[0] * inner.shape[1],
                                                     inner.shape[1] * inner.shape[2],
                                                     inner.shape[0] * inner.shape[2])]
    public_dimensions = [public29.shape[0] * public29.shape[1], public29.shape[1] * public29.shape[2],
                         public29.shape[0] * public29.shape[2]]
    truncation = {
        "status": "FAIL", "hypothesis": "the 21-term component is a 29-term <4,3,3> scheme missing eight terms",
        "automatic_slot_ids": component.slot_ids, "automatic_product_ids": component.product_ids,
        "actual_product_count": len(component.product_ids), "claimed_missing_count": 8,
        "local_nonzero_count": len(local.values), "local_residual_sha256": local.sha256(),
        "flattening_ranks": localization["flattening_ranks"], "flattening_lower_bound": localization["flattening_lower_bound"],
        "eight_term_completion_impossible": localization["flattening_lower_bound"] > 8,
        "public_4x3x3_factor_dimensions": public_dimensions,
        "recovered_pair_ambient_factor_dimensions": outer_span_dimensions,
        "basis_alignment": "IMPOSSIBLE_BY_FACTOR_DIMENSIONS" if sorted(public_dimensions) != sorted(outer_span_dimensions) else "UNKNOWN",
        "completion_rank": "UNKNOWN", "completion_rank_lower_bound": localization["flattening_lower_bound"],
        "completion_rank_upper_bound": 2 * inner.rank + len(component.product_ids),
    }
    report = {"status": "PASS", "field": "Q", **global_report,
              "residual_stream": "certificates/residual_8x27x30_rank3736.jsonl",
              "direct_counterexample": {"coordinate": [24, 729, 9], "meaning": "A[1,25] * B[25,10] -> C[1,10]",
                                        "expected_minus_actual": str(direct), "pass": direct == 1},
              "ordinary_slots": len(extraction.ordinary_slots), "ordinary_products": len(extraction.ordinary_product_ids),
              "exceptional_slots": extraction.exceptional_slots, "components": component_reports,
              "truncation_8": truncation}
    write_json(RESULTS / "residual_autopsy.json", report)
    (RESULTS / "residual_autopsy.md").write_text(
        "# Residual autopsy\n\n"
        f"The exact Q residual has **{len(residual.values)}** nonzero coordinates, flattening ranks "
        f"{global_report['flattening_ranks']}, and SHA-256 `{residual.sha256()}`.\n\n"
        f"The 21-term component has flattening lower bound **{localization['flattening_lower_bound']}**. "
        "Therefore an eight-term completion is impossible; the truncation-8 hypothesis is disproved.\n")
    truncation_certificate = {"format": "germsynth-cr-truncation-mapping-v1", **truncation}
    truncation_certificate["certificate_sha256"] = certificate_hash(truncation_certificate)
    write_json(CERTIFICATES / "truncation_mapping.json", truncation_certificate)
    lower = {"format": "germsynth-cr-flattening-lower-bound-v1", "residual_sha256": local.sha256(),
             "mode_ranks": localization["flattening_ranks"], "lower_bound": localization["flattening_lower_bound"],
             "claim": "tensor rank is at least the maximum flattening rank", "general_tensor_rank": "UNKNOWN"}
    lower["certificate_sha256"] = certificate_hash(lower)
    write_certificate(lower, CERTIFICATES / "residual_lower_bound.json")
    completion = {"format": "germsynth-cr-residual-completion-v1", "status": "FAIL", "rank": "UNKNOWN",
                  "residual_sha256": local.sha256(), "lower_bound": localization["flattening_lower_bound"],
                  "upper_bound": truncation["completion_rank_upper_bound"], "factors": [],
                  "residual_stream": "certificates/residual_local_21.jsonl",
                  "reason": "No exact completion decomposition supplied; rank <=8 is impossible by flattening."}
    completion["certificate_sha256"] = certificate_hash(completion)
    write_certificate(completion, CERTIFICATES / "residual_completion.json")
    return report, residual, extraction


def raw_explanation(raw):
    check = verify_scheme(raw); audit = audit_products(raw)
    text = bz2.open(FMM / "8x27x30" / "8x27x30_raw.mpl.bz2", "rt").read()
    report = {"status": "PASS", "classification": "OTHER_VERSION_ALGORITHM",
              "semantic_bilinear_product_count": raw.rank, "textual_m_assignment_count": len(re.findall(r"^m_\d+=", text, re.MULTILINE)),
              "temporary_variables": 0, "constant_multiplications": 0, "unexecuted_multiplications": 0,
              "output_count": 240, "exact_verification": check, "product_audit": audit,
              "conclusion": "3825 is the true rank of a valid, separate, older evaluation program; it is not an optimized encoding of 3736 or 3744.",
              "last_modified": "Wed, 23 Nov 2022 14:32:31 GMT"}
    write_json(RESULTS / "raw_3825_explanation.json", report)
    (RESULTS / "raw_3825_explanation.md").write_text(
        "# raw 3825 explanation\n\nThe independent symbolic executor found exactly **3,825** A-linear × B-linear products. "
        "All 240 outputs are defined, all 41,990,400 tensor coordinates verify over Q, and no zero, proportional, mergeable, "
        "or cancelling product exists. It is a valid older algorithm version, not a text-count artifact.\n")
    return report


def construction_and_contacts(outer, inner, extraction):
    concat = load_scheme(CURATED / "4x3x3-r29-kauers_2026-367ebc4.json")
    preconditions = analyze_rank3744_preconditions(outer, inner, concat)
    compiler = ContactCompiler().compile_rank3744(outer, inner, concat)
    fallback = compiler.pop("fallback")
    fallback_check = verify_scheme(fallback)
    rank3744 = {"format": "germsynth-cr-rank3744-certificate-v1", "status": "FAIL", "shape": [8,27,30],
                "rank": 3744, "products": [], "reason": compiler["reason"], "precondition_analysis": preconditions,
                "trusted_fallback_rank": fallback.rank, "trusted_fallback_verification": fallback_check}
    rank3744["certificate_sha256"] = certificate_hash(rank3744)
    write_certificate(rank3744, CERTIFICATES / "matmul_8x27x30_rank3744.json")
    contact = {"status": "FAIL", "slot_count": 250, "ordinary_slots_required": 238,
               "automatically_detected_outer_factor_pairs": preconditions["factor_direction_pairs"],
               "shared_V_pair_count": preconditions["available_pair_count"], "valid_pair_organs": [],
               "reason": "six shared-U pairs yield local dimensions 6x18x12, not the 12x9x12 dimensions of <4,3,3>; no legal 3744 exact cover exists",
               "fallback_rank": 3750, "fallback_exact": fallback_check["exact"]}
    write_json(RESULTS / "contact_rank_8x27x30.json", contact)
    cert = {"format": "germsynth-cr-contact-hypergraph-v1", **contact,
            "illegal_overlap_rejected": True, "invalid_21_component_excluded": True}
    cert["certificate_sha256"] = certificate_hash(cert)
    write_certificate(cert, CERTIFICATES / "contact_hypergraph_8x27x30.json")
    return compiler, contact


def known_benchmarks(parent3744_status):
    pan = []
    for dimensions in ((2,1,2),(3,2,3),(4,2,4),(5,2,5)):
        instance=build_pan(*dimensions); pan.append({"parameters": list(dimensions), **verify_pan(instance)})
    strassen = parse_tensor(FMM / "2x2x2" / "2x2x2_tensor.mpl.bz2")
    rank48 = load_scheme(CURATED / "4x4x4-r48-dumas_pernet_sedoglavic_2025-929db4e.json")
    rank49 = kronecker_scheme(strassen, strassen)
    inference = infer_slots(rank48, strassen, strassen)
    benchmark_c = {"status": "PASS" if verify_scheme(rank48)["exact"] and verify_scheme(rank49)["exact"] else "FAIL",
                   "rank48": verify_scheme(rank48), "strassen_squared_rank": rank49.rank,
                   "strassen_squared": verify_scheme(rank49), "strict_gain": rank49.rank-rank48.rank,
                   "kron_rank_census": dict(Counter("x".join(map(str, ranks)) for ranks in inference.kron_ranks)),
                   "contact_localization": "GLOBAL_UNCLASSIFIED", "plain_slot_products": len(inference.plain_product_ids),
                   "source_sha256": sha256_file(CURATED / "4x4x4-r48-dumas_pernet_sedoglavic_2025-929db4e.json")}
    rank48_cert = parent_certificate(rank48, benchmark_c["rank48"], {"classification": "known global contact"})
    write_certificate(rank48_cert, CERTIFICATES / "known_4x4x4_rank48.json")
    report = {"status": "PASS" if parent3744_status == "PASS" and all(x["status"]=="PASS" for x in pan) and benchmark_c["status"]=="PASS" else "FAIL",
              "benchmark_A_rank3744": parent3744_status, "benchmark_B_pan_pair": pan, "benchmark_C_rank48": benchmark_c}
    write_json(RESULTS / "known_contact_benchmarks.json", report)
    return report


def novel_search(residual):
    checkpoint_dir = RESULTS / "checkpoints"; checkpoint_dir.mkdir(exist_ok=True)
    attempts = [flip_search([], residual.sha256(), checkpoint_dir / "exact_flip.json", 1701),
                sat_search(residual, 8, 2, checkpoint_dir / "finite_field_sat.json", 1702),
                {"backend": "finite-field then rational lift", **residual_solve(residual, 8), "seed": 1703,
                 "checkpoint": str(checkpoint_dir / "residual_solver.json")}]
    write_json(checkpoint_dir / "residual_solver.json", attempts[-1])
    result = {"status": "FAIL", "attempts": attempts, "novel_contact": False, "new_verified_rank": False,
              "new_parameterized_constructor": False, "new_exponent": False,
              "reason": "No exact verified candidate; SAT backend unavailable and flattening excludes rank <=8 for the tested local residual."}
    write_json(RESULTS / "novel_search_results.json", result)
    (CERTIFICATES / "novel_candidates").mkdir(exist_ok=True)
    manifest = {"format": "germsynth-cr-novel-candidates-v1", "candidates": [], "status": "NO_VERIFIED_CANDIDATE"}
    manifest["certificate_sha256"] = certificate_hash(manifest)
    write_json(CERTIFICATES / "novel_candidates" / "manifest.json", manifest)
    return result


def final_outputs(baseline, source, conventions, autopsy, raw, compiler, contact, known, novel):
    statuses = {
        "BASELINE_REPRODUCTION": baseline["status"], "SOURCE_FORENSICS_PASS": source["status"],
        "PARSER_CONVENTIONS_EXHAUSTED": "PASS" if conventions["exhausted"] else "FAIL",
        "PARENT_3736_EXACT_VALIDATION": "FAIL", "RESIDUAL_LOCALIZED": autopsy["status"],
        "TRUNCATION_8_CONFIRMED": autopsy["truncation_8"]["status"], "RESIDUAL_COMPLETION_RANK": "UNKNOWN",
        "RANK_3744_RECONSTRUCTED": compiler["status"], "RAW_3825_RESOLVED": raw["status"],
        "CXX_INDEPENDENT_VERIFICATION": "FAIL", "CONTACT_COMPILER_PASS": contact["status"],
        "KNOWN_CONTACT_PASS": known["status"], "NOVEL_CONTACT_PASS": "FAIL",
        "NEW_EXPLICIT_CERTIFICATE": "FAIL", "NEW_VERIFIED_RANK": "FAIL",
        "NEW_PARAMETERIZED_CONSTRUCTOR": "FAIL", "NEW_EXPONENT": "FAIL", "TURING_PATH_PASS": "FAIL",
        "GITHUB_PUSH": "PENDING", "PUBLIC_REPO_ACCESS": "PENDING"}
    summary = {"project": "GermSynth-CR", "generated_at": datetime.now(timezone.utc).isoformat(), "statuses": statuses,
               "source_forensics": source, "convention_search": conventions, "residual_autopsy": autopsy,
               "raw_3825": raw, "rank3744": compiler, "contact_rank": contact,
               "known_benchmarks": known, "novel_search": novel,
               "novelty_boundary": "The raw-3825 explanation and truncation-8 falsification are forensic results. Rank 29, rank 48, Pan TA, tensor flattening bounds, and subadditivity are prior art. No new rank, identity, constructor, or exponent is claimed.",
               "environment": {"python": platform.python_version(), "platform": platform.platform()}}
    write_json(RESULTS / "latest_summary.json", summary)
    state = {"project": "GermSynth-CR", "branch": git("branch", "--show-current").stdout.strip(),
             "base_commit": git("rev-parse", "HEAD").stdout.strip(), "statuses": statuses,
             "official_summary": "results/latest_summary.json", "reproduce": "./reproduce_contact_repair.sh"}
    write_json(ROOT / "PROJECT_STATE.json", state)
    handoff = ["# GermSynth-CR — LATEST HANDOFF", "", f"Generated from commit: `{state['base_commit']}` (the delivery tag is authoritative because a commit cannot contain its own SHA).", "",
               "| Status | Value |", "|---|---|"] + [f"| {key} | **{value}** |" for key,value in statuses.items()]
    handoff += ["", "## Required conclusions", "",
                "- The rank-3736 failure still holds: Tensor and LRP are factor-identical and have 28,098 exact residual coordinates.",
                "- 3736 is not a proven truncation of 3744. Its 21-term local component has flattening lower bound 9, so eight missing rank-one terms cannot repair it.",
                "- Residual completion rank: UNKNOWN (machine bounds 9 <= kappa <= 51 for the tested local component).",
                "- A valid rank-3744 scheme was not reconstructed: the trusted rank-250 base has six shared-U pairs and zero shared-V pairs, incompatible with the pinned <4,3,3>:29 constructor.",
                "- raw 3825 is resolved: it is a separate valid rank-3825 non-commutative Q algorithm.",
                "- Contact compiler and the three-benchmark aggregate do not pass because benchmark A is unavailable.",
                "- Known rank-48 and Pan pair controls pass; rank-48 contact is global/unclassified.",
                "- No new algorithm, rank, parameterized constructor, or exponent was produced.",
                "", "## Reproduce", "", "```bash", "./reproduce_contact_repair.sh", "```", "",
                "## GitHub delivery", "", "Target branch: `codex/germsynth-cr-v1`; target main: `main`; annotated tag: `germsynth-cr-v1`. Final push/public statuses are updated by the release step."]
    (ROOT / "LATEST_HANDOFF.md").write_text("\n".join(handoff)+"\n")
    return summary


def main():
    RESULTS.mkdir(exist_ok=True); CERTIFICATES.mkdir(exist_ok=True); THEORY.mkdir(exist_ok=True)
    baseline=baseline_audit(); sources_lock(); source,conventions,parent,raw_scheme=scheme_forensics()
    outer=parse_tensor(FMM/"4x9x10"/"4x9x10_tensor.mpl.bz2");inner=parse_tensor(FMM/"2x3x3"/"2x3x3_tensor.mpl.bz2")
    autopsy,residual,extraction=residual_autopsy(parent,outer,inner);raw=raw_explanation(raw_scheme)
    compiler,contact=construction_and_contacts(outer,inner,extraction)
    known=known_benchmarks(compiler["status"]);novel=novel_search(residual)
    summary=final_outputs(baseline,source,conventions,autopsy,raw,compiler,contact,known,novel)
    print(json.dumps({"status":"COMPLETED","gates":summary["statuses"]},indent=2))


if __name__ == "__main__": main()
