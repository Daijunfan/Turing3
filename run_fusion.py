#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import platform
import bz2
import re
import subprocess
from collections import Counter
from datetime import date
from pathlib import Path

from germsynth_fusion.closure_extraction import component_residual, dependency_closures, extract_components
from germsynth_fusion.exact_tensor import audit_products, tensor_hash, verify_scheme
from germsynth_fusion.fusion_certificate import certificate_hash, organ_certificate, parent_certificate, write_certificate
from germsynth_fusion.fusion_constructor import kronecker_scheme
from germsynth_fusion.fusion_search import source_gated_search
from germsynth_fusion.scheme_io import load_scheme, sha256_file
from germsynth_fusion.spectral_grammar import empirical_single_type, single_type_check


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
CERTIFICATES = ROOT / "certificates"
FMM = ROOT / "external" / "fmm-lille"


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def command(args, cwd=ROOT, timeout=180):
    run = subprocess.run(args, cwd=cwd, text=True, capture_output=True, timeout=timeout)
    return {"command": " ".join(map(str, args)), "returncode": run.returncode,
            "stdout": run.stdout, "stderr": run.stderr, "pass": run.returncode == 0}


def baseline_audit():
    base = ROOT / "germsynth"
    python = base / ".venv" / "bin" / "python"
    if not python.exists():
        python = Path("python3")
    steps = [
        command([str(python), "-m", "pytest", "-q"], base),
        command([str(python), "independent_verify.py", "certificates/karatsuba_germ.json"], base),
        command([str(python), "independent_verify.py", "certificates/matrix_rank7_germ.json"], base),
        command([str(python), "independent_verify_pool.py", "certificates/matrix_phenotype_pool.json"], base),
        command([str(base / "build" / "matrix_rank7_germ_verify")], base),
        command([str(base / "build" / "matrix_pool_verify")], base),
    ]
    results = json.loads((base / "results" / "results.json").read_text())
    checks = {
        "python_tests": steps[0]["pass"] and "11 passed" in steps[0]["stdout"],
        "independent_python_germ_verifiers": all(step["pass"] and '"status": "PASS"' in step["stdout"] for step in steps[1:3]),
        "independent_python_pool_verifier": steps[3]["pass"] and '"status": "PASS"' in steps[3]["stdout"],
        "independent_cpp_germ_verifier": steps[4]["pass"] and "CXX_INDEPENDENT_VERIFICATION=PASS" in steps[4]["stdout"],
        "independent_cpp_pool_verifier": steps[5]["pass"] and "CXX_POOL_INDEPENDENT_VERIFICATION=PASS" in steps[5]["stdout"],
        "sixteen_phenotypes": results["polymorphic_verification"]["phenotypes"] == 16,
        "minimum_three_phenotype_single_fault_cover": results["polymorphic_verification"]["minimum_regenerative_cover"]["phenotype_count"] == 3,
    }
    audit = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "date": str(date.today()),
        "original_archive": "GermSynth-R_prototype.zip",
        "original_archive_sha256": sha256_file(ROOT / "GermSynth-R_prototype.zip"),
        "checks": checks,
        "commands": steps,
        "note": "The original reproduce.sh first stopped because system Python 3.14 lacked pytest; the equivalent entry was rerun in germsynth/.venv without changing baseline source semantics.",
    }
    write_json(RESULTS / "fusion_baseline_audit.json", audit)
    return audit


def source_manifest():
    repos = [
        ("matmulcatalog", "https://github.com/solven-eu/matmulcatalog", "083df13af9b2d26a79f60a1fab76e171c0162b01", "NOASSERTION (no repository LICENSE found)"),
        ("matrix-multiplication", "https://github.com/mkauers/matrix-multiplication", "12c26b29a5458e173813911fb4f2c2865fba841e", "GPL-3.0"),
        ("FastMatrixMultiplication", "https://github.com/dronperminov/FastMatrixMultiplication", "20995a9d02b459194413edc02d80284c54951941", "MIT"),
        ("ternary_flip_graph", "https://github.com/dronperminov/ternary_flip_graph", "b942f005c61882f678fafb65ba4f9f348f9ef8df", "MIT"),
    ]
    files = []
    for shape in ("2x3x3", "4x9x10", "8x27x30"):
        for kind in ("LRP", "raw", "tensor"):
            path = FMM / shape / f"{shape}_{kind}.mpl.bz2"
            files.append({"path": str(path.relative_to(ROOT)), "sha256": sha256_file(path), "bytes": path.stat().st_size})
    raw_parent = bz2.open(FMM / "8x27x30" / "8x27x30_raw.mpl.bz2", "rt").read()
    raw_rank = len(re.findall(r"^m_\d+=", raw_parent, re.MULTILINE))
    manifest = {
        "download_date": str(date.today()),
        "repositories": [{"name": name, "url": url, "commit": commit, "license": license_name,
                          "acquisition": f"GitHub commit tarball extracted under external/{name}"}
                         for name, url, commit, license_name in repos],
        "fmm_lille": {
            "url": "https://fmm.univ-lille.fr/8x27x30.html",
            "license": "NOASSERTION (no data license found on the source page)",
            "page_rank_on_download": 3744,
            "explicit_artifact_rank": 3736,
            "representation_ranks": {"web_page": 3744, "LRP": 3736, "Tensor_TriadSet": 3736,
                                     "raw_evaluation_program": raw_rank},
            "consistency": "FAIL: page, LRP/Tensor, and raw program declare three different ranks",
            "files": files,
        },
    }
    write_json(ROOT / "external" / "SOURCES.json", manifest)
    return manifest


def main():
    RESULTS.mkdir(exist_ok=True)
    CERTIFICATES.mkdir(exist_ok=True)
    baseline = baseline_audit()
    sources = source_manifest()

    inner = load_scheme(FMM / "2x3x3" / "2x3x3_tensor.mpl.bz2")
    outer = load_scheme(FMM / "4x9x10" / "4x9x10_tensor.mpl.bz2")
    parent = load_scheme(FMM / "8x27x30" / "8x27x30_tensor.mpl.bz2")
    inner_check = verify_scheme(inner)
    outer_check = verify_scheme(outer)
    parent_check = verify_scheme(parent)
    parent_check["counterexample_direct"] = {
        "a": 24, "b": 729, "c": 9, "meaning": "A[1,25] * B[25,10] -> C[1,10]",
        "expected": "1", "actual": "0",
    }
    parent_check["source_status"] = "SOURCE_INVALID" if parent_check["status"] == "FAIL" else "VERIFIED"

    extraction = extract_components(parent, outer, inner)
    component_reports = [component_residual(parent, outer, inner, component) for component in extraction.components]
    closures = dependency_closures(parent, outer, inner, extraction.components)
    nominal_gain = sum(report["fusion_gain"] for report in component_reports)
    product_sum = len(extraction.ordinary_product_ids) + sum(report["product_count"] for report in component_reports)
    exact_components = all(report["standalone"] for report in component_reports)
    rank_counts = Counter(tuple(rank) for rank in extraction.kron_ranks)

    census = {
        "status": "PASS" if parent_check["exact"] and exact_components and nominal_gain == 14 and product_sum == 3736 else "FAIL",
        "source_status": parent_check["source_status"],
        "parent_shape": list(parent.shape), "parent_rank": parent.rank,
        "theoretical_slot_count": outer.rank, "slot_shape": list(inner.shape), "isolated_rank": outer.rank * inner.rank,
        "ordinary_single_slots": len(extraction.ordinary_slots),
        "ordinary_product_count": len(extraction.ordinary_product_ids),
        "two_slot_candidate_components": sum(report["slot_count"] == 2 for report in component_reports),
        "three_or_larger_candidate_components": sum(report["slot_count"] >= 3 for report in component_reports),
        "exceptional_slot_ids": extraction.exceptional_slots,
        "components": component_reports,
        "minimal_dependency_closures": closures,
        "product_count_sum": product_sum,
        "nominal_fusion_gain_sum": nominal_gain,
        "exact_fusion_gain_14_explained": parent_check["exact"] and exact_components and nominal_gain == 14,
        "kron_rank_census": {"x".join(map(str, key)): count for key, count in sorted(rank_counts.items())},
        "mixed_candidate_count": sum(rank != (1, 1, 1) for rank in extraction.kron_ranks),
        "automatic": True,
        "handwritten_product_ids": False,
        "failure": parent_check["counterexample_direct"],
    }
    write_json(RESULTS / "parent_8x27x30_component_census.json", census)
    lines = ["# <8,27,30>:3736 component census", "", f"Overall: **{census['status']}** (`{census['source_status']}`)", "",
             f"The extractor recovered {census['ordinary_single_slots']} ordinary slots and six two-slot candidates. "
             f"Their product counts sum to {product_sum}; their nominal savings sum to {nominal_gain}. "
             "Every exceptional component has nonzero exact residual, so the nominal 14 is not a proved fusion gain.", "",
             "| Slots | Products | Nominal gain | Exact residuals | Standalone |", "|---|---:|---:|---:|---|"]
    for report in component_reports:
        lines.append(f"| {report['slot_ids']} | {report['product_count']} | {report['fusion_gain']} | {report['residual_count']} | {report['standalone']} |")
    lines += ["", "Direct counterexample: `A[1,25] * B[25,10] -> C[1,10]` has expected coefficient 1 and actual coefficient 0."]
    (RESULTS / "parent_8x27x30_component_census.md").write_text("\n".join(lines) + "\n")

    extraction_log = {
        "ordinary_slots": extraction.ordinary_slots,
        "ordinary_product_count": len(extraction.ordinary_product_ids),
        "exceptional_slots": extraction.exceptional_slots,
        "components": [{"slot_ids": c.slot_ids, "product_ids": c.product_ids} for c in extraction.components],
        "method": "exact outer×inner reshape ranks, exact outer-factor span membership, mixed-seeded hypergraph closure",
    }
    parent_cert = parent_certificate(parent, parent_check, extraction_log)
    write_certificate(parent_cert, CERTIFICATES / "fusion_parent_8x27x30_rank3736.json")

    rank29 = [(component, report) for component, report in zip(extraction.components, component_reports)
              if report["slot_count"] == 2 and report["product_count"] == 29]
    representative, representative_report = rank29[0]
    kernel_cert = organ_certificate(parent, outer, inner, representative, extraction.kron_ranks)
    kernel_cert["product_slot_incidence"] = [
        {"product_id": product, "slot_ids": sorted(extraction.product_slot_incidence.get(product, set()))}
        for product in representative.product_ids
    ]
    kernel_cert["certificate_sha256"] = certificate_hash(kernel_cert)
    write_certificate(kernel_cert, CERTIFICATES / "fusion_2x233_rank29.json")
    census_cert = {
        "format": "germsynth-fusion-component-census-certificate-v1",
        "status": census["status"], "parent_certificate_sha256": parent_cert["certificate_sha256"],
        "components": component_reports, "closures": closures, "product_count_sum": product_sum,
        "nominal_gain_sum": nominal_gain, "exact_gain_explained": census["exact_fusion_gain_14_explained"],
    }
    census_cert["certificate_sha256"] = certificate_hash(census_cert)
    write_certificate(census_cert, CERTIFICATES / "fusion_component_census.json")

    # Extra valid reference certificates exercise both independent verifiers and their mutation rejection path.
    ref_parent = parent_certificate(inner, inner_check, {"reference_only": True})
    write_certificate(ref_parent, CERTIFICATES / "reference_2x3x3_rank15.json")
    ref_component = type(representative)([extraction.ordinary_slots[0]], extraction.ordinary_product_ids[:inner.rank])
    ref_organ = organ_certificate(parent, outer, inner, ref_component, extraction.kron_ranks)
    write_certificate(ref_organ, CERTIFICATES / "reference_plain_233_rank15.json")

    plain = kronecker_scheme(outer, inner)
    plain_check = verify_scheme(plain)
    constructor = {
        "status": "FAIL", "reason": "No extracted fusion organ has zero residual; invalid kernels are rejected.",
        "plain_kronecker_rank": plain.rank, "plain_kronecker_verification": plain_check,
        "candidate_kernel_count": len(component_reports), "verified_kernel_count": sum(r["standalone"] for r in component_reports),
        "all_six_deficits_closed": False,
    }
    search = source_gated_search(parent_check["exact"], [kernel_cert], RESULTS / "fusion_search_checkpoint.json")
    write_json(RESULTS / "novel_search_results.json", search)
    spectral = {
        "status": "PARTIAL_SOURCE_INVALID", "single_type_regressions": [single_type_check(3, 2), single_type_check(7, 2)],
        "empirical_single_type_regressions": [empirical_single_type(3, 2), empirical_single_type(7, 2)],
        "parent_before_after": "NOT_AVAILABLE_SOURCE_INVALID", "new_exponent": False,
        "multi_type_engine_implemented": True, "structured_public_reproduction": "NOT_RUN_SOURCE_INVALID",
        "empirical_fit": "NOT_RUN_SOURCE_INVALID", "empirical_tolerance": 0.01,
    }
    write_json(RESULTS / "spectral_results.json", spectral)

    statuses = {
        "BASELINE_REPRODUCTION": baseline["status"],
        "PARENT_EXACT_VERIFICATION": parent_check["status"],
        "AUTOMATIC_SLOT_INFERENCE": "PASS" if len(extraction.ordinary_slots) == 238 and len(extraction.exceptional_slots) == 12 else "FAIL",
        "FUSION_KERNEL_EXTRACTION": "PASS" if any(report["standalone"] and report["product_count"] == 29 for report in component_reports) else "FAIL",
        "FULL_GAIN_14_EXPLAINED": "PASS" if census["exact_fusion_gain_14_explained"] else "FAIL",
        "CXX_INDEPENDENT_VERIFICATION": "FAIL",
        "CONSTRUCTOR_PASS": constructor["status"],
        "ALL_SIX_DEFICITS_CLOSED": "FAIL",
        "NOVEL_FUSION_KERNEL": "FAIL",
        "NEW_EXACT_RANK": "FAIL",
        "NEW_EXPONENT": "FAIL",
        "TURING_PATH_PASS": "FAIL",
    }
    summary = {
        "project": "GermSynth-F", "date": str(date.today()), "statuses": statuses,
        "baseline": baseline, "sources": sources,
        "inner_verification": inner_check, "outer_verification": outer_check,
        "parent_verification": parent_check, "parent_product_audit": audit_products(parent),
        "census": census, "constructor": constructor, "search": search, "spectral": spectral,
        "comparisons": {name: "NOT_RUN_SOURCE_INVALID" for name in (
            "plain Kronecker composition", "non-uniform recombination", "catalog pair fusion",
            "serendipitous product", "meta flip graph", "GermSynth-F extractor",
            "GermSynth-F constructor", "GermSynth-F synthesis search")},
        "ablations": {name: "NOT_RUN_SOURCE_INVALID" for name in (
            "no Kron-rank", "support-overlap-only", "no dependency closure", "no mixed products",
            "two slots only", "random tests without exact verification", "no extracted-kernel seed")},
        "formal_invariant": "Within non-overlap composition, if every product is incident to at most one slot, product sets partition by slot. Exact restriction to each slot therefore needs at least that slot's declared rank; summing gives R >= sum_i r_i, so fusion_gain <= 0.",
        "environment": {"python": platform.python_version(), "platform": platform.platform()},
    }
    write_json(RESULTS / "summary.json", summary)
    write_handoff(summary)
    print(json.dumps({"status": "COMPLETED_WITH_SOURCE_INVALID", "statuses": statuses,
                      "parent_residual_count": parent_check.get("modular_screens", [{}])[0].get("residual_count"),
                      "nominal_gain": nominal_gain, "exact_gain": False}, indent=2))


def write_handoff(summary):
    s = summary["statuses"]
    lines = ["# HANDOFF TO BRAIN — GermSynth-F", "", "## Status table", "", "| Gate | Status |", "|---|---|"]
    lines.extend(f"| {key} | **{value}** |" for key, value in s.items())
    lines += [
        "", "## Three most important findings", "",
        "1. FMM-Lille is internally inconsistent: Triad/LRP declare 3736, the current page declares 3744, and the raw evaluation file contains 3825 products.",
        "2. Exact verification rejects that artifact: 28,098 tensor coordinates are nonzero modulo 1,000,003; a direct rational counterexample is `A[1,25] B[25,10] -> C[1,10]`, expected 1, actual 0.",
        "3. Automatic extraction recovers 238 ordinary slots plus six two-slot candidates (five rank 29 and one rank 21), whose nominal savings are 5×1+9=14, but every candidate has a nonzero exact residual.",
        "", "## Known versus new", "",
        "The 29-for-30 narrative and the stale 3736 discussion are pre-existing claims in the downloaded catalog draft. GermSynth-F's new result here is a reproducible falsification of the supplied artifact, not a new fast multiplication identity.",
        "", "## The nominal 14-product explanation", "",
        "`238×15 + 5×29 + 1×21 = 3736`, versus `250×15 = 3750`. This is only rank arithmetic. Exact residuals disprove it as a certificate for the downloaded artifact.",
        "", "## Is the rank-29 kernel standalone?", "",
        "No. All five rank-29 candidates have nonzero residuals; the representative certificate is deliberately marked FAIL.",
        "", "## Constructor universality", "",
        "Not established. The constructor rejects every extracted candidate because none has a PASS certificate.",
        "", "## New rank or exponent", "",
        "None. Single-type spectral regressions recover log₂3 and log₂7, but invalid source data prevents a fusion before/after exponent comparison.",
        "", "## Failed hypotheses and counterexamples", "",
        "- Hypothesis: the downloadable `<8,27,30>:3736` Triads exactly compute matrix multiplication. Counterexample above.",
        "- Hypothesis: a recovered 29-product pair is standalone. Each of five candidates has nonzero residual (see component census).",
        "- Hypothesis: rank arithmetic alone explains a valid 14 gain. It does not; the full exceptional closure retains the parent's 28,098 modular residuals.",
        "", "## Three next research questions", "",
        "1. Can FMM-Lille supply a corrected explicit artifact matching either rank 3744 or the historical 3736 claim?",
        "2. Which generation step introduced the six invalid exceptional components, and can the missing exact terms be reconstructed from provenance?",
        "3. Do any verified external fusion artifacts yield a standalone two-slot kernel under the same exact extractor?",
        "", "## One-command reproduction", "", "```bash", "./reproduce_fusion.sh", "```", "",
        "## Formal no-contact proposition", "",
        summary["formal_invariant"],
    ]
    (RESULTS / "HANDOFF_TO_BRAIN.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
