from pathlib import Path

from germsynth_cr.contact_compiler import analyze_rank3744_preconditions
from germsynth_cr.exact_tensor import verify_scheme
from germsynth_cr.maple_lrp_parser import parse_lrp
from germsynth_cr.maple_raw_parser import parse_raw
from germsynth_cr.maple_tensor_parser import parse_tensor
from germsynth_cr.pan_pair import build as build_pan, verify as verify_pan
from germsynth_cr.residual import compute_residual
from germsynth_cr.residual_localization import flattening_ranks
from germsynth_fusion.scheme_io import load_scheme


ROOT=Path(__file__).resolve().parents[1];FMM=ROOT/"external"/"fmm-lille";CURATED=ROOT/"external"/"curated"


def test_independent_parsers_on_controls():
    for name,rank in (("2x2x2",7),("2x3x3",15),("3x3x4",29)):
        shape=tuple(map(int,name.split("x")))
        assert verify_scheme(parse_tensor(FMM/name/f"{name}_tensor.mpl.bz2"))["status"]=="PASS"
        assert verify_scheme(parse_lrp(FMM/name/f"{name}_LRP.mpl.bz2",shape))["status"]=="PASS"
        raw=parse_raw(FMM/name/f"{name}_raw.mpl.bz2")
        assert raw.rank==rank and verify_scheme(raw)["status"]=="PASS"


def test_parent_representations_and_raw_semantics():
    tensor=parse_tensor(FMM/"8x27x30"/"8x27x30_tensor.mpl.bz2")
    lrp=parse_lrp(FMM/"8x27x30"/"8x27x30_LRP.mpl.bz2",(8,27,30))
    assert tensor.U==lrp.U and tensor.V==lrp.V and tensor.W==lrp.W
    assert verify_scheme(tensor)["modular_screens"][0]["residual_count"]==28098
    raw=parse_raw(FMM/"8x27x30"/"8x27x30_raw.mpl.bz2")
    assert raw.rank==3825 and verify_scheme(raw)["status"]=="PASS"


def test_21_component_eight_term_repair_is_flattening_impossible():
    data=__import__("json").loads((ROOT/"results"/"residual_autopsy.json").read_text())
    assert data["truncation_8"]["flattening_ranks"]==[6,9,6]
    assert data["truncation_8"]["eight_term_completion_impossible"] is True


def test_rank3744_preconditions_fail_on_pinned_base():
    outer=parse_tensor(FMM/"4x9x10"/"4x9x10_tensor.mpl.bz2")
    inner=parse_tensor(FMM/"2x3x3"/"2x3x3_tensor.mpl.bz2")
    concat=load_scheme(CURATED/"4x3x3-r29-kauers_2026-367ebc4.json")
    result=analyze_rank3744_preconditions(outer,inner,concat)
    assert len(result["factor_direction_pairs"]["U"])==6
    assert result["factor_direction_pairs"]["V"]==[]
    assert result["preconditions_met"] is False


def test_pan_profitability_boundary():
    assert verify_pan(build_pan(4,2,4))["gain_vs_naive"]==0
    result=verify_pan(build_pan(5,2,5))
    assert result["status"]=="PASS" and result["gain_vs_naive"]==5 and result["mixed_product_count"]>0
