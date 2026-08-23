import pandas as pd
from heartdelta.analysis import complete_level_cohort


def test_complete_level_qc_keeps_only_all_three_levels():
    rows=[]
    for case,levels in (("ok",("basal","mid","apical")),("missing",("basal","mid"))):
        for level in levels:
            rows.append({"case_id":case,"timepoint":"fu1","level":level,"segment":pd.NA,"hu_median":50,"fraction_hu_0_100":.9,"baseline_hu_median":50,"baseline_fraction_hu_0_100":.9})
    kept,audit=complete_level_cohort(pd.DataFrame(rows))
    assert set(kept.case_id)=={"ok"}
    assert audit.set_index("case_id").loc["missing","complete_three_level_qc"] == False
