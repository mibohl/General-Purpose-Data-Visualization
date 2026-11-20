"""Dataset list and easy/hard labels, fixed in data/datasets.md before any experiment."""

EASY = [
    "hela_hf_psms_aligned",
    "hela_lumos_2h_psms_aligned",
    "pancreas_psms_aligned",
    "yeast_120min_psms_aligned",
    "yeast_60min_psms_aligned",
    "arabidopsis_psms_aligned",
    "proteometools_library",
    "dia_fixed_mods",
    "unmod_fixed_mods",
]
HARD = [
    "ATLANTIS_SILICA_fixed_mods",
    "LUNA_SILICA_fixed_mods",
    "LUNA_HILIC_fixed_mods",
    "SCX_fixed_mods",
    "Xbridge_fixed_mods",
    "prosit_ptm_2020",
    "mod_fixed_mods",
]
ALL = EASY + HARD
GROUP = {**{d: "easy" for d in EASY}, **{d: "hard" for d in HARD}}
