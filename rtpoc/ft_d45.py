"""DeepLC 4.5 fine-tuning through its shipped transfer-learning API.

`deeplc.finetune_and_predict` = `deeplc.finetune` (new low-rank setup head,
rank + 2 = 66 parameters, encoder and pretrained setups frozen, early stopping
on a validation set, best weights restored) followed by `deeplc.calibrate` on
the reference (MultiHeadRidgeCalibration over the single fine-tuned column).
`finetune_and_predict` does not take a validation set, so this module calls the
same two public functions in the same order, passing our own 15% validation
split of the context to `finetune(psm_list_validation=...)`. No training loop
of our own.
"""
from __future__ import annotations

import logging
import math
import time

import numpy as np
import pandas as pd
import torch
from psm_utils import PSM, PSMList
from psm_utils.io.peptide_record import peprec_to_proforma

import deeplc
from deeplc import _model_ops
from deeplc.core import DEFAULT_MODEL

VAL_FRAC = 0.15

# count validation passes so we can report epochs actually run
_EPOCHS = {"n": 0}
_orig_validate = _model_ops._validate_epoch


def _counting_validate(*a, **k):
    _EPOCHS["n"] += 1
    return _orig_validate(*a, **k)


_model_ops._validate_epoch = _counting_validate
for name in ("deeplc", "deeplc.core", "deeplc._model_ops", "deeplc.calibration"):
    logging.getLogger(name).setLevel(logging.ERROR)


def to_psms(df: pd.DataFrame, idx, with_rt: bool) -> PSMList:
    rows = df.iloc[np.asarray(idx)]
    return PSMList(psm_list=[
        PSM(peptidoform=peprec_to_proforma(s, m if m else "-"), spectrum_id=str(i),
            retention_time=float(t) if with_rt else None)
        for i, s, m, t in zip(rows.index, rows.seq, rows.modifications, rows.tr)])


def val_split(ctx, seed):
    """Fixed 15% validation split of the context (depends on seed and context only)."""
    perm = np.random.default_rng(10_000 + seed).permutation(len(ctx))
    nv = max(2, math.ceil(VAL_FRAC * len(ctx)))
    return ctx[perm[nv:]], ctx[perm[:nv]]


def finetune_predict(df, ctx, predict_idx, seed, cfg, threads=4):
    """Fine-tune on ctx (85/15 train/val), calibrate on all of ctx, predict predict_idx.

    Returns calibrated predictions, raw fine-tuned-head predictions, seconds for
    the fine-tune itself, seconds in total, and epochs run.
    """
    torch.manual_seed(seed)
    torch.set_num_threads(threads)
    tr_idx, va_idx = val_split(np.asarray(ctx), seed)
    ref_all = to_psms(df, ctx, True)
    t0 = time.time()
    _EPOCHS["n"] = 0
    model = deeplc.finetune(
        psm_list_reference=to_psms(df, tr_idx, True),
        psm_list_validation=to_psms(df, va_idx, True),
        model=DEFAULT_MODEL,
        train_kwargs=dict(epochs=cfg["epochs"], learning_rate=cfg["lr"],
                          patience=cfg["patience"], show_progress=False,
                          num_threads=threads, device="cpu"),
    )
    t_ft = time.time() - t0
    epochs_run = max(0, _EPOCHS["n"] - 1)  # first validation scores the start point
    cal = deeplc.calibrate(psm_list_reference=ref_all, model=model)
    raw = np.asarray(deeplc.predict(to_psms(df, predict_idx, False), model=model,
                                    return_matrix=True), dtype=float)
    pred = np.asarray(cal.transform(raw), dtype=float)
    return pred, raw.reshape(len(predict_idx), -1)[:, 0], t_ft, time.time() - t0, epochs_run
