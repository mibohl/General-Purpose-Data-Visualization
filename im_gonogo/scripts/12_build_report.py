"""Fill report placeholders from tables. Usage: python -I 12_build_report.py <work_out_dir> <deliverable_dir>"""
import sys
import numpy as np, pandas as pd
O, D = sys.argv[1], sys.argv[2]
sys.path.insert(0, __file__.rsplit('/', 1)[0])
head = open(f'{D}/report_head.md').read(); body = open(f'{D}/report_body.md').read()
step4 = open(f'{D}/report_step4.md').read()
t = lambda k: open(f'{D}/_tbl_{k}.md').read()
dec = pd.read_csv(f'{D}/decision_numbers.csv')
rows = []
for setting, g in dec.groupby('setting', sort=False):
    r = g[g.oracle == 'oracle_isotonic'].iloc[0]
    rng = f"{g.A_minus_O_over_A.min():.2f} to {g.A_minus_O_over_A.max():.2f}"
    rows.append(f"| {setting} | {r.F:.4g} | {r.D:.4g} | {r.A:.4g} | {r.A50:.4g} | {r.A_all_insample:.4g} | {r.O:.4g} | {r.A_over_F:.2f} | "
                f"{'pass' if r.crit1_pass else '**fail**'} | {r.A_minus_O_over_A:.2f} ({rng}) | {'pass' if r.crit2_pass else '**fail**'} | "
                f"{'NO-GO' if not (r.crit1_pass and r.crit2_pass) else 'go to Step 5'} |")
dec_md = ("| setting | F | D | A (affine, 20) | affine, 50 | affine, all (in-sample) | O (isotonic) | A/F | rule 1: A > 2F | "
          "(A−O)/A (range over 3 oracle variants) | rule 2: ≥ 0.25 | verdict |\n|" + "---|" * 12 + "\n" + "\n".join(rows) + "\n\n"
          "Units: Vs/cm² for 1/K0 rows, Å² for CCS rows. Oracle variants: in-sample isotonic (pre-registered), 5-fold cross-fitted isotonic, per-charge isotonic.")
p3 = pd.read_csv(f'{O}/step3_pooled.csv')
meth = ['raw', 'shift_20', 'affine_20', 'shift_50', 'affine_50', 'affine_all_insample', 'oracle_isotonic', 'oracle_isotonic_cv5', 'oracle_isotonic_per_charge']
lines = ['| space | subset | method | median abs dev | seed range | IQR |', '|---|---|---|---|---|---|']
for (sp, sub), g in p3.groupby(['space', 'subset']):
    g = g.set_index('method')
    for m in meth:
        r = g.loc[m]; f = '{:.4f}' if sp == 'k0' else '{:.2f}'
        sr = '' if r.med_seed_min == r.med_seed_max else f'{f.format(r.med_seed_min)} to {f.format(r.med_seed_max)}'
        lines.append(f"| {'1/K0' if sp == 'k0' else 'CCS'} | {sub} | {m} | {f.format(r.med)} | {sr} | {f.format(r.q25)} to {f.format(r.q75)} |")
rep = {'{{DECISION_TABLE}}': dec_md, '{{PROVENANCE}}': t('provenance'), '{{F_TABLE}}': t('F'), '{{FV_TABLE}}': t('Fv'),
       '{{STEP3_POOLED}}': '\n'.join(lines), '{{STEP3_TABLE}}': t('step3'), '{{STEP5_TABLE}}': t('step5'), '{{META_TABLE}}': t('meta'),
       '{{STEP4}}': step4}
out = head + body
for k, v in rep.items():
    assert k in out, k; out = out.replace(k, v)
hl = open(f'{D}/report_step4_headline.md').read().strip()
out = out.replace('{{STEP4_HEADLINE}}', hl)
assert '{{' not in out
open(f'{D}/report.md', 'w').write(out)
print('report.md written', len(out))
