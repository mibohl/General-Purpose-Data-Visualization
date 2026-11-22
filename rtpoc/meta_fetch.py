"""Fetch PRIDE Archive v3 project records (JSON) for a list of PXD accessions.

Saves data/metadata/pride/<PXD>.json (skips existing). Only public metadata is
requested; nothing is sent but the accession.
Usage: python -I rtpoc/meta_fetch.py <file with one PXD per line or task names> [out subdir]
"""
import json, os, re, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "metadata", sys.argv[2] if len(sys.argv) > 2 else "pride")
os.makedirs(OUT, exist_ok=True)
pxds = sorted(set(re.findall(r"PXD\d{6}", open(sys.argv[1]).read())))
print(len(pxds), "accessions", flush=True)


def get(p):
    f = os.path.join(OUT, p + ".json")
    if os.path.exists(f):
        return p, "cached"
    for attempt in range(4):
        try:
            with urllib.request.urlopen(f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{p}",
                                        timeout=60) as r:
                d = json.load(r)
            json.dump(d, open(f, "w"))
            return p, "ok"
        except Exception as e:  # noqa: BLE001
            err = str(e)
            time.sleep(2 * (attempt + 1))
    return p, "fail: " + err


with ThreadPoolExecutor(4) as ex:
    res = list(ex.map(get, pxds))
fails = [r for r in res if r[1].startswith("fail")]
print("done", len(res), "failed", len(fails), fails[:10], flush=True)
