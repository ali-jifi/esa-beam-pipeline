# carry labels from a backup of candidate jsonls into regenerated ones
# join on candidate_id, fallback probe+t_center rounded to 0.1s
import argparse
import json
from pathlib import Path

CANDIDATES = Path(__file__).resolve().parents[1] / "candidates"
LAB = ("label", "label_source", "label_confidence")


def main():
    p = argparse.ArgumentParser(description="migrate labels after a batch regen")
    p.add_argument("--backup", default=str(CANDIDATES / "pre_regen"),
                   help="dir holding the labeled pre-regen jsonls")
    args = p.parse_args()
    old_dir = Path(args.backup)

    carried_id = carried_time = orphaned = 0
    for new_path in sorted(CANDIDATES.glob("*_survey.jsonl")):
        old_path = old_dir / new_path.name
        if not old_path.exists():
            continue
        old_by_id, old_by_time = {}, {}
        for line in open(old_path):
            r = json.loads(line)
            if r.get("label"):
                old_by_id[r["candidate_id"]] = r
                old_by_time[(r["probe"], round(r["t_center"], 1))] = r
        if not old_by_id:
            continue
        recs = [json.loads(l) for l in open(new_path)]
        matched = set()
        for r in recs:
            src = old_by_id.get(r["candidate_id"])
            via_time = False
            if src is None:
                src = old_by_time.get((r["probe"], round(r["t_center"], 1)))
                via_time = src is not None
            if src is None:
                continue
            for k in LAB:
                r[k] = src.get(k)
            matched.add(src["candidate_id"])
            carried_time += via_time
            carried_id += not via_time
        for cid, src in old_by_id.items():
            if cid not in matched:
                orphaned += 1
                print(f"[orphan] {cid} ({src.get('label')}) in {new_path.name}")
        tmp = new_path.with_suffix(".jsonl.tmp")
        with open(tmp, "w") as fh:
            for r in recs:
                fh.write(json.dumps(r) + "\n")
        tmp.replace(new_path)

    print(f"carried by id: {carried_id}, by time: {carried_time}, "
          f"orphaned: {orphaned}")


if __name__ == "__main__":
    main()
