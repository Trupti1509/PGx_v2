# STEP 1: Convert raw DNA to VCF
#!/usr/bin/env python3
import os
import sys
import time
import re
from datetime import datetime

RS_REGEX = re.compile(r'^(rs\d+)')

PHARMCAT_VCF = "pharmcat_data/pharmcat_positions_3.2.0.vcf"
OUTPUT_DIR = "results"

# ── Supplemental SNPs (not in pharmcat_positions.vcf) ─────────────────────────
# These are clinically relevant PGx variants that PharmCAT does not track but
# that appear in consumer DNA raw files.  They are extracted alongside PharmCAT
# targets, written to the step1 VCF, and later read by step4 to build synthetic
# gene rows for the GSI merge.

SUPPLEMENTAL_SNPS: dict = {
    "rs1799963": {"chrom": "chr11", "pos": 46761055, "ref": "G", "alt": "A"},
    "rs6025":    {"chrom": "chr1",  "pos": 169549811, "ref": "C", "alt": "T"},
}

CHROM_ORDER = {
    **{f"chr{i}": i for i in range(1, 23)},
    "chrX": 23,
    "chrY": 24,
    "chrM": 25,
    "chrMT": 25,
}

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
BAD_ALLELES = {"-", ".", "0", None, ""}

def die(msg, code=1):
    print(msg, file=sys.stderr)
    sys.exit(code)

def _valid_alleles(a1, a2):
    return a1 not in BAD_ALLELES and a2 not in BAD_ALLELES

def _is_header_or_meta(s: str) -> bool:
    low = s.lower()
    if s.startswith("#"):
        return True
    if "markername" in low and "chrom" in low and "pos" in low and "gt" in low:
        return True
    if "rsid" in low and ("chromosome" in low or "position" in low):
        return True
    return False

def load_pharmcat_positions(path):
    # print(f"[INFO] Loading PharmCAT reference: {path}")
    if not os.path.exists(path):
        bgz = os.path.join(
            os.path.dirname(path) or ".",
            "pharmcat_positions_3.2.0.vcf.bgz",
        )
        if os.path.exists(bgz):
            die(f"[ERROR] PharmCAT reference VCF missing.\nRun: gunzip -c {bgz} > {path}")
        die(f"[ERROR] Missing PharmCAT reference: {path}")

    meta_lines = []
    targets = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("##"):
                if not line.startswith("##fileformat="):
                    meta_lines.append(line.rstrip("\n"))
                continue
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split()
            if len(parts) < 5:
                continue
            chrom, pos, rsid, ref, alt = parts[:5]
            try:
                pos_i = int(pos)
            except ValueError:
                continue

            rsids = [r.strip() for r in rsid.split(";")] if ";" in rsid else [rsid]
            for sid in rsids:
                if sid.startswith("rs"):
                    targets[sid] = {
                        "chrom": chrom,
                        "pos": pos_i,
                        "ref": ref,
                        "alt": alt,
                    }

    # print(f"[INFO] Loaded {len(targets)} PharmCAT target rsIDs from {path}")

    # Inject supplemental SNPs that are not tracked by PharmCAT but are needed
    added = []
    for rsid, info in SUPPLEMENTAL_SNPS.items():
        if rsid not in targets:
            targets[rsid] = info
            added.append(rsid)
    if added:
        pass
        # print(f"[INFO] Injected {len(added)} supplemental SNP(s) into target list: {', '.join(added)}")
    return meta_lines, targets

def compute_vcf_gt(a1, a2, ref, alt_str):
    a1, a2, ref, alt_str = a1.upper(), a2.upper(), ref.upper(), alt_str.upper()

    if a1 in ("I", "D") or a2 in ("I", "D"):
        if len(alt_str) == len(ref):
            return None
        alt_ins = len(alt_str) > len(ref)

        def map_id(x):
            if x == "I":
                return "1" if alt_ins else "0"
            if x == "D":
                return "0" if alt_ins else "1"
            return None

        i1, i2 = map_id(a1), map_id(a2)
        if i1 is None or i2 is None:
            return None
        return "/".join(sorted((i1, i2)))

    alleles = [ref] + alt_str.split(",")

    def idx(b):
        return str(alleles.index(b)) if b in alleles else None

    def try_pair(x, y):
        i1, i2 = idx(x), idx(y)
        return "/".join(sorted((i1, i2))) if i1 and i2 else None

    gt = try_pair(a1, a2)
    if gt:
        return gt

    c1, c2 = COMPLEMENT.get(a1), COMPLEMENT.get(a2)
    return try_pair(c1, c2) if c1 and c2 else None

def parse_raw_line(line):
    s = line.strip()
    if not s or _is_header_or_meta(s):
        return None, None, None

    sep = "\t" if "\t" in s else ","
    parts = [p.strip().strip('"').strip("'") for p in s.split(sep)]
    if len(parts) < 2:
        return None, None, None

    rsid = None
    for p in parts:
        if p.startswith("rs") and "_" not in p and ":" not in p:
            _m = RS_REGEX.match(p)
            if _m:
                rsid = _m.group(1)
            break
    if not rsid:
        return None, None, None

    a1 = a2 = None

    if len(parts) >= 4:
        gt_field = parts[3].upper()
        if "/" in gt_field and len(gt_field) == 3:
            x, y = gt_field.split("/")
            if len(x) == 1 and len(y) == 1 and x.isalpha() and y.isalpha():
                a1, a2 = x, y
        elif len(parts) >= 5:
            c3, c4 = parts[3].upper(), parts[4].upper()
            if len(c3) == 1 and len(c4) == 1 and c3.isalpha() and c4.isalpha():
                a1, a2 = c3, c4

    if a1 is None:
        last = parts[-1].upper()
        if len(last) == 2 and last.isalpha():
            a1, a2 = last[0], last[1]
        elif len(last) == 1 and last.isalpha():
            a1 = a2 = last

    if not _valid_alleles(a1, a2):
        return None, None, None

    return rsid, a1, a2

def _print_progress(fh, size, total_lines, matched):
    pct = (fh.tell() / size) * 100 if size > 0 else 0
    sys.stdout.write(f"\r[PROGRESS] {pct:5.1f}%  lines={total_lines}  matches={matched}")
    sys.stdout.flush()

def convert_raw_to_vcf(raw_path, meta_lines, targets, out_path, sample_id):
    t0 = time.time()
    matches = []
    total_lines = matched = 0

    size = os.path.getsize(raw_path)
    # print(f"[INFO] Raw file: {raw_path} ({size/1024/1024:.2f} MB)")
    # print("[INFO] Scanning and matching variants")

    with open(raw_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            total_lines += 1
            rsid, a1, a2 = parse_raw_line(line)
            if not rsid:
                continue
            info = targets.get(rsid)
            if not info:
                continue

            gt = compute_vcf_gt(a1, a2, info["ref"], info["alt"])
            if not gt:
                continue

            matches.append(
                {
                    "chrom": info["chrom"],
                    "pos": info["pos"],
                    "id": rsid,
                    "ref": info["ref"],
                    "alt": info["alt"],
                    "gt": gt,
                }
            )
            matched += 1

            if total_lines % 10_000 == 0:
                _print_progress(f, size, total_lines, matched)

        _print_progress(f, size, total_lines, matched)
        sys.stdout.write("\n")

    if not matches:
        die("[ERROR] No variants matched. Check raw file format and contents.")

    # print(f"[INFO] Matched {len(matches)} variants. Sorting.")
    matches.sort(key=lambda m: (CHROM_ORDER.get(m["chrom"], 99), m["pos"]))

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    # print(f"[INFO] Writing VCF: {out_path}")

    with open(out_path, "w", encoding="ascii", newline="\n") as out:
        out.write("##fileformat=VCFv4.2\n")
        for h in meta_lines:
            if not (h.startswith("##fileformat=") or h.startswith("##FORMAT=<ID=GT")):
                out.write(h + "\n")

        out.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
        out.write("##source=raw_to_pharmcat_vcf_multi\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write(f"##sample={sample_id}\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_id + "\n")

        for m in matches:
            out.write(f"{m['chrom']}\t{m['pos']}\t{m['id']}\t{m['ref']}\t{m['alt']}\t.\tPASS\t.\tGT\t{m['gt']}\n")

    elapsed = time.time() - t0


def run_step1(raw_path: str):
    if not os.path.exists(raw_path):
        die(f"[ERROR] Input file missing: {raw_path}")
    meta_lines, targets = load_pharmcat_positions(PHARMCAT_VCF)
    sample_id = os.path.splitext(os.path.basename(raw_path))[0]
    out_path = os.path.join(OUTPUT_DIR, f"step1_{sample_id}.vcf")
    convert_raw_to_vcf(raw_path, meta_lines, targets, out_path, sample_id)

def main():
    if len(sys.argv) < 2:
        die("Usage: python3 step1_convert_rawdna_to_vcf.py <raw_dna_file>")
    raw_path = sys.argv[1]
    run_step1(raw_path)

if __name__ == "__main__":
    main()