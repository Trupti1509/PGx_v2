#!/usr/bin/env python3
"""
Generate a small Warfarin calculator input CSV from:
  1) the Step 1 VCF, and
  2) the PharmCAT 3.2.0 JSON report.

This DOES NOT calculate a warfarin dose.
It only captures:
- CYP2C9 / CYP4F2 / VKORC1 PharmCAT summary calls
- the nine warfarin-relevant marker calls
- marker availability / no-call status
- simple coverage flags for later calculator eligibility logic
"""

import csv
import json
import os
import re
import sys
from pathlib import Path

PHARMCAT_VERSION = "3.2.0"

WARFARIN_MARKERS = {
    "rs1799853": {"gene": "CYP2C9", "meaning": "CYP2C9 *2",  "group": "core"},
    "rs1057910": {"gene": "CYP2C9", "meaning": "CYP2C9 *3",  "group": "core"},
    "rs28371686": {"gene": "CYP2C9", "meaning": "CYP2C9 *5", "group": "expanded"},
    "rs9332131": {"gene": "CYP2C9", "meaning": "CYP2C9 *6",  "group": "expanded"},
    "rs7900194": {"gene": "CYP2C9", "meaning": "CYP2C9 *8",  "group": "expanded"},
    "rs28371685": {"gene": "CYP2C9", "meaning": "CYP2C9 *11","group": "expanded"},
    "rs9923231": {"gene": "VKORC1", "meaning": "VKORC1 -1639G>A", "group": "core"},
    "rs2108622": {"gene": "CYP4F2", "meaning": "CYP4F2 *3", "group": "additional"},
    "rs12777823": {"gene": "CYP2C cluster", "meaning": "rs12777823", "group": "additional"},
}

CORE_MARKERS = {"rs1799853", "rs1057910", "rs9923231"}
EXPANDED_CYP2C9_MARKERS = {"rs28371686", "rs9332131", "rs7900194", "rs28371685"}
ADDITIONAL_MARKERS = {"rs2108622", "rs12777823"}

RS_RE = re.compile(r"(rs\d+)", re.I)


def as_list(x):
    return x if isinstance(x, list) else []


def as_dict(x):
    return x if isinstance(x, dict) else {}


def clean_text(x):
    if x is None:
        return ""
    if isinstance(x, str):
        return " ".join(x.split()).strip()
    if isinstance(x, (int, float, bool)):
        return str(x)
    if isinstance(x, dict):
        for key in ("name", "label", "id", "rsid", "rsID", "hgvs", "variant", "gene"):
            val = x.get(key)
            if isinstance(val, str) and val.strip():
                return " ".join(val.split()).strip()
        return json.dumps(x, ensure_ascii=False, sort_keys=True)
    if isinstance(x, list):
        return "; ".join(clean_text(i) for i in x if clean_text(i))
    return str(x).strip()


def join_values(x, sep="; "):
    vals = []
    for item in as_list(x):
        val = clean_text(item)
        if val:
            vals.append(val)
    return sep.join(vals)


def extract_gene_summary(pharmcat_json, gene):
    info = as_dict(as_dict(pharmcat_json.get("genes")).get(gene))
    dips = as_list(info.get("recommendationDiplotypes"))
    if not dips:
        dips = as_list(info.get("sourceDiplotypes"))

    if not dips:
        return {
            "diplotype": "UNAVAILABLE",
            "phenotype": "UNAVAILABLE",
            "activity_score": "",
        }

    diplotypes = []
    phenotypes = []
    activities = []

    for dip in dips:
        dip = as_dict(dip)
        label = clean_text(dip.get("label"))
        phenotype = join_values(dip.get("phenotypes"))
        activity = clean_text(dip.get("activityScore"))

        if label:
            diplotypes.append(label)
        if phenotype:
            phenotypes.append(phenotype)
        if activity:
            activities.append(activity)

    return {
        "diplotype": " OR ".join(sorted(set(diplotypes))) if diplotypes else "UNAVAILABLE",
        "phenotype": "; ".join(sorted(set(phenotypes))) if phenotypes else "UNAVAILABLE",
        "activity_score": "; ".join(sorted(set(activities))),
    }


def decode_gt(ref, alts, gt):
    gt = (gt or "").replace("|", "/")
    if gt in {"", ".", "./."}:
        return "", "NO_CALL"

    allele_table = [ref] + alts
    idxs = gt.split("/")

    alleles = []
    for idx in idxs:
        if idx == ".":
            return "", "NO_CALL"
        try:
            i = int(idx)
        except ValueError:
            return "", "INVALID_GT"
        if i < 0 or i >= len(allele_table):
            return "", "INVALID_GT"
        alleles.append(allele_table[i])

    return "/".join(alleles), "CALLED"


def parse_warfarin_markers(vcf_path):
    found = {}

    with open(vcf_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue

            chrom, pos, id_field, ref, alt = parts[:5]
            sample_field = parts[9]

            rsids = [m.group(1).lower() for m in RS_RE.finditer(id_field)]
            wanted = [rs for rs in rsids if rs in WARFARIN_MARKERS]
            if not wanted:
                continue

            gt = sample_field.split(":", 1)[0].strip()
            alts = [] if alt in {"", "."} else alt.split(",")
            genotype, status = decode_gt(ref, alts, gt)

            for rsid in wanted:
                found[rsid] = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "vcf_gt": gt,
                    "genotype": genotype,
                    "status": status,
                }

    rows = {}
    for rsid, meta in WARFARIN_MARKERS.items():
        call = found.get(rsid)
        if call is None:
            rows[rsid] = {
                **meta,
                "chrom": "",
                "pos": "",
                "ref": "",
                "alt": "",
                "vcf_gt": "",
                "genotype": "",
                "status": "UNAVAILABLE",
            }
        else:
            rows[rsid] = {**meta, **call}

    return rows


def vkorc1_clinical_notation(marker_rows):
    row = marker_rows["rs9923231"]
    if row["status"] != "CALLED" or not row["genotype"]:
        return "UNAVAILABLE"

    alleles = row["genotype"].split("/")
    comp = {"C": "G", "T": "A", "G": "C", "A": "T"}

    try:
        converted = [comp[a] for a in alleles]
    except KeyError:
        return "UNAVAILABLE"

    if set(converted) == {"G", "A"}:
        return "G/A"
    return "/".join(converted)


def all_called(marker_rows, marker_set):
    return all(marker_rows[rs]["status"] == "CALLED" for rs in marker_set)


def generate_warfarin_csv(vcf_path, pharmcat_json_path, sample_id=None, output_dir="results/warfarin_inputs"):
    vcf_path = Path(vcf_path)
    pharmcat_json_path = Path(pharmcat_json_path)

    if not vcf_path.exists():
        raise FileNotFoundError(f"Step 1 VCF not found: {vcf_path}")
    if not pharmcat_json_path.exists():
        raise FileNotFoundError(f"PharmCAT JSON not found: {pharmcat_json_path}")

    if sample_id is None:
        sample_id = vcf_path.stem

    display_sample_id = sample_id.removeprefix("step1_")

    with open(pharmcat_json_path, encoding="utf-8") as fh:
        pharmcat = json.load(fh)

    marker_rows = parse_warfarin_markers(vcf_path)

    cyp2c9 = extract_gene_summary(pharmcat, "CYP2C9")
    cyp4f2 = extract_gene_summary(pharmcat, "CYP4F2")
    vkorc1 = extract_gene_summary(pharmcat, "VKORC1")

    core_complete = all_called(marker_rows, CORE_MARKERS)
    expanded_complete = all_called(marker_rows, EXPANDED_CYP2C9_MARKERS)
    additional_complete = all_called(marker_rows, ADDITIONAL_MARKERS)

    os.makedirs(output_dir, exist_ok=True)
    out_path = Path(output_dir) / f"Warfarin_Input_{display_sample_id}.csv"

    fieldnames = [
        "file_type",
        "schema_version",
        "sample_id",
        "pharmcat_version",
        "cyp2c9_diplotype",
        "cyp2c9_phenotype",
        "cyp2c9_activity_score",
        "cyp4f2_diplotype",
        "cyp4f2_phenotype",
        "vkorc1_pharmcat_diplotype",
        "vkorc1_clinical_1639",
        "core_complete",
        "expanded_cyp2c9_complete",
        "additional_markers_complete",
        "marker",
        "gene",
        "meaning",
        "group",
        "status",
        "vcf_gt",
        "genotype",
        "chrom",
        "pos",
        "ref",
        "alt",
    ]

    with open(out_path, "w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for rsid in WARFARIN_MARKERS:
            row = marker_rows[rsid]
            writer.writerow({
                "file_type": "XCODE_WARFARIN_PGX_INPUT",
                "schema_version": "1.0",
                "sample_id": display_sample_id,
                "pharmcat_version": PHARMCAT_VERSION,
                "cyp2c9_diplotype": cyp2c9["diplotype"],
                "cyp2c9_phenotype": cyp2c9["phenotype"],
                "cyp2c9_activity_score": cyp2c9["activity_score"],
                "cyp4f2_diplotype": cyp4f2["diplotype"],
                "cyp4f2_phenotype": cyp4f2["phenotype"],
                "vkorc1_pharmcat_diplotype": vkorc1["diplotype"],
                "vkorc1_clinical_1639": vkorc1_clinical_notation(marker_rows),
                "core_complete": str(core_complete).upper(),
                "expanded_cyp2c9_complete": str(expanded_complete).upper(),
                "additional_markers_complete": str(additional_complete).upper(),
                "marker": rsid,
                "gene": row["gene"],
                "meaning": row["meaning"],
                "group": row["group"],
                "status": row["status"],
                "vcf_gt": row["vcf_gt"],
                "genotype": row["genotype"],
                "chrom": row["chrom"],
                "pos": row["pos"],
                "ref": row["ref"],
                "alt": row["alt"],
            })

    if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
        print(f"[WARFARIN] Input CSV created: {out_path}")
        print(f"[WARFARIN] Core markers complete: {core_complete}")
        print(f"[WARFARIN] Expanded CYP2C9 markers complete: {expanded_complete}")
        print(f"[WARFARIN] Additional markers complete: {additional_complete}")

    return str(out_path)


def main():
    if len(sys.argv) not in (3, 4):
        print(
            "Usage:\n"
            "  python step2b_warfarin_input.py <step1.vcf> <pharmcat.report.json> [sample_id]"
        )
        sys.exit(1)

    vcf_path = sys.argv[1]
    json_path = sys.argv[2]
    sample_id = sys.argv[3] if len(sys.argv) == 4 else None

    generate_warfarin_csv(vcf_path, json_path, sample_id=sample_id)


if __name__ == "__main__":
    main()
