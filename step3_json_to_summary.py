#!/usr/bin/env python3
import os
import sys
import glob
import json
import re
import html
import argparse
from datetime import datetime
import pandas as pd

REPORT_DIR = "results/reports"
OUT_DIR = "results"
DEFAULT_LOE = "LOE-CPIC.xlsx"

GENOTYPE_ONLY_GENES = {"CYP4F2", "VKORC1", "IFNL3"}

def _map_genotype_key(gene: str, diplotype_label: str) -> str:
    return str(diplotype_label).strip().lower()

def latest_report(d=REPORT_DIR):
    files = glob.glob(os.path.join(d, "*.report.json"))
    if not files:
        raise FileNotFoundError(f"No *.report.json in {d}")
    return max(files, key=os.path.getctime)

def load_json(p):
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)

def clean(x):
    if not isinstance(x, str):
        return x
    x = re.sub(r"<[^>]+>", "", x)
    return " ".join(html.unescape(x).split()).strip()

def as_list(x):
    return x if isinstance(x, list) else []

def as_dict(x):
    return x if isinstance(x, dict) else {}

def s(x):
    if x is None:
        return ""
    if isinstance(x, str):
        return clean(x)
    if isinstance(x, (int, float, bool)):
        return str(x)
    if isinstance(x, dict):
        for k in ("name", "label", "id", "rsid", "rsID", "hgvs", "variant", "gene"):
            if isinstance(x.get(k), str) and x[k].strip():
                return clean(x[k])
        return json.dumps(x, ensure_ascii=False, sort_keys=True)
    if isinstance(x, list):
        parts = [s(i) for i in x if s(i)]
        return "; ".join(parts)
    return clean(str(x))

def join(x, sep="; "):
    parts = []
    for i in as_list(x):
        v = s(i)
        if v:
            parts.append(v)
    return sep.join(parts)

def norm_drug(x):
    return (x or "").strip().lower()

def norm_gene(x):
    return (x or "").strip().upper()

def first_nonempty(d, keys):
    for k in keys:
        v = d.get(k, None)
        if isinstance(v, str) and clean(v):
            return clean(v)
        if v not in (None, "", [], {}):
            vv = s(v)
            if isinstance(vv, str) and vv.strip() and vv not in ("{}", "[]"):
                return vv.strip()
    return ""

def drop_all_empty_cols(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    out = out.replace(r"^\s*$", pd.NA, regex=True)
    return out.dropna(axis=1, how="all")

def fill_na_activity(df: pd.DataFrame, col="Activity Score"):
    if df.empty or col not in df.columns:
        return df
    out = df.copy()
    mask = out[col].astype(str).str.strip().isin(["", "nan"]) | out[col].isna()
    out.loc[mask, col] = "N/A"
    return out

def source_rank(source: str) -> int:
    t = (source or "").lower()
    if "cpic" in t:
        return 1
    if "dpwg" in t:
        return 2
    if "fda" in t:
        return 3
    return 9

def build_drug_citation_map(data: dict):
    m = {}
    for source, drug_dict in as_dict(data.get("drugs")).items():
        for drug, info in as_dict(drug_dict).items():
            info = as_dict(info)
            pmids = []
            for cite in as_list(info.get("citations")):
                cite = as_dict(cite)
                pmid = cite.get("pmid")
                if pmid:
                    val = str(pmid).strip()
                    if val:
                        pmids.append(val)
            pmids = "; ".join(sorted(set(pmids)))
            m[(source, norm_drug(drug))] = pmids
    return m

def extract_genotypes(data):
    rows = []
    gene_summary = {}

    for gene, info in as_dict(data.get("genes")).items():
        info = as_dict(info)
        rec_dips = as_list(info.get("recommendationDiplotypes"))
        dips = rec_dips if rec_dips else as_list(info.get("sourceDiplotypes"))

        msgs = join(info.get("messages"))
        uncalled = join(info.get("uncalledHaplotypes"))
        undoc = s(info.get("hasUndocumentedVariations"))
        voi = join(info.get("variantsOfInterest"))

        drugs = ", ".join(
            [
                as_dict(d).get("name", "")
                for d in as_list(info.get("relatedDrugs"))
                if as_dict(d).get("name")
            ]
        )

        if not dips:
            gene_summary[gene] = {
                "phenotype": "Indeterminate",
                "diplotype": "No Call",
                "activity": "",
            }
            rows.append(
                {
                    "Gene": gene,
                    "Diplotype": "No Call",
                    "Phenotype": "Indeterminate",
                    "Activity Score": "",
                    "Related Drugs": drugs,
                    "Variants Of Interest": voi,
                    "Messages": msgs,
                    "Uncalled": uncalled,
                    "Undocumented": undoc,
                }
            )
            continue

        phenos_for_gene = []
        dips_for_gene = []
        acts_for_gene = []

        for dip in dips:
            dip = as_dict(dip)
            raw_p = join(dip.get("phenotypes")) or ""
            d_lbl = s(dip.get("label", "Unknown"))
            a = s(dip.get("activityScore", ""))

            if gene in GENOTYPE_ONLY_GENES:
                geno_key = _map_genotype_key(gene, d_lbl)
                p = geno_key if geno_key else (raw_p or "Indeterminate")
            else:
                p = raw_p or "Indeterminate"

            # Preserve original PharmCAT phenotype (before genotype-key override)
            # so genotype summary can display the actual clinical phenotype.
            pharmcat_pheno = raw_p or "Indeterminate"

            a1 = as_dict(dip.get("allele1"))
            a2 = as_dict(dip.get("allele2"))

            rows.append(
                {
                    "Gene": gene,
                    "Diplotype": d_lbl,
                    "Phenotype": p,
                    "PharmCAT Phenotype": pharmcat_pheno,
                    "Allele 1": s(a1.get("name")),
                    "Allele 1 Func": s(a1.get("function")),
                    "Allele 2": s(a2.get("name")),
                    "Allele 2 Func": s(a2.get("function")),
                    "Activity Score": a,
                    "Related Drugs": drugs,
                    "Variants Of Interest": voi,
                    "Messages": msgs,
                    "Uncalled": uncalled,
                    "Undocumented": undoc,
                }
            )
            phenos_for_gene.append(p)
            dips_for_gene.append(d_lbl)
            acts_for_gene.append(a)

        gene_summary[gene] = {
            "phenotype": "; ".join(sorted(set([x for x in phenos_for_gene if x]))) or "Indeterminate",
            "diplotype": " OR ".join(sorted(set([x for x in dips_for_gene if x]))) or "Unknown",
            "activity": "; ".join(sorted(set([x for x in acts_for_gene if x]))),
        }

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("Gene", kind="stable")
    return df, gene_summary

def resolve_genes_for_annotation(ann, guide, gene_summary):
    phen = ann.get("phenotypes")
    if isinstance(phen, dict) and phen:
        return sorted([g for g in phen.keys() if g])

    raw = (guide.get("geneSymbol") or guide.get("name") or "")
    toks = re.findall(r"\b[A-Z0-9]{2,}\b", str(raw).upper())
    toks = [t for t in toks if t not in ("CPIC", "DPWG", "FDA", "GUIDELINE")]
    toks = [t for t in toks if t in gene_summary]
    return toks or ["UNKNOWN"]

def extract_recs(data, gene_summary):
    rows = []
    drug_cite_map = build_drug_citation_map(data)

    for source, drug_dict in as_dict(data.get("drugs")).items():
        for drug, info in as_dict(drug_dict).items():
            for guide in as_list(as_dict(info).get("guidelines")):
                guide = as_dict(guide)
                for ann in as_list(guide.get("annotations")):
                    ann = as_dict(ann)
                    rec = clean(ann.get("drugRecommendation", ""))
                    genes = resolve_genes_for_annotation(ann, guide, gene_summary)

                    dosing_struct = first_nonempty(ann, ["dosingInformation", "dosingRecommendation", "doseRecommendation", "dosingGuidance"])
                    alt_struct = first_nonempty(ann, ["alternativeDrug", "alternateDrug", "alternativeDrugs", "alternateDrugs"])
                    other_struct = first_nonempty(ann, ["otherPrescribingGuidance", "prescribingGuidance", "warnings", "testingRecommendations", "commentary", "notes"])

                    has_dose = bool(ann.get("dosingInformationAvailable")) or bool(dosing_struct)
                    has_alt = bool(ann.get("alternativeDrugAvailable")) or bool(ann.get("alternateDrugAvailable")) or bool(alt_struct)
                    has_other = bool(ann.get("otherPrescribingGuidanceAvailable")) or bool(other_struct)

                    dosing = dosing_struct if has_dose else ""
                    alt_drug = alt_struct if has_alt else ""
                    other = other_struct if has_other else ""

                    if has_dose and not dosing: dosing = "See Recommendation"
                    if has_alt and not alt_drug: alt_drug = "See Recommendation"
                    if has_other and not other: other = "See Recommendation"

                    for gene in genes:
                        g_info = gene_summary.get(gene, {})
                        ph_dict = as_dict(ann.get("phenotypes"))

                        rows.append({
                            "Drug": drug,
                            "Gene": gene,
                            "Patient Phenotype": s(ph_dict.get(gene)) or g_info.get("phenotype", ""),
                            "Patient Diplotype": g_info.get("diplotype", ""),
                            "Activity Score": g_info.get("activity", ""),
                            "Source": source,
                            "Source Priority": source_rank(source),
                            "Classification": s(ann.get("classification")),
                            "Implication": join(ann.get("implications")),
                            "Recommendation": rec,
                            "Tag: Dosing Info": "Yes" if has_dose else "No",
                            "Tag: Alternate Drug": "Yes" if has_alt else "No",
                            "Tag: Other Guidance": "Yes" if has_other else "No",
                            "Dosing Info": dosing,
                            "Alternative Drug": alt_drug,
                            "Other Guidance": other,
                            "Guideline Name": clean(guide.get("name", "")),
                            "Guideline URL": guide.get("url", ""),
                            "Drug Citation PMIDs": drug_cite_map.get((source, norm_drug(drug)), ""),
                        })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["Drug", "Gene", "Source Priority"], kind="stable").drop(columns=["Source Priority"])
    return df

def extract_citations(data):
    rows, seen = [], set()

    def add(context, source, pmid, title="", journal="", year="", url=""):
        if not pmid: return
        pmid = str(pmid).strip()
        if not pmid or pmid in seen: return
        seen.add(pmid)
        rows.append({
            "Context": context,
            "Source": source or "",
            "PMID": pmid,
            "Title": clean(title) if title else "",
            "Journal": journal or "",
            "Year": year or "",
            "PubMed Link": url or f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
        })

    for source, drug_dict in as_dict(data.get("drugs")).items():
        for drug, info in as_dict(drug_dict).items():
            info = as_dict(info)
            for cite in as_list(info.get("citations")):
                cite = as_dict(cite)
                add(
                    context=f"{drug}",
                    source=source,
                    pmid=cite.get("pmid"),
                    title=cite.get("title", ""),
                    journal=cite.get("journal", ""),
                    year=cite.get("year", ""),
                    url=cite.get("url", ""),
                )

    for ref in as_list(data.get("references")):
        ref = as_dict(ref)
        add(
            context="PharmCAT Method",
            source="PharmCAT",
            pmid=ref.get("pmid"),
            title=ref.get("title", ""),
            journal=ref.get("journal", ""),
            year=ref.get("year", ""),
            url=ref.get("url", ""),
        )

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("Year", ascending=False, kind="stable")
    return df

def load_cpic_loe(path):
    if not os.path.exists(path):
        return pd.DataFrame()

    df = pd.read_excel(path, sheet_name="CPIC Gene-Drug Pairs")
    if "Drug" not in df.columns or "Gene" not in df.columns:
        return pd.DataFrame()

    rename = {}
    for c in df.columns:
        cl = c.lower().strip()
        if cl == "cpic level": rename[c] = "CPIC Level"
        elif "pharmgkb level of evidence" in cl: rename[c] = "PharmGKB LoE"
        elif "pgx on fda label" in cl: rename[c] = "PGx on FDA Label"
        elif cl == "guideline": rename[c] = "CPIC Guideline URL"
        elif "cpic publications" in cl: rename[c] = "CPIC PMIDs"
    df = df.rename(columns=rename)

    df["_d"] = df["Drug"].map(norm_drug)
    df["_g"] = df["Gene"].map(norm_gene)

    cols = ["_d", "_g", "CPIC Level", "PharmGKB LoE", "PGx on FDA Label", "CPIC Guideline URL", "CPIC PMIDs"]
    cols = [c for c in cols if c in df.columns]
    return df[cols].drop_duplicates(subset=["_d", "_g"])

def attach_evidence(df_recs, loe):
    if df_recs.empty or loe.empty:
        return df_recs

    out = df_recs.copy()
    out["_d"] = out["Drug"].map(norm_drug)
    out["_g"] = out["Gene"].map(norm_gene)

    merged = out.merge(loe, how="left", on=["_d", "_g"]).drop(columns=["_d", "_g"], errors="ignore")

    is_cpic = merged["Source"].astype(str).str.contains("CPIC", case=False, na=False)
    for c in ["CPIC Level", "CPIC Guideline URL", "CPIC PMIDs"]:
        if c in merged.columns:
            merged.loc[~is_cpic, c] = ""
    return merged

def get_actionable(df):
    if df.empty:
        return df

    cls = df["Classification"].fillna("").astype(str)
    imp = df["Implication"].fillna("").astype(str)
    rec = df["Recommendation"].fillna("").astype(str)

    mask_strong = cls.str.contains(r"Strong|Moderate|Actionable", case=False, na=False)
    mask_imp_ok = ~imp.str.contains(r"normal|standard|no recommendation", case=False, na=False)
    mask_rec_ok = ~rec.str.contains(r"initiate standard dosing|recommend standard dosing", case=False, na=False)

    out = df[mask_strong & mask_imp_ok & mask_rec_ok].copy()

    def rank(x):
        x = str(x).lower()
        if "strong" in x: return 1
        if "moderate" in x: return 2
        if "actionable" in x: return 3
        return 99

    out["_rank"] = out["Classification"].map(rank)
    out = out.sort_values(["_rank", "Drug", "Gene"], kind="stable").drop(columns=["_rank"])

    cols = [
        "Drug", "Gene", "Patient Phenotype", "Patient Diplotype", "Activity Score",
        "Classification", "CPIC Level", "PharmGKB LoE", "PGx on FDA Label",
        "Implication", "Recommendation", "Tag: Dosing Info", "Dosing Info",
        "Tag: Alternate Drug", "Alternative Drug", "Tag: Other Guidance", "Other Guidance",
        "Guideline URL", "Drug Citation PMIDs",
    ]
    cols = [c for c in cols if c in out.columns]
    return out[cols]

def get_qc(df):
    if df.empty:
        return df
    d = df.copy()
    for col in ["Messages", "Uncalled"]:
        if col not in d.columns:
            d[col] = ""
    bad = d[
        d["Diplotype"].astype(str).str.contains("No Call", case=False, na=False) |
        d["Phenotype"].astype(str).str.contains("Indeterminate", case=False, na=False) |
        (d["Messages"].astype(str).str.strip() != "") |
        (d["Uncalled"].astype(str).str.strip() != "")
    ].copy()
    cols = ["Gene", "Diplotype", "Phenotype", "Messages", "Uncalled", "Undocumented", "Variants Of Interest"]
    cols = [c for c in cols if c in bad.columns]
    return bad[cols]

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", default="", help="Input *.report.json (default: latest)")
    ap.add_argument("--loe", default=DEFAULT_LOE, help="CPIC LOE Excel")

    if argv is None:
        args = ap.parse_args()
    else:
        args = ap.parse_args(argv)

    report = args.report.strip() or latest_report()
    print(f"[STEP 3] Processing: {report}")

    data = load_json(report)

    df_geno, gene_summary = extract_genotypes(data)
    df_recs = extract_recs(data, gene_summary)
    df_cites = extract_citations(data)

    loe = load_cpic_loe(args.loe)
    df_recs = attach_evidence(df_recs, loe)

    df_recs = fill_na_activity(df_recs, "Activity Score")
    df_act = get_actionable(df_recs)
    df_act = fill_na_activity(df_act, "Activity Score")
    df_qc = get_qc(df_geno)

    df_geno = drop_all_empty_cols(df_geno)
    df_recs = drop_all_empty_cols(df_recs)
    df_act = drop_all_empty_cols(df_act)
    df_qc = drop_all_empty_cols(df_qc)
    df_cites = drop_all_empty_cols(df_cites)

    sample = os.path.basename(report).replace(".report.json", "")
    out_file = os.path.join(OUT_DIR, f"step3_{sample}_MASTER.xlsx")
    os.makedirs(OUT_DIR, exist_ok=True)

    with pd.ExcelWriter(out_file, engine="xlsxwriter") as w:
        def write(df, name, color):
            if df is None or df.empty:
                df = pd.DataFrame([{"Status": "No Data"}])
            df.to_excel(w, sheet_name=name, index=False)
            ws = w.sheets[name]
            fmt = w.book.add_format({"bold": True, "bg_color": color, "font_color": "white", "border": 1})
            wrap = w.book.add_format({"text_wrap": True, "valign": "top"})
            for i, col in enumerate(df.columns):
                ws.write(0, i, col, fmt)
                width = 18
                if any(k in col.lower() for k in ["recommendation", "implication", "guidance", "dosing", "alternative", "phenotype"]):
                    width = 55
                if "url" in col.lower() or "link" in col.lower():
                    width = 35
                if col.lower() == "title":
                    width = 70
                ws.set_column(i, i, width, wrap)
            ws.freeze_panes(1, 0)

        meta = pd.DataFrame([{
            "Generated": datetime.now().isoformat(timespec="seconds"),
            "JSON Report": report,
            "LOE File": args.loe if os.path.exists(args.loe) else "NOT FOUND",
            "Rows: Actionable": len(df_act),
            "Rows: All Recommendations": len(df_recs),
            "Rows: Genotype Details": len(df_geno),
            "Rows: QC Flags": len(df_qc),
            "Rows: Citations": len(df_cites),
        }])
        meta.to_excel(w, sheet_name="0_METADATA", index=False)

        write(df_act, "1_ACTIONABLE", "#C00000")
        write(df_recs, "2_ALL_RECOMMENDATIONS", "#2F5597")
        write(df_geno, "3_GENOTYPE_DETAILS", "#548235")
        write(df_qc, "4_QC_FLAGS", "#BF8F00")
        write(df_cites, "5_CITATIONS", "#7030A0")

    print(f"[STEP 3] Success: {out_file}")

if __name__ == "__main__":
    main()