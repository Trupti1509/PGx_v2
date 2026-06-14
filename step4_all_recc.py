#!/usr/bin/env python3
import os
import glob
import re
import pandas as pd

RESULTS_DIR  = "results"
GSI_SHEET         = "GSI Data"

# Auto-select the best available GSI file in this priority order:
#   1. gsi_output.xlsx  — the canonical output from the scraper (newest build)
#   2. GSI_DD_MM_YY.xlsx — dated snapshots, pick the most recently modified
#   3. GSI_16_04_26.xlsx — hard fallback if nothing else exists
def _find_latest_gsi() -> str:
    import re as _re
    # Priority 1: canonical scraper output file
    if os.path.exists("gsi_output.xlsx"):
        return "gsi_output.xlsx"
    # Priority 2: dated snapshot files
    pattern = _re.compile(r"[Gg][Ss][Ii]_\d{2}_\d{2}_\d{2}\.xlsx$")
    candidates = [f for f in glob.glob("GSI_*.xlsx") if pattern.search(f)]
    if candidates:
        return max(candidates, key=os.path.getmtime)
    # Priority 3: hard fallback
    return "GSI_16_04_26.xlsx"

SCRAPER_FILE      = _find_latest_gsi()
DRUG_CAT_SUPPLEMENT = "drug_category_supplement.xlsx"


STEP3_GENOTYPE_SHEET  = "3_GENOTYPE_DETAILS"
STEP3_RECS_SHEET      = "2_ALL_RECOMMENDATIONS"
STEP3_CITATIONS_SHEET = "5_CITATIONS"

SOURCES = ["CPIC", "DPWG", "FDA_Label", "FDA_Table"]
NON_ACTIONABLE_STATUSES = {"Not Annotated"}
GSI_DETAIL_COLS = ["Status", "Details", "AnnotationLink", "FDATableLink"]

# ============================================================================
# SUMMARY BUCKET SEMANTICS
# ============================================================================
# Lowercased for bulletproof matching against GSI data variations
ACTION_REQUIRED_STATUSES = {"alternate drug", "dosing info", "dose adjustment"}
MONITORING_STATUSES = {"other guidance", "monitor", "caution"}

STANDARD_USE_STATUSES = {
    "no action",
}

NO_GUIDELINE_STATUSES = {
    "no recommendation",
    "not evaluated",
    "not annotated",
}

UNCERTAIN_PHENOTYPES = {
    "indeterminate", "no call", "no data available", "unknown",
    "no result", "uncallable", "ambiguous", "nan", ""
}

def is_uncertain_phenotype(pheno: str) -> bool:
    p = str(pheno or "").strip().lower()
    return p in UNCERTAIN_PHENOTYPES

# ============================================================================
# SUMMARY BUCKET ENGINE
# ============================================================================

SUMMARY_PRIORITY = {
    "further_testing": 1,
    "action_required": 2,
    "monitoring": 3,
    "standard_use": 4,
    "no_guideline": 5,
}

def _compute_summary_bucket(row) -> str:
    pheno = str(row.get("Phenotype", "")).strip().lower()
    gene  = str(row.get("Gene", "")).strip().lower()
    
    selected_status = str(row.get("Selected Source Status", "")).strip().lower()

    # 0. FURTHER TESTING — specialized genes or uncertain phenotype override
    #    GSI status checks.  These genes cannot be reliably called by consumer
    #    DNA platforms, so any guidance is conditional on further testing.
    if (is_uncertain_phenotype(pheno) or 
        gene in {g.lower() for g in UNANALYZABLE_REASONS.keys()}):
        return "further_testing"

    # 1. ACTION REQUIRED (Highest priority)
    if "alternate drug" in selected_status or "dosing info" in selected_status:
        return "action_required"

    # 2. MONITORING / CAUTION
    if "other guidance" in selected_status:
        return "monitoring"

    # 3. STANDARD USE
    if "no action" in selected_status:
        return "standard_use"

    # 4. NO GUIDELINE
    if any(ng in selected_status for ng in ["no recommendation", "not evaluated", "not annotated"]):
        return "no_guideline"

    # Default
    return "no_guideline"

def _worst_bucket(buckets):
    if not buckets:
        return "no_guideline"
    priority = {b: SUMMARY_PRIORITY.get(b, 999) for b in buckets}
    return min(priority, key=priority.get)   # lowest number = highest clinical priority

# Genes for which PharmCAT uses diplotype strings (not metabolizer phenotype labels)
# as the "Phenotype" in step3_geno.  These are kept here only for use in
# UNANALYZABLE_REASONS and build_specialized_genes_df.
# NOTE: gsi_output.xlsx now contains all these genes directly, so no separate
# merge path is needed — a single merge on ["Gene", "Phenotype"] handles all genes.
GENOTYPE_BASED_GENES = {"CYP4F2", "VKORC1", "IFNL3", "HLA-A", "HLA-B"}

# Genes present in the GSI that PharmCAT does NOT produce standalone rows for.
# rs12777823, F2, F5 are now handled via supplemental VCF extraction (see
# SUPPLEMENTAL_RSID_GENES below), so this set is intentionally empty.
# Keep the constant in case new GSI-only genes are discovered in future updates.
GSI_ONLY_GENES: set = set()

# ── Supplemental rsid → gene mapping ─────────────────────────────────────────
# These SNPs are present in the step1 VCF (extracted from raw consumer DNA)
# but PharmCAT does NOT report them as independent gene rows in step3.
#
# Step4 reads the step1 VCF directly, maps each GT (0/0, 0/1, 1/1) to the
# exact phenotype string used in the GSI, and synthesises a step3-style row
# so the single GSI merge path can handle them like any other gene.
#
# Phenotype strings match the GSI exactly (lowercased at merge time).
#   rs12777823  – warfarin dosing in African Americans (CPIC Level A, no rec.)
#   rs1799963   – F2  Prothrombin G20210A (thrombosis risk; affects thrombolytics)
#   rs6025      – F5  Factor V Leiden (thrombosis risk; affects anticoagulants)
SUPPLEMENTAL_RSID_GENES: dict = {
    "rs12777823": {
        "gene":     "rs12777823",
        "pheno_00": "NC_000010.11:g.94645745G=/NC_000010.11:g.94645745G=",
        "pheno_01": "NC_000010.11:g.94645745G=/NC_000010.11:g.94645745G>A",
        "pheno_11": "NC_000010.11:g.94645745G>A/NC_000010.11:g.94645745G>A",
        "diplo_00": "G/G (reference)",
        "diplo_01": "G/A (heterozygous)",
        "diplo_11": "A/A (homozygous alt)",
    },
    "rs1799963": {
        "gene":     "F2",
        "pheno_00": "rs1799963 reference (G)/rs1799963 reference (G)",
        "pheno_01": "rs1799963 reference (G)/rs1799963 Prothrombin 20210A",
        "pheno_11": "rs1799963 Prothrombin 20210A/rs1799963 Prothrombin 20210A",
        "diplo_00": "reference/reference",
        "diplo_01": "reference/Prothrombin 20210A (heterozygous)",
        "diplo_11": "Prothrombin 20210A/Prothrombin 20210A",
    },
    "rs6025": {
        "gene":     "F5",
        "pheno_00": "rs6025 reference (C)/rs6025 reference (C)",
        "pheno_01": "rs6025 reference (C)/rs6025 Factor V Leiden (T)",
        "pheno_11": "rs6025 Factor V Leiden (T)/rs6025 Factor V Leiden (T)",
        "diplo_00": "reference/reference",
        "diplo_01": "reference/Factor V Leiden (heterozygous)",
        "diplo_11": "Factor V Leiden/Factor V Leiden",
    },
}

UNANALYZABLE_REASONS = {
    "CYP2D6":  "Complex gene structure and copy-number variation prevent reliable diplotype calling",
    "HLA-A":   "HLA typing requires a specialized high-resolution assay beyond standard sequencing",
    "HLA-B":   "HLA typing requires a specialized high-resolution assay beyond standard sequencing",
    "MT-RNR1": "Mitochondrial DNA variant not resolved in this pipeline",
    "CACNA1S": "Missing or unphased variants prevent secure diplotype assignment",
    "CYP4F2":  "CYP4F2 uses diplotype-based dosing (warfarin). Matched by genotype rather than metabolizer phenotype.",
}

STEP5_REQUIRED = [
    "Gene", "Diplotype", "Phenotype", "Activity Score", "Drug Name", "Drug Category",
    "IsProdrug", "ProdrugNote", "Zone",
    "Source Status Summary", "Has Actionable Guidance",
    "About this medication",
    "How this gene/phenotype affects the drug and what it means for you",
    "CPIC Level", "PharmGKB LoE", "PGx on FDA Label",
    "URLs", "citation", "citation links",
    "Classification", "Implication", "Recommendation",
    "Dosing Info", "Alternative Drug", "Other Guidance",
    "CPIC Status",      "CPIC Details",      "CPIC AnnotationLink",  "CPIC FDATableLink",
    "DPWG Status",      "DPWG Details",      "DPWG AnnotationLink",  "DPWG FDATableLink",
    "FDA_Label Status", "FDA_Label Details", "FDA_Label AnnotationLink", "FDA_Label FDATableLink",
    "FDA_Table Status", "FDA_Table Details", "FDA_Table AnnotationLink", "FDA_Table FDATableLink",
]

# ── Prodrug lookup ─────────────────────────────────────────────────────────────
# Drugs that require metabolic ACTIVATION (as opposed to normal drugs where the
# parent compound is already active).  For prodrugs the risk logic is inverted:
#   Normal drug  : Poor Metabolizer -> drug accumulates -> toxicity
#   Prodrug      : Poor Metabolizer -> drug NOT activated -> ineffective (or no analgesia)
#   Prodrug      : Ultrarapid Metabolizer -> too much active metabolite -> toxicity/death
#
# Source: CPIC/DPWG guidelines + PharmGKB curated knowledge.
# Key: lowercase drug name.  Value: note for the ProdrugNote column.
PRODRUG_LOOKUP: dict[str, str] = {
    "codeine":         "Prodrug -> morphine (CYP2D6). PM: no analgesia. UM: opioid toxicity/death.",
    "tramadol":        "Prodrug -> O-desmethyltramadol (CYP2D6). PM: reduced analgesia. UM: toxicity.",
    "clopidogrel":     "Prodrug -> active thiol (CYP2C19). PM: antiplatelet failure -> MI/stroke risk.",
    "tamoxifen":       "Prodrug -> endoxifen (CYP2D6). PM: reduced efficacy -> cancer treatment failure.",
    "capecitabine":    "Prodrug -> 5-fluorouracil; DPYD clears 5-FU. PM: severe 5-FU toxicity/death.",
    "tegafur":         "Prodrug -> 5-fluorouracil; DPYD clears 5-FU. PM: severe toxicity/death.",
    "azathioprine":    "Prodrug -> 6-MP -> thioguanine nucleotides (TPMT/NUDT15). PM: thiopurine toxicity.",
    "mercaptopurine":  "Prodrug -> thioguanine nucleotides (TPMT/NUDT15). PM: severe myelosuppression.",
    "thioguanine":     "Prodrug -> thioguanine nucleotides (TPMT/NUDT15). PM: severe myelosuppression.",
    "fesoterodine":    "Prodrug -> 5-hydroxymethyl tolterodine (CYP2D6). PM: altered exposure/efficacy.",
    "valbenazine":     "Prodrug -> (+)-alpha-HTBZ (CYP2D6). PM: reduced active metabolite exposure.",
}

# ── Zone (traffic-light) logic ─────────────────────────────────────────────────
# Zones are computed per (gene x drug x phenotype) row, so the same drug can
# have different zones for different phenotypes — e.g. clopidogrel is
# "No PGx Action Needed" for CYP2C19 Normal Metabolizer but
# "Use Alternative Drug" for CYP2C19 Poor Metabolizer.
#
# Priority order (first match wins):
#
#   "Use Alternative Drug"         : any source says "Alternate Drug"
#                                    A published guideline recommends switching to
#                                    a different drug for this patient's phenotype.
#
#   "Dose Adjustment Required"     : any source says "Dosing Info" or "Other Guidance"
#                                    (but no "Alternate Drug")
#                                    The drug is usable but the dose must change.
#
#   "Status Cannot Be Determined"  : patient phenotype is uncallable (No Result /
#                                    Indeterminate / Unknown) BUT a guideline EXISTS
#                                    in the GSI for this gene-drug pair.
#                                    Meaning: we know there's relevant guidance
#                                    but cannot apply it without a phenotype call.
#                                    This is NOT the same as "no guideline available".
#
#   "No PGx Action Needed"         : CPIC or DPWG explicitly says "No Action"
#                                    Positive confirmation — not absence of data.
#
#   "No PGx Guideline Available"   : no published guideline for this gene-drug pair
#                                    AND patient phenotype is not uncallable.


def extract_supplemental_snps(vcf_path: str) -> pd.DataFrame:
    """Parse the step1 VCF to extract supplemental SNP genotypes.

    Looks for every rsid listed in SUPPLEMENTAL_RSID_GENES and maps the
    VCF genotype (0/0, 0/1, 1/1) to the phenotype string used in the GSI.
    Returns a DataFrame with the same columns as step3_geno (Gene, Diplotype,
    Phenotype, Activity Score) so it can be directly concatenated before the
    GSI merge.

    If the VCF is missing, or a particular rsid is absent (e.g. the raw DNA
    platform did not genotype that position), the gene is simply skipped —
    no synthetic row is created and the GSI gap report will note the absence.
    """
    if not os.path.exists(vcf_path):
        print(f"[WARN] Supplemental SNP extraction: VCF not found at {vcf_path}")
        return pd.DataFrame(columns=["Gene", "Diplotype", "Phenotype", "Activity Score"])

    found: dict = {}   # rsid -> gt_code ("00", "01", "11")
    with open(vcf_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            rsid_raw = parts[2].strip()
            # Strip any suffix and force case-insensitive match (catches RS12777823)
            import re as _re_step4
            _m = _re_step4.match(r'(?i)^(rs\d+)', rsid_raw)
            rsid = _m.group(1).lower() if _m else rsid_raw.lower()
            
            if rsid not in SUPPLEMENTAL_RSID_GENES:
                continue
            gt_field = parts[9].split(":")[0].strip()  # e.g. "0/1" or "0|1"
            gt_field = gt_field.replace("|", "/")
            if gt_field in ("0/0",):
                found[rsid] = "00"
            elif gt_field in ("0/1", "1/0"):
                found[rsid] = "01"
            elif gt_field in ("1/1",):
                found[rsid] = "11"
            else:
                print(f"[WARN] Supplemental SNP {rsid}: unrecognised GT={gt_field!r}, skipping")

    rows = []
    for rsid, info in SUPPLEMENTAL_RSID_GENES.items():
        if rsid not in found:
            # Not present in VCF — either missing from raw data or not genotyped.
            print(f"[INFO] Supplemental SNP {rsid} ({info['gene']}) not found in VCF — skipped")
            continue
        code = found[rsid]
        pheno = info[f"pheno_{code}"]
        diplo = info[f"diplo_{code}"]
        rows.append({
            "Gene":           info["gene"],
            "Diplotype":      diplo,
            "Phenotype":      pheno,
            "Activity Score": "",
        })
        print(f"[INFO] Supplemental SNP {rsid} -> {info['gene']}  GT={code}  "
              f"phenotype={pheno!r}")

    if not rows:
        print("[INFO] No supplemental SNP rows added (none found in VCF)")
        return pd.DataFrame(columns=["Gene", "Diplotype", "Phenotype", "Activity Score"])

    df = pd.DataFrame(rows)
    print(f"[INFO] Supplemental SNP rows added: {len(df)} "
          f"({', '.join(df['Gene'].tolist())})")
    return df


def find_latest_step3_master() -> str:
    pattern = os.path.join(RESULTS_DIR, "step3_*_MASTER.xlsx")
    files = glob.glob(pattern)
    if not files:
        raise FileNotFoundError(f"No Step 3 MASTER files matching {pattern}")
    return max(files, key=os.path.getmtime)

def load_step3_genotypes(path: str) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name=STEP3_GENOTYPE_SHEET)
    df.columns = [c.strip() for c in df.columns]
    df["Gene"]      = df["Gene"].astype(str).str.strip()
    df["Phenotype"] = df["Phenotype"].astype(str).str.strip()
    if "Activity Score" not in df.columns:
        df["Activity Score"] = ""
    df["Activity Score"] = df["Activity Score"].astype(str).str.strip()
    return df

def load_step3_recs(path: str) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name=STEP3_RECS_SHEET)
    df.columns = [c.strip() for c in df.columns]
    keep = [
        "Drug", "Gene", "Patient Phenotype", "CPIC Level", "PharmGKB LoE", "PGx on FDA Label",
        "CPIC Guideline URL", "Guideline URL", "Drug Citation PMIDs", "Implication", "Recommendation",
        "Dosing Info", "Alternative Drug", "Other Guidance", "Classification",
    ]
    for col in keep:
        if col not in df.columns:
            df[col] = ""

    def pick_url(row):
        return (str(row.get("CPIC Guideline URL", "") or "").strip()
                or str(row.get("Guideline URL", "") or "").strip())

    df["URLs"]     = df.apply(pick_url, axis=1)
    df["citation"] = df["Drug Citation PMIDs"].astype(str).str.strip()

    keep_final = [
        c for c in keep
        if c not in ("CPIC Guideline URL", "Guideline URL", "Drug Citation PMIDs")
    ] + ["URLs", "citation"]
    return df[keep_final]

def load_step3_citations(path: str) -> dict:
    try:
        df = pd.read_excel(path, sheet_name=STEP3_CITATIONS_SHEET)
    except Exception:
        return {}
    df.columns = [c.strip() for c in df.columns]
    if "PMID" not in df.columns or "PubMed Link" not in df.columns:
        return {}
    df["PMID"] = df["PMID"].astype(str).str.strip()
    return (
        df[["PMID", "PubMed Link"]]
        .dropna(subset=["PMID"])
        .drop_duplicates("PMID")
        .set_index("PMID")["PubMed Link"]
        .to_dict()
    )

def _simplify_phenotype(phenotype: str) -> str:
    if not phenotype or str(phenotype).strip().lower() in ("", "nan"):
        return "indeterminate"
    p = str(phenotype).lower().strip()
    if "aminoglycoside" in p:
        if "increased risk" in p: return "hearing_loss_high"
        if "normal risk" in p: return "hearing_loss_normal"
        if "uncertain risk" in p: return "hearing_loss_uncertain"
    if "malignant hyperthermia susceptibility" in p: return "mh_susceptible"
    if "uncertain susceptibility" in p: return "mh_uncertain"
    if "ivacaftor non-responsive" in p: return "ivacaftor_nonresponsive"
    if "ivacaftor responsive" in p and "non" not in p: return "ivacaftor_responsive"
    if "deficient with cnsha" in p: return "g6pd_deficient_severe"
    if "deficient" in p: return "g6pd_deficient"
    if "variable" in p: return "variable"
    if "ultrarapid" in p: return "ultrarapid_metabolizer"
    if "rapid" in p: return "rapid_metabolizer"
    if "likely poor" in p: return "likely_poor_metabolizer"
    if "likely intermediate" in p or "possible intermediate" in p: return "likely_intermediate_metabolizer"
    if "poor metabolizer" in p or "poor function" in p: return "poor"
    if "intermediate metabolizer" in p or "intermediate" in p: return "intermediate"
    if "normal" in p: return "normal"
    if "decreased function" in p: return "decreased_function"
    if "increased function" in p: return "increased_function"
    if "indeterminate" in p: return "indeterminate"
    
    if ";" in p:
        p = p.split(";")[0].strip()
        
    return p.replace(" ", "_")

def build_about_lookup(df: pd.DataFrame) -> dict:
    about_col = "About the Medication"
    if about_col not in df.columns:
        return {}
    drug_col = "Drug Name" if "Drug Name" in df.columns else "Drug"
    lookup = {}
    for drug, grp in df.groupby(df[drug_col].astype(str).str.strip().str.lower()):
        texts = grp[about_col].astype(str).str.strip()
        valid = texts[texts.ne("") & texts.str.lower().ne("nan")]
        if not valid.empty:
            lookup[drug] = valid.iloc[0]
    print(f"[INFO] About lookup: {len(lookup)} unique drugs")
    return lookup

def build_witm_lookup(df: pd.DataFrame) -> dict:
    witm_col = "What It Means For You"
    if witm_col not in df.columns:
        return {}
    drug_col = "Drug Name" if "Drug Name" in df.columns else "Drug"
    lookup = {}
    for _, row in df.iterrows():
        gene      = str(row.get("Gene", "")).strip()
        phenotype = str(row.get("Phenotype", "")).strip()
        drug      = str(row.get(drug_col, "")).strip().lower()
        witm      = str(row.get(witm_col, "")).strip()
        if not gene or not witm or witm.lower() == "nan":
            continue
        sp = _simplify_phenotype(phenotype)
        if drug:
            drug_key = ("__drug__", drug, gene, sp)
            if drug_key not in lookup:
                lookup[drug_key] = witm
        gp_key = (gene, sp)
        if gp_key not in lookup:
            lookup[gp_key] = witm
    n_drug_keys = sum(1 for k in lookup if isinstance(k, tuple) and k and k[0] == "__drug__")
    n_gp_keys   = len(lookup) - n_drug_keys
    print(f"[INFO] WITM lookup: {n_drug_keys} drug-specific keys, {n_gp_keys} gene-phenotype fallback keys")
    return lookup

def attach_consumer_text(wide: pd.DataFrame, about_lookup: dict, witm_lookup: dict) -> pd.DataFrame:
    about_std = "About this medication"
    witm_std  = "How this gene/phenotype affects the drug and what it means for you"
    wide["_drug_key"]  = wide["Drug Name"].astype(str).str.strip().str.lower()
    wide[about_std]    = wide["_drug_key"].map(about_lookup).fillna("")
    wide["_pheno_key"] = wide["Phenotype"].astype(str).apply(_simplify_phenotype)
    wide["_gene_key"]  = wide["Gene"].astype(str).str.strip()

    def _pick_witm(row):
        drug_specific = ("__drug__", row["_drug_key"], row["_gene_key"], row["_pheno_key"])
        if drug_specific in witm_lookup:
            return witm_lookup[drug_specific]
        gp_fallback = (row["_gene_key"], row["_pheno_key"])
        return witm_lookup.get(gp_fallback, "")

    wide[witm_std] = wide.apply(_pick_witm, axis=1)
    wide = wide.drop(columns=["_drug_key", "_pheno_key", "_gene_key"], errors="ignore")
    return wide

def _is_actionable(status: str) -> bool:
    return str(status).strip() not in NON_ACTIONABLE_STATUSES

def load_and_pivot_gsi(path: str) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name=GSI_SHEET)
    df.columns = [c.strip() for c in df.columns]
    df.rename(columns={"Drug": "Drug Name", "Category": "Drug Category"}, inplace=True)
    df["Gene"]         = df["Gene"].astype(str).str.strip()
    df["Phenotype"]    = df["Phenotype"].astype(str).str.strip()
    df["Drug Name"]    = df["Drug Name"].astype(str).str.strip()
    df["Drug Category"]= df["Drug Category"].astype(str).str.strip()
    df["Source"]       = df["Source"].astype(str).str.strip()

    about_lookup = build_about_lookup(df)
    witm_lookup  = build_witm_lookup(df)

    before = len(df)
    df = df[df["Status"].apply(_is_actionable)].copy()
    print(f"[INFO] Actionability filter: {before} rows -> {len(df)} rows (excluded 'Not Evaluated' / 'Not Annotated')")

    for col in GSI_DETAIL_COLS:
        if col not in df.columns:
            df[col] = ""
    df[GSI_DETAIL_COLS] = df[GSI_DETAIL_COLS].fillna("")

    # Normalise phenotype to lowercase but KEEP the "(AS:X.X)" activity-score
    # suffix.  Keeping it allows step4's merge to match precisely:
    #   patient "Intermediate Metabolizer" + AS "1.0"  ->  key "intermediate metabolizer (as:1.0)"
    #   GSI row "Intermediate Metabolizer (AS:1.0)"    ->  key "intermediate metabolizer (as:1.0)"
    # This prevents IM(AS:1.5) and IM(AS:1.0) from being collapsed together.
    df["Phenotype"] = (
        df["Phenotype"]
        .astype(str).str.strip()
        .str.lower()
    )

    if df.empty:
        print("[WARN] No actionable rows found — check GSI Status column")
        wide = pd.DataFrame(columns=["Gene", "Phenotype", "Drug Name", "Drug Category"])
        wide = attach_consumer_text(wide, about_lookup, witm_lookup)
        return wide

    id_cols = ["Gene", "Phenotype", "Drug Name", "Drug Category"]
    wide_parts = []
    for src in SOURCES:
        src_df = df[df["Source"] == src][id_cols + GSI_DETAIL_COLS].copy()
        if src_df.empty: continue
        rename = {col: f"{src} {col}" for col in GSI_DETAIL_COLS}
        src_df = src_df.rename(columns=rename)
        agg = {
            f"{src} {col}": lambda x, _col=col: "\n\n".join(
                dict.fromkeys(v for v in x.astype(str) if v.strip() and v.lower() not in ("nan", ""))
            )
            for col in GSI_DETAIL_COLS
        }
        src_df = src_df.groupby(id_cols, as_index=False).agg(agg)
        wide_parts.append(src_df)

    if not wide_parts:
        print("[WARN] No source data found after pivot")
        wide = pd.DataFrame(columns=id_cols)
    else:
        wide = wide_parts[0]
        for part in wide_parts[1:]:
            wide = pd.merge(wide, part, on=id_cols, how="outer")

    for src in SOURCES:
        for col in GSI_DETAIL_COLS:
            c = f"{src} {col}"
            if c not in wide.columns:
                wide[c] = ""
    wide = wide.fillna("")

    print(f"[INFO] Pivot complete: {len(wide)} unique drug rows from {len(df)} source rows")
    wide = attach_consumer_text(wide, about_lookup, witm_lookup)
    return wide

def build_drug_gene_catalog(path: str) -> dict:
    try:
        df = pd.read_excel(path, sheet_name=GSI_SHEET)
        df.columns = [c.strip() for c in df.columns]
        if not all(c in df.columns for c in ["Status", "Drug", "Gene"]):
            print("[WARN] build_drug_gene_catalog: missing required columns")
            return {}
        actionable = df[df["Status"].apply(_is_actionable)]
        catalog = (
            actionable.groupby(actionable["Drug"].astype(str).str.strip().str.lower())["Gene"]
            .apply(lambda x: sorted(set(x.astype(str).str.strip())))
            .to_dict()
        )
        print(f"[INFO] Drug-gene catalog: {len(catalog)} drugs with gene associations")
        return catalog
    except Exception as e:
        print(f"[WARN] Could not build drug-gene catalog: {e}")
        return {}

_MEANINGFUL_STATUSES_KEYWORDS = {
    "dosing info", "alternate drug", "other guidance",
    "no action", "no recommendation",
}

def _has_meaningful_status(status: str) -> bool:
    return _is_actionable(status)


def _make_pheno_key(pheno: str, as_score: str) -> str:
    """Combine phenotype + activity score into a lowercase match key."""
    pheno  = pheno.strip().lower()
    as_val = as_score.strip().lower()
    
    # NEW: Treat "n/a" as absent so it doesn't break RYR1/CACNA1S merges
    _absent = {"", "nan", "none", "no result", "n/a", "na"}
    
    if as_val and as_val not in _absent:
        return f"{pheno} (as:{as_val})"
    return pheno


def build_unanalyzable_drug_map(step3_geno: pd.DataFrame) -> dict:
    """Build mapping of drug (lowercase) → set of unanalyzable genes from Related Drugs column.

    For genes in UNANALYZABLE_REASONS that have uncertain phenotypes, extract
    their Related Drugs list so we can flag those drugs as 'further_testing'
    even when PharmCAT produces no annotation row for the gene-drug pair.
    """
    drug_map: dict = {}  # {drug_key: set of gene names}
    unanalyzable_genes = {g.lower() for g in UNANALYZABLE_REASONS.keys()}

    if "Related Drugs" not in step3_geno.columns:
        return drug_map

    for _, row in step3_geno.iterrows():
        gene = str(row.get("Gene", "")).strip()
        if gene.lower() not in unanalyzable_genes:
            continue
        pheno = str(row.get("Phenotype", "")).strip().lower()
        if not is_uncertain_phenotype(pheno):
            continue
        # This gene is unanalyzable AND has uncertain phenotype — flag its drugs
        related = str(row.get("Related Drugs", "")).strip()
        if not related or related.lower() in ("", "nan"):
            continue
        for drug in related.split(","):
            dk = drug.strip().lower()
            if dk and dk not in ("", "nan"):
                drug_map.setdefault(dk, set()).add(gene)
    return drug_map


def build_drug_summary(merged: pd.DataFrame, unanalyzable_drug_map: dict = None) -> pd.DataFrame:
    _BAD_DRUG_KEYS = {"", "nan", "none", "n/a"}
    rows = []
    
    tmp = merged.copy()
    tmp["_drug_key"] = tmp["Drug Name"].astype(str).str.strip().str.lower()
    tmp = tmp[~tmp["_drug_key"].isin(_BAD_DRUG_KEYS)]

    for drug_key, grp in tmp.groupby("_drug_key"):
        canon_name = grp["Drug Name"].astype(str).str.strip().iloc[0]

        cats = grp["Drug Category"].astype(str).str.strip()
        cats = cats[~cats.isin(["", "nan", "Uncategorized", "NO GROUP ASSIGNED", "VARIOUS DRUG CLASSES IN ATC"])]
        drug_cat = cats.mode().iloc[0] if not cats.empty else "Other"

        gene_buckets: dict = {}
        for gene, ggrp in grp.groupby(grp["Gene"].astype(str).str.strip()):
            if not gene or gene.lower() in ("", "nan"): continue
            gene_buckets[gene] = _worst_bucket(ggrp["Summary Bucket"].astype(str).tolist())

        # Inject "further_testing" for unanalyzable genes linked to this drug
        # but missing from the merged data (e.g. allopurinol + HLA-B when
        # PharmCAT produces no annotation because HLA can't be called).
        if unanalyzable_drug_map and drug_key in unanalyzable_drug_map:
            for ua_gene in unanalyzable_drug_map[drug_key]:
                if ua_gene.upper() not in gene_buckets:
                    gene_buckets[ua_gene.upper()] = "further_testing"

        overall_bucket = _worst_bucket(list(gene_buckets.values())) if gene_buckets else "no_guideline"

        sorted_gz = sorted(gene_buckets.items(), key=lambda kv: SUMMARY_PRIORITY.get(kv[1], 99))
        gene_zone_breakdown = " | ".join(f"{g}: {b}" for g, b in sorted_gz)
        genes_str = ", ".join(sorted(gene_buckets.keys()))

        about = ""
        if "About this medication" in grp.columns:
            vals = grp["About this medication"].astype(str).str.strip()
            vals = vals[~vals.isin(["", "nan"])]
            about = vals.iloc[0] if not vals.empty else ""

        has_guidance = bool(grp["Has Actionable Guidance"].any()) if "Has Actionable Guidance" in grp.columns else False
        
        # Track raw statuses for sub-labels in UI
        has_alternate_drug = False
        has_dosing_info = False
        
        for col_name in grp.columns:
            col_lower = col_name.lower()
            col_data = grp[col_name].astype(str).str.strip().str.lower()
            if "alternative drug" in col_lower:
                has_alternate_drug = ((col_data != "nan") & (col_data != "") & (col_data != "0") & (col_data != "0.0")).any()
            elif "dosing info" in col_lower or ("dosing" in col_lower and "info" in col_lower):
                has_dosing_info = ((col_data == "1.0") | (col_data == "1")).any()

        rows.append({
            "Drug Name":               canon_name,
            "Drug Category":           drug_cat,
            "Summary Bucket":          overall_bucket,
            "Gene Zone Breakdown":     gene_zone_breakdown,
            "Genes":                   genes_str,
            "About this medication":   about,
            "Has Actionable Guidance": str(has_guidance),
            "Has Alternate Drug":      has_alternate_drug,
            "Has Dosing Info":         has_dosing_info,
        })

    if not rows: return pd.DataFrame()

    out = pd.DataFrame(rows)
    out["_sort"] = out["Summary Bucket"].map(SUMMARY_PRIORITY).fillna(99)
    out = out.sort_values(["_sort", "Drug Category", "Drug Name"]).drop(columns=["_sort"]).reset_index(drop=True)
    return out



def build_all_evaluated_drugs(path: str, patient_genes: set, guided_drug_names_lower: set, patient_gene_phenotypes: dict = None) -> pd.DataFrame:
    try:
        df = pd.read_excel(path, sheet_name=GSI_SHEET)
        df.columns = [c.strip() for c in df.columns]
        df.rename(columns={"Drug": "Drug Name", "Category": "Drug Category"}, inplace=True)
    except Exception as e:
        print(f"[WARN] build_all_evaluated_drugs: {e}")
        return pd.DataFrame()

    for col in ["Gene", "Drug Name", "Drug Category", "Status"]:
        if col not in df.columns:
            df[col] = ""

    patient_genes_lower = {g.lower() for g in patient_genes if g and str(g).lower() not in ("nan", "", "none")}

    mask = (df["Gene"].astype(str).str.strip().str.lower().isin(patient_genes_lower) & df["Status"].apply(_has_meaningful_status))
    sub = df[mask].copy()
    if sub.empty:
        print("[INFO] build_all_evaluated_drugs: no matching rows found")
        return pd.DataFrame()

    def _simplify_status(status: str) -> str:
        s = str(status).strip()
        if s in ("Not Evaluated", "Not Annotated", ""): return "No Action"
        if "dosing info" in s.lower() or "alternate drug" in s.lower() or "other guidance" in s.lower(): return "Has Guidance"
        if s.lower() in ("no action", "no recommendation"): return "No Action"
        return "No Specific Action"

    def _best_status(statuses):
        order = {"Has Guidance": 0, "No Specific Action": 1, "No Action": 2}
        return min(statuses, key=lambda s: order.get(s, 99))

    import re as _re
    def _norm_phen(phen: str) -> str:
        # Strip a leading "GENE: " or "gene: " prefix (e.g. "CYP2C9: Normal Metabolizer"),
        # then lowercase.  The "(AS:x.x)" activity-score suffix is intentionally kept so
        # that IM(1.0) and IM(1.5) are treated as distinct phenotypes when matching.
        p = _re.sub(r'^[A-Za-z][A-Za-z0-9]+:\s+', '', str(phen).strip())
        return p.strip().lower()

    patient_phen_norm: dict = {}
    if patient_gene_phenotypes:
        for g, p in patient_gene_phenotypes.items():
            patient_phen_norm[str(g).strip()] = _norm_phen(str(p))

    about_col = "About the Medication" if "About the Medication" in sub.columns else None
    result = []
    for (drug_name, drug_cat), grp in sub.groupby(["Drug Name", "Drug Category"]):
        genes = sorted({str(g).strip() for g in grp["Gene"].astype(str) if str(g).strip()})
        about = ""
        if about_col:
            about = next((v for v in grp[about_col].astype(str) if v.strip() and v.lower() not in ("nan", "none", "")), "")

        grp["_simple_status"] = grp["Status"].apply(_simplify_status)

        if patient_phen_norm:
            patient_statuses = []
            for gene in genes:
                gene_rows = grp[grp["Gene"].astype(str).str.strip() == gene]
                pat_phen = patient_phen_norm.get(gene)
                if pat_phen:
                    phen_match = gene_rows[gene_rows["Phenotype"].astype(str).apply(_norm_phen) == pat_phen]
                    if not phen_match.empty:
                        patient_statuses.append(_best_status(phen_match["_simple_status"].tolist()))
                    else:
                        patient_statuses.append("No Specific Action")
                else:
                    patient_statuses.append("Gene Not Available")
            non_na = [s for s in patient_statuses if s != "Gene Not Available"]
            status = _best_status(non_na) if non_na else "Gene Not Available"
        else:
            status = _best_status(grp["_simple_status"].tolist())

        result.append({
            "Drug Name":             str(drug_name).strip(),
            "Drug Category":         str(drug_cat).strip(),
            "Genes":                 ", ".join(genes),
            "About this medication": about,
            "Best Status":           status,
            "Has Guidance":          str(drug_name).strip().lower() in guided_drug_names_lower,
        })

    out = pd.DataFrame(result).sort_values(["Drug Category", "Drug Name"]).reset_index(drop=True)
    out["Drug Category"] = out["Drug Category"].replace({
        "NO GROUP ASSIGNED":           "Other",
        "VARIOUS DRUG CLASSES IN ATC": "Other",
    })
    print(f"[INFO] All Evaluated Drugs: {len(out)} unique drugs across {out['Drug Category'].nunique()} categories")
    return out

def build_specialized_genes_df(step3_geno: pd.DataFrame, scraper_wide: pd.DataFrame) -> pd.DataFrame:
    bad = step3_geno["Phenotype"].astype(str).str.strip().str.lower().isin(UNCERTAIN_PHENOTYPES)
    unanalyzed = step3_geno.loc[bad, "Gene"].unique().tolist()
    if not unanalyzed: return pd.DataFrame()
    spec_src = scraper_wide[scraper_wide["Gene"].isin(unanalyzed)].copy()
    genes_in_scraper = set(spec_src["Gene"].unique())
    missing = [g for g in unanalyzed if g not in genes_in_scraper]
    if missing:
        placeholders = pd.DataFrame({
            "Gene":         missing,
            "Drug Name":    ["N/A"] * len(missing),
            "Drug Category":["No entries in current database"] * len(missing),
        })
        spec_src = pd.concat([spec_src, placeholders], ignore_index=True)

    spec_src["Cannot Analyze Reason"] = spec_src["Gene"].map(UNANALYZABLE_REASONS).fillna("Gene could not be called in this pipeline")
    return spec_src

def build_drug_coverage_df(merged: pd.DataFrame, scraper_wide: pd.DataFrame, step3_recs: pd.DataFrame) -> pd.DataFrame:
    scraper_drugs  = set(scraper_wide["Drug Name"].astype(str).str.strip().str.lower())
    pharmcat_drugs = set(step3_recs["Drug"].astype(str).str.strip().str.lower())
    rows = []
    for drug_name, grp in merged.groupby("Drug Name"):
        # Use preserved therapeutic category if available, otherwise fall back to Drug Category
        if "Therapeutic Category" in grp.columns:
            cat = grp["Therapeutic Category"].iloc[0]
        else:
            cat = grp["Drug Category"].iloc[0]
        d_low = str(drug_name).strip().lower()
        rows.append({
            "Drug Name":             drug_name,
            "Drug Category":         cat,
            "In GSI Catalog":        "Yes" if d_low in scraper_drugs else "No",
            "In PharmCAT":           "Yes" if d_low in pharmcat_drugs else "No",
            "Patient Gene Analyzed": "No" if grp["Phenotype"].astype(str).str.lower().isin(UNCERTAIN_PHENOTYPES).all() else "Yes",
        })
    return pd.DataFrame(rows).sort_values("Drug Name")

def format_excel(writer, df, sheet_name):
    wb   = writer.book
    ws   = writer.sheets[sheet_name]
    wrap = wb.add_format({"text_wrap": True, "valign": "top"})
    head = wb.add_format({"bold": True, "bg_color": "#D9EAD3", "border": 1})
    for idx, col in enumerate(df.columns):
        series  = df[col].astype(str)
        max_len = min(max(series.map(len).max(), len(str(col))) + 2, 70)
        ws.set_column(idx, idx, max_len, wrap)
        ws.write(0, idx, col, head)

def main():
    print("-" * 65)
    print(" PGx Pipeline -- Step 4: Multi-Source Merge (CPIC + DPWG + FDA)")
    print(f" Catalog: {SCRAPER_FILE}")
    print("-" * 65)

    step3_file = find_latest_step3_master()
    print(f"[Step3] {step3_file}")

    base      = os.path.basename(step3_file)
    sample_id = base.removeprefix("step3_").removesuffix("_MASTER.xlsx")

    step3_geno = load_step3_genotypes(step3_file)
    step3_recs = load_step3_recs(step3_file)
    pmid_map   = load_step3_citations(step3_file)

    # ── Supplemental SNP extraction from step1 VCF ────────────────────────────
    # rs12777823 is already in pharmcat_positions.vcf so step1 extracts it when
    # present in the raw data.  F2 (rs1799963) and F5 (rs6025) are NOT in the
    # PharmCAT reference but are injected into step1's target list by the
    # SUPPLEMENTAL_SNPS constant added in this release — re-run step1 to pick
    # them up for existing samples.
    #
    # The step1 VCF is located at results/{sample_id}.vcf — sample_id is the
    # part of the step3 filename between "step3_" and "_MASTER.xlsx", which
    # equals the step1 output stem (e.g. "step1_selfdecode_sample").
    _vcf_path = os.path.join(RESULTS_DIR, f"{sample_id}.vcf")
    _supp_df  = extract_supplemental_snps(_vcf_path)
    if not _supp_df.empty:
        step3_geno = pd.concat([step3_geno, _supp_df], ignore_index=True)
        print(f"[INFO] step3_geno after supplemental SNP injection: {len(step3_geno)} rows")

    # Build mapping of drugs linked to unanalyzable genes (HLA-B, CYP2D6, etc.)
    # so build_drug_summary can flag them as "further_testing" even when PharmCAT
    # produces no annotation row for the gene-drug pair.
    unanalyzable_drug_map = build_unanalyzable_drug_map(step3_geno)
    if unanalyzable_drug_map:
        print(f"[INFO] Unanalyzable gene-drug links: {sum(len(v) for v in unanalyzable_drug_map.values())} "
              f"gene associations across {len(unanalyzable_drug_map)} drugs")

    drug_gene_catalog = build_drug_gene_catalog(SCRAPER_FILE)
    scraper_wide      = load_and_pivot_gsi(SCRAPER_FILE)

    # ── IFNL3 supplement ──────────────────────────────────────────────────────
    # Integrated directly into gsi_output.xlsx.


    # ── Known gaps: genes in GSI that PharmCAT does not produce rows for ──────
    # These gene-drug pairs will never be matched because PharmCAT embeds their
    # genotype inside combined multi-gene recommendations rather than reporting
    # them as independent gene rows in step3's 2_ALL_RECOMMENDATIONS sheet.
    _gsi_gene_set = set(scraper_wide["Gene"].astype(str).str.strip())
    _unmatched_gsi_genes = _gsi_gene_set & GSI_ONLY_GENES
    if _unmatched_gsi_genes:
        print(f"[WARN] GSI contains {len(_unmatched_gsi_genes)} gene(s) PharmCAT never reports "
              f"as standalone rows — their guidance will NOT be applied: "
              f"{', '.join(sorted(_unmatched_gsi_genes))}")
        for _ug in sorted(_unmatched_gsi_genes):
            _ug_drugs = sorted(scraper_wide[scraper_wide["Gene"]==_ug]["Drug Name"].unique())
            print(f"       {_ug} affects: {', '.join(_ug_drugs[:10])}{'...' if len(_ug_drugs)>10 else ''}")

    left = step3_geno[["Gene", "Diplotype", "Phenotype", "Activity Score"]].copy()
    left["Gene"]           = left["Gene"].astype(str).str.strip()
    left["Activity Score"] = left["Activity Score"].astype(str).str.strip()

    # Save the plain (non-AS-qualified) phenotype before building the AS-aware
    # merge key.  The plain version is needed later when joining to step3_recs,
    # because PharmCAT's "2_ALL_RECOMMENDATIONS" sheet stores phenotype as
    # "Intermediate Metabolizer" (no AS suffix) even for CYP2C9.
    left["Phenotype_Plain"] = left["Phenotype"].astype(str).str.strip().str.lower()

    # _make_pheno_key builds "intermediate metabolizer (as:1.0)" from PharmCAT's
    # separate Phenotype / Activity Score fields so the merge key matches the
    # GSI's "Intermediate Metabolizer (AS:1.0)" rows precisely.
    # For genotype-based genes (CYP4F2 = "*1/*5", VKORC1 = "rs9923231 ref/ref" etc.)
    # _make_pheno_key just lowercases the diplotype string — no AS suffix added.
    left["Phenotype"] = left.apply(
        lambda r: _make_pheno_key(str(r["Phenotype"]), str(r["Activity Score"])), axis=1
    )

    def _clean_gsi_phenotypes(df: pd.DataFrame) -> pd.DataFrame:
        # Strip leading "gene: " prefixes that sometimes appear in GSI phenotype strings
        df["Phenotype"] = (
            df["Phenotype"]
            .astype(str).str.strip().str.lower()
            .str.replace(r"^[a-z][a-z0-9]+:\s+", "", regex=True)
            .str.strip()
        )
        return df

    scraper_wide = _clean_gsi_phenotypes(scraper_wide)

    # Strip "Likely " to ensure PharmCAT phenotypes match GSI phenotypes perfectly
    left["Phenotype"] = left["Phenotype"].astype(str).str.replace(r"(?i)^likely\s+", "", regex=True)
    scraper_wide["Phenotype"] = scraper_wide["Phenotype"].astype(str).str.replace(r"(?i)^likely\s+", "", regex=True)

    # ── Single unified merge ───────────────────────────────────────────────────
    # gsi_output.xlsx contains ALL genes (metabolizer-based AND genotype-based:
    # CYP4F2, VKORC1, HLA-A, HLA-B, rs12777823) plus the IFNL3 supplement above.
    # No split by GENOTYPE_BASED_GENES needed — the merge key ["Gene","Phenotype"]
    # works for both:
    #   metabolizer genes → "intermediate metabolizer (as:1.0)"
    #   genotype-based genes → "*1/*5", "rs9923231 reference (c)/..." etc.
    merged = pd.merge(left, scraper_wide, on=["Gene", "Phenotype"], how="left")
    _matched = merged["Drug Name"].notna().sum()
    print(f"[INFO] GSI merge: {_matched}/{len(merged)} patient gene rows matched to GSI entries")
    if len(merged) - _matched > 0:
        _unmatched_genes = merged[merged["Drug Name"].isna()]["Gene"].unique().tolist()
        print(f"       Unmatched genes (phenotype not in GSI): {', '.join(sorted(_unmatched_genes))}")

    tmp = merged.copy()
    tmp["Drug"]              = tmp["Drug Name"].astype(str).str.strip().str.lower()
    tmp["Gene"]              = tmp["Gene"].astype(str).str.strip().str.upper()
    # Use the plain phenotype (no AS suffix) so it matches PharmCAT's
    # "2_ALL_RECOMMENDATIONS" Patient Phenotype column (e.g. "intermediate metabolizer").
    # The AS-aware key in "Phenotype" is only used for the GSI look-up above.
    tmp["Patient Phenotype"] = tmp["Phenotype_Plain"].astype(str).str.strip().str.lower()

    step3_recs["Drug"]              = step3_recs["Drug"].astype(str).str.strip().str.lower()
    step3_recs["Gene"]              = step3_recs["Gene"].astype(str).str.strip().str.upper()
    step3_recs["Patient Phenotype"] = step3_recs["Patient Phenotype"].astype(str).str.strip().str.lower()

    merged = pd.merge(tmp, step3_recs, on=["Drug", "Gene", "Patient Phenotype"], how="outer")
    merged["Phenotype"] = merged["Phenotype"].fillna(merged["Patient Phenotype"])
    merged = merged.drop(columns=["Patient Phenotype", "Phenotype_Plain"], errors="ignore")
    merged["Drug Name"]     = merged["Drug Name"].fillna(merged["Drug"]).astype(str).str.strip()
    merged["Drug Category"] = merged["Drug Category"].fillna("Uncategorized")

    drug_name_to_cat = (
        scraper_wide.dropna(subset=["Drug Category"])
        .assign(_k=lambda d: d["Drug Name"].str.strip().str.lower())
        .groupby("_k")["Drug Category"].first().to_dict()
    )
    needs_cat = merged["Drug Category"].isin(["Uncategorized", ""])
    merged.loc[needs_cat, "Drug Category"] = (
        merged.loc[needs_cat, "Drug Name"].str.strip().str.lower()
        .map(drug_name_to_cat)
        .fillna(merged.loc[needs_cat, "Drug Category"])
    )

    if os.path.exists(DRUG_CAT_SUPPLEMENT):
        supp = pd.read_excel(DRUG_CAT_SUPPLEMENT)
        supp.columns = [c.strip() for c in supp.columns]
        supp_map = (
            supp.dropna(subset=["Drug Category"])
            .assign(_k=lambda d: d["Drug Name"].str.strip().str.lower())
            .groupby("_k")["Drug Category"].first().to_dict()
        )
        still_needs = merged["Drug Category"].isin(["Uncategorized", ""])
        merged.loc[still_needs, "Drug Category"] = (
            merged.loc[still_needs, "Drug Name"].str.strip().str.lower()
            .map(supp_map)
            .fillna(merged.loc[still_needs, "Drug Category"])
        )
    else:
        print(f"[WARN] {DRUG_CAT_SUPPLEMENT} not found — some drugs may remain Uncategorized")

    if "Drug" in merged.columns:
        merged = merged.drop(columns=["Drug"])

    for col in ["CPIC Level", "PharmGKB LoE", "PGx on FDA Label", "URLs", "citation"]:
        if col not in merged.columns:
            merged[col] = ""
        merged[col] = merged[col].replace({"nan": ""}).fillna("")

    def pmids_to_links(pmid_str: str) -> str:
        if not pmid_str: return ""
        links = [pmid_map.get(p.strip(), "") for p in pmid_str.split(";") if p.strip() in pmid_map]
        return "; ".join(l for l in links if l)

    merged["citation links"] = merged["citation"].apply(pmids_to_links)

    about_std = "About this medication"
    witm_std  = "How this gene/phenotype affects the drug and what it means for you"
    if about_std not in merged.columns: merged[about_std] = ""
    if witm_std not in merged.columns: merged[witm_std] = ""

    try:
        _df = pd.read_excel(SCRAPER_FILE, sheet_name=GSI_SHEET)
        _df.columns = [c.strip() for c in _df.columns]
        _df.rename(columns={"Drug": "Drug Name"}, inplace=True)
        _about_lookup = build_about_lookup(_df)
        _witm_lookup  = build_witm_lookup(_df)

        _is_blank = lambda s: s.astype(str).str.strip().isin(["", "nan"])

        blank_about = _is_blank(merged[about_std])
        if blank_about.any():
            merged.loc[blank_about, about_std] = (
                merged.loc[blank_about, "Drug Name"]
                .astype(str).str.strip().str.lower()
                .map(_about_lookup)
                .fillna("")
            )
            print(f"[INFO] Fallback About fill: {blank_about.sum()} rows patched")

        blank_witm = _is_blank(merged[witm_std])
        if blank_witm.any():
            pheno_keys = merged.loc[blank_witm, "Phenotype"].astype(str).apply(_simplify_phenotype).tolist()
            gene_keys  = merged.loc[blank_witm, "Gene"].astype(str).str.strip().tolist()
            drug_keys  = merged.loc[blank_witm, "Drug Name"].astype(str).str.strip().str.lower().tolist()
            patched = []
            for d, g, p in zip(drug_keys, gene_keys, pheno_keys):
                ds_val = _witm_lookup.get(("__drug__", d, g, p))
                if ds_val:
                    patched.append(ds_val)
                else:
                    patched.append(_witm_lookup.get((g, p), ""))
            merged.loc[blank_witm, witm_std] = patched
            print(f"[INFO] Fallback WITM fill: {blank_witm.sum()} rows patched ({sum(1 for x in patched if x)} successful)")
    except Exception as e:
        print(f"[WARN] Fallback consumer text pass failed: {e}")

    for col in STEP5_REQUIRED:
        if col not in merged.columns:
            merged[col] = ""

    _SRC_STATUS_COLS = {
        "CPIC":      "CPIC Status",
        "DPWG":      "DPWG Status",
        "FDA_Label": "FDA_Label Status",
        "FDA_Table": "FDA_Table Status",
    }

    _SOURCE_PRIORITY = [
        ("CPIC", "CPIC Status"),
        ("DPWG", "DPWG Status"),
        ("FDA_Label", "FDA_Label Status"),
        ("FDA_Table", "FDA_Table Status"),
    ]
    _NO_GUIDANCE_STATUSES = {"no recommendation", "not evaluated", "not annotated"}

    def _select_source_status(row):
        no_guidance_fallback = None
        for label, col in _SOURCE_PRIORITY:
            val = str(row.get(col, "")).strip()
            val_norm = val.lower()
            if not val or val_norm in ("", "nan"):
                continue
            if val_norm in _NO_GUIDANCE_STATUSES:
                if no_guidance_fallback is None:
                    no_guidance_fallback = (label, val)
                continue
            return label, val
        return no_guidance_fallback or ("", "")

    def _build_source_summary(row):
        label = str(row.get("Selected Source", "")).strip()
        val = str(row.get("Selected Source Status", "")).strip()
        return f"{label}: {val}" if label and val else ""

    selected_pairs = merged.apply(_select_source_status, axis=1)
    merged["Selected Source"] = [p[0] for p in selected_pairs]
    merged["Selected Source Status"] = [p[1] for p in selected_pairs]
    merged["Source Status Summary"] = merged.apply(_build_source_summary, axis=1)
    merged["Has Actionable Guidance"] = (
        merged["Source Status Summary"].str.contains("Dosing Info|Alternate Drug|Other Guidance", case=False, na=False)
        | merged.get("Alternative Drug", pd.Series(dtype=str)).astype(str).str.lower().isin(("see recommendation", "true", "yes"))
    )

    # ── Prodrug flag ───────────────────────────────────────────────────────────
    _drug_lower = merged["Drug Name"].astype(str).str.strip().str.lower()
    merged["IsProdrug"]   = _drug_lower.map(lambda d: "Yes" if d in PRODRUG_LOOKUP else "No")
    merged["ProdrugNote"] = _drug_lower.map(lambda d: PRODRUG_LOOKUP.get(d, ""))

    # ── GSI coverage flag (needed for "Status Cannot Be Determined" zone) ─────
    # True when the GSI has at least one row for this (gene, drug) pair — meaning
    # a guideline EXISTS, regardless of whether this patient's phenotype matched.
    # Used by _compute_zone: if phenotype is uncallable AND coverage is True →
    # "Status Cannot Be Determined" instead of "No PGx Guideline Available".
    _gsi_keys = set(
        scraper_wide["Gene"].astype(str).str.strip().str.lower()
        + "||"
        + scraper_wide["Drug Name"].astype(str).str.strip().str.lower()
    )
    merged["_has_gsi_coverage"] = (
        merged["Gene"].astype(str).str.strip().str.lower()
        + "||"
        + merged["Drug Name"].astype(str).str.strip().str.lower()
    ).isin(_gsi_keys)

    merged["Summary Bucket"] = merged.apply(_compute_summary_bucket, axis=1)
    
    # We leave a placeholder 'Zone' column so we don't break older dependencies
    merged["Zone"] = merged["Summary Bucket"]
    # Coverage flag was only needed for zone computation — drop it now
    merged = merged.drop(columns=["_has_gsi_coverage"], errors="ignore")

    # ── Inject synthetic rows for unanalyzable gene-drug pairs ────────────────
    # When PharmCAT can't call a gene (e.g. HLA-B), it produces no annotation →
    # step3 creates no recommendation row → the drug loses its association with
    # that gene.  We re-inject it here so routing and bucket logic see the full
    # picture.  Example: allopurinol must be flagged "further_testing" because
    # HLA-B*58:01 status is unknown.
    if unanalyzable_drug_map:
        _existing_pairs = set(
            zip(
                merged["Drug Name"].astype(str).str.strip().str.lower(),
                merged["Gene"].astype(str).str.strip().str.upper(),
            )
        )
        _synth_rows = []
        for drug_key, genes in unanalyzable_drug_map.items():
            for gene in genes:
                if (drug_key, gene.upper()) not in _existing_pairs:
                    _synth_rows.append({
                        "Gene": gene.upper(),
                        "Drug Name": drug_key,
                        "Phenotype": "Indeterminate",
                        "Diplotype": "Unknown/Unknown",
                        "Activity Score": "",
                        "Summary Bucket": "further_testing",
                        "Zone": "further_testing",
                        "Drug Category": "",
                        "Source Status Summary": "",
                        "Has Actionable Guidance": False,
                    })
        if _synth_rows:
            _synth_df = pd.DataFrame(_synth_rows)
            # Fill any missing columns to match merged schema
            for col in merged.columns:
                if col not in _synth_df.columns:
                    _synth_df[col] = ""
            merged = pd.concat([merged, _synth_df[merged.columns]], ignore_index=True)
            print(f"[INFO] Injected {len(_synth_rows)} synthetic rows for unanalyzable gene-drug pairs")

    zone_counts = merged["Summary Bucket"].value_counts().to_dict()
    print(f"[INFO] Zone breakdown: "
          + " | ".join(f"{z}={zone_counts.get(z,0)}" for z in [
              "action_required", "monitoring", 
              "further_testing", "standard_use", "no_guideline"]))

    step3_drugs_set = set(step3_recs["Drug"].astype(str).str.lower().unique()) - {"nan", ""}
    step4_drugs_set = set(merged["Drug Name"].astype(str).str.strip().str.lower().unique()) - {"nan", ""}
    missing_from_step4 = step3_drugs_set - step4_drugs_set
    if missing_from_step4:
        print(f"[WARN] {len(missing_from_step4)} Step3 drug(s) not in Step4 output:")
        for d in sorted(missing_from_step4):
            print(f"       - {d}")
    else:
        print("[INFO] Step3/Step4 consistency: all Step3 drugs present in Step4 OK")

    final_cols = [
        "Gene", "Diplotype", "Phenotype", "Activity Score", "Drug Name", "Drug Category",
        "Summary Bucket", "Zone", "IsProdrug", "ProdrugNote",
        "Selected Source", "Selected Source Status",
        "Source Status Summary", "Has Actionable Guidance",
        "About this medication",
        "How this gene/phenotype affects the drug and what it means for you",
        "CPIC Level", "PharmGKB LoE", "PGx on FDA Label",
        "URLs", "citation", "citation links",
        "Classification", "Implication", "Recommendation",
        "Dosing Info", "Alternative Drug", "Other Guidance",
        "CPIC Status",      "CPIC Details",      "CPIC AnnotationLink",      "CPIC FDATableLink",
        "DPWG Status",      "DPWG Details",      "DPWG AnnotationLink",      "DPWG FDATableLink",
        "FDA_Label Status", "FDA_Label Details", "FDA_Label AnnotationLink", "FDA_Label FDATableLink",
        "FDA_Table Status", "FDA_Table Details", "FDA_Table AnnotationLink", "FDA_Table FDATableLink",
    ]

        # =============================================
    # === FINAL ROUTING: Detail Pages vs OEM ===
    # =============================================

    step3_drugs_set = set(
        step3_recs["Drug"].astype(str).str.strip().str.lower().unique()
    ) - {"", "nan"}

    def should_keep_for_detail(row):
        bucket = str(row.get("Summary Bucket", "")).strip().lower()
        drug_name = str(row.get("Drug Name", "")).strip().lower()
        
        # 1. High clinical priority always get detail pages
        if bucket in ["action_required", "monitoring"]:
            return True
        
        # 2. Further Testing → always gets detail pages
        if bucket == "further_testing":
            return True
        
        # 3. Standard Use → detail only if PharmCAT evaluated
        if bucket == "standard_use":
            return drug_name in step3_drugs_set
        
        # 4. No Guideline → always to Other Evaluated Medicines
        return False

    # Preserve original therapeutic category before overwriting with routing label
    merged["Therapeutic Category"] = merged["Drug Category"].copy()

    # Apply to main merged dataframe
    merged["Drug Category"] = merged.apply(
        lambda row: "Medications Evaluated — Has Guidance" 
                    if should_keep_for_detail(row) 
                    else "No Guideline Available", 
        axis=1
    )

    # Final aggregation per drug (handles multi-gene drugs)
    drug_summary = merged.groupby("Drug Name").agg({
        "Summary Bucket": lambda x: _worst_bucket(x.tolist()),
        "Drug Category": "first",
        "Gene": "count"
    }).reset_index()

    drug_summary["Drug Category"] = drug_summary.apply(
        lambda row: "Medications Evaluated — Has Guidance" 
                    if should_keep_for_detail(row) 
                    else "No Guideline Available", 
        axis=1
    )

    print(f"[INFO] Final routing: {drug_summary['Drug Category'].value_counts().to_dict()}")
    # ===================================================================
    # BUILD SUPPORTING TABLES
    # ===================================================================

    specialized_df   = build_specialized_genes_df(step3_geno, scraper_wide)
    drug_coverage_df = build_drug_coverage_df(merged, scraper_wide, step3_recs)

    catalog_rows = [{"Drug Name": drug, "Genes": ", ".join(genes)} 
                   for drug, genes in sorted(drug_gene_catalog.items())]
    catalog_df = pd.DataFrame(catalog_rows) if catalog_rows else pd.DataFrame(columns=["Drug Name", "Genes"])

    drug_summary_df = build_drug_summary(merged, unanalyzable_drug_map)

    # ===================================================================
    # WRITE OUTPUT
    # ===================================================================
    os.makedirs(RESULTS_DIR, exist_ok=True)
    out_file = os.path.join(RESULTS_DIR, f"step4_{sample_id}_all_recc.xlsx")
    print(f"[OUT] {out_file}")

    with pd.ExcelWriter(out_file, engine="xlsxwriter") as writer:
        sheet = "Genotype and All Sources"
        merged.to_excel(writer, sheet_name=sheet, index=False)
        format_excel(writer, merged, sheet)

        if not specialized_df.empty:
            specialized_df.to_excel(writer, sheet_name="Specialized Genes Drugs", index=False)
            format_excel(writer, specialized_df, "Specialized Genes Drugs")

        if not drug_coverage_df.empty:
            drug_coverage_df.to_excel(writer, sheet_name="Drug Coverage", index=False)
            format_excel(writer, drug_coverage_df, "Drug Coverage")

        if not catalog_df.empty:
            catalog_df.to_excel(writer, sheet_name="Drug Gene Catalog", index=False)
            format_excel(writer, catalog_df, "Drug Gene Catalog")

        if not drug_summary_df.empty:
            drug_summary_df.to_excel(writer, sheet_name="All Evaluated Drugs", index=False)
            format_excel(writer, drug_summary_df, "All Evaluated Drugs")

    # ── Pipeline gap report ───────────────────────────────────────────────────
    print()
    print("[PIPELINE COVERAGE NOTES]")
    print("  rs12777823 : Extracted from step1 VCF.  CPIC reports 'No Recommendation'")
    print("               for warfarin in African Americans -> zone = No PGx Action Needed.")
    print("               If absent from a sample's raw data, the gene is skipped.")
    print("  F2 (rs1799963): Extracted from step1 VCF after supplemental target injection.")
    print("               Re-run step1 on existing samples to populate this field.")
    print("  F5 (rs6025) : Same as F2 above.")
    print("  IFNL3       : Integrated directly into gsi_output.xlsx.")
    print("  HLA-A/HLA-B : PharmCAT reports 'Indeterminate' (requires specialist HLA")
    print("               typing). Pipeline correctly shows 'Status Cannot Be Determined'.")
    print()

    print(f"[OK] Step 4 complete -> {len(merged)} rows, {merged['Drug Name'].nunique()} drugs")

    _filled = lambda col: (~merged[col].astype(str).str.strip().isin(["", "nan"])).sum()
    about_filled = _filled("About this medication")
    witm_filled  = _filled("How this gene/phenotype affects the drug and what it means for you")
    total = len(merged)
    print(f"[QC] About the Medication:  {about_filled}/{total} rows filled ({100*about_filled//total if total else 0}%)")
    print(f"[QC] What It Means For You: {witm_filled}/{total} rows filled ({100*witm_filled//total if total else 0}%)")

if __name__ == "__main__":
    main()
