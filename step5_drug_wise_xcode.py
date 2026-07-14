#!/usr/bin/env python3
import os
import glob
import sys
import re
import pandas as pd

from HELPER.htmls_drug_wise import (
    styles, welcome_template, how_to_read_template,
    faqs_template, doctor_page_template,
    drug_detail_template, genotype_summary_template,
    pharmcat_no_guidance_template,
    toc_template, specialized_genes_template,
)

from HELPER.generation import generate_report

ROOT_DIR        = os.path.abspath(os.getcwd())
RESULTS_DIR     = os.path.join(ROOT_DIR, "results")
DEPS_DIR        = os.path.join(ROOT_DIR, "02_deps")
OUTPUT_DIR      = os.path.join(RESULTS_DIR, "reports_drugwise_pdf")
TEMP_DIR        = os.path.join(RESULTS_DIR, "temp")

FRONT_COVER     = os.path.join(DEPS_DIR, "covers", "pgx_front.pdf")
BACK_COVER      = os.path.join(DEPS_DIR, "covers", "pgx_back.pdf")

def _find_ghostscript() -> str:
    if sys.platform != "win32":
        return "gs"
    import glob as _glob
    for pattern in [
        r"C:\Program Files\gs\gs*\bin\gswin64c.exe",
        r"C:\Program Files (x86)\gs\gs*\bin\gswin64c.exe",
    ]:
        matches = sorted(_glob.glob(pattern), reverse=True)
        if matches:
            return matches[0]
    return os.path.join(DEPS_DIR, "gswin64c.exe")

GHOSTSCRIPT_BIN = _find_ghostscript()
TOC_PAGE_BUDGET = 860

def _toc_entry_height(kind: str) -> int:
    """Approximate printed height for each TOC row type."""
    if kind == "chapter":
        return 36
    if kind in {"category", "numbered_category"}:
        return 32
    return 26

def _split_toc_entries(toc_entries):
    chunks = []
    current = []
    used = 0

    for entry in toc_entries:
        kind = entry[0]
        height = _toc_entry_height(kind)
        starts_new_group = kind == "chapter"

        if current and (used + height > TOC_PAGE_BUDGET or starts_new_group and used > TOC_PAGE_BUDGET * 0.85):
            chunks.append(current)
            current = []
            used = 0

        current.append(entry)
        used += height

    if current:
        chunks.append(current)

    return chunks

# Drugs that should be forced to "No Guideline Available" regardless of available data
FORCED_NO_GUIDELINE_DRUGS = {
    'irinotecan',
    'quetiapine',
    'tegafur',
}

SUMMARY_BUCKET_PRIORITY = {
    "further_testing": 1,
    "action_required": 2,
    "monitoring": 3,
    "standard_use": 4,
    "no_guideline": 5,
}

def _normalize_summary_bucket(value: str) -> str:
    bucket = str(value or "").strip().lower().replace(" / ", "_")
    bucket = bucket.replace(" ", "_").replace("-", "_")
    if bucket in SUMMARY_BUCKET_PRIORITY:
        return bucket
    if bucket in {"use_with_caution_monitoring", "use_with_caution", "monitoring"}:
        return "monitoring"
    if bucket in {"action_required", "actionrequired"}:
        return "action_required"
    if bucket in {"further_testing", "furthertesting", "further_testing_required"}:
        return "further_testing"
    if bucket in {"standard_use", "standarduse"}:
        return "standard_use"
    if bucket in {"no_guideline", "noguideline"}:
        return "no_guideline"
    return "no_guideline"

def _worst_summary_bucket(values) -> str:
    buckets = [_normalize_summary_bucket(v) for v in values if str(v).strip()]
    if not buckets:
        return "no_guideline"
    return min(buckets, key=lambda b: SUMMARY_BUCKET_PRIORITY.get(b, 999))

def _evaluation_group_for_buckets(values) -> str:
    buckets = {_normalize_summary_bucket(v) for v in values if str(v).strip()}
    if not buckets or buckets.issubset({"no_guideline"}):
        return "No Current Guidance"
    return "Has Guidance"

def _load_therapeutic_category_map() -> dict:
    for candidate in ("gsi_output.xlsx", os.path.join(ROOT_DIR, "gsi_output.xlsx")):
        if not os.path.exists(candidate):
            continue
        try:
            gsi = pd.read_excel(candidate)
            gsi.columns = [c.strip() for c in gsi.columns]
            drug_col = next((c for c in gsi.columns if str(c).strip().lower() in {"drug", "drug name"}), None)
            cat_col = next((c for c in gsi.columns if str(c).strip().lower() in {"category", "drug category"}), None)
            if not drug_col or not cat_col:
                continue
            mapping = {}
            for row in gsi.to_dict("records"):
                drug = str(row.get(drug_col, "")).strip().lower()
                category = str(row.get(cat_col, "")).strip()
                if drug and category and category.lower() not in {"nan", "none"}:
                    mapping.setdefault(drug, category)
            if mapping:
                return mapping
        except Exception:
            continue
    return {}

def _load_therapeutic_category_order() -> list:
    for candidate in ("gsi_output.xlsx", os.path.join(ROOT_DIR, "gsi_output.xlsx")):
        if not os.path.exists(candidate):
            continue
        try:
            gsi = pd.read_excel(candidate)
            gsi.columns = [c.strip() for c in gsi.columns]
            cat_col = next((c for c in gsi.columns if str(c).strip().lower() in {"category", "drug category"}), None)
            if not cat_col:
                continue
            order = [c for c in dict.fromkeys(gsi[cat_col].dropna().astype(str).tolist()) if c and c.lower() not in {"nan", "none"}]
            if order:
                return order
        except Exception:
            continue
    return []

NON_STANDARD_CATEGORIES = {
    "NO GROUP ASSIGNED",
    "VARIOUS DRUG CLASSES IN ATC",
}

def find_latest_step4():
    for pattern in [f"{RESULTS_DIR}/step4_*_all_recc.xlsx",
                    f"{RESULTS_DIR}/step4_*_all_recc.csv"]:
        files = glob.glob(pattern)
        if files:
            return max(files, key=os.path.getmtime)
    raise FileNotFoundError("No Step 4 output found in results/")

def assign_summary_bucket(row):
    clf = str(row.get("Classification", "")).strip().lower()
    alt = str(row.get("Alternate Drug", "")).strip().lower()
    dos = str(row.get("Dosing Info", "")).strip().lower()
    oth = str(row.get("Other Guidance", "")).strip().lower()
    pheno = str(row.get("Phenotype", "")).strip().lower()
    full_status = str(row.get("Selected Source Status", "")).strip().lower()
    if not full_status:
        statuses = [
            str(row.get("CPIC Status", "")).strip().lower(),
            str(row.get("DPWG Status", "")).strip().lower(),
            str(row.get("FDA_Label Status", "")).strip().lower(),
            str(row.get("FDA_Table Status", "")).strip().lower(),
        ]
        full_status = " | ".join(statuses)
    
    # Priority 1: Action Required
    if alt in ["yes", "true", "y", "see recommendation"] or dos in ["yes", "true", "y", "1", "1.0"] or "alternate drug" in full_status or "dosing info" in full_status:
        return "action_required"
    # Priority 2: Monitoring
    if oth in ["yes", "true", "y"] or "other guidance" in full_status or "monitoring" in clf or "use with caution" in clf:
        return "monitoring"
    # Priority 3: Standard Use ("No Action")
    if clf == "no action" or "no action" in full_status:
        return "standard_use"
    # Priority 4: No Guideline Available
    if clf in ["no recommendation", "not evaluated", "not annotated"]:
        return "no_guideline"
    # Priority 5: Further Testing
    # Do not treat PharmCAT's RYR1 "Uncertain Susceptibility" phenotype as a
    # no-call result. It is a valid reportable phenotype with guidance.
    if any(term in pheno for term in ["indeterminate", "no call", "no result", "unknown", "uncallable"]) or "specialized" in pheno or "further testing" in clf:
        return "further_testing"
    
    # Fallback to existing or no_guideline
    existing = str(row.get("Summary Bucket", "")).strip()
    if existing:
        return _normalize_summary_bucket(existing)
    return "no_guideline"

def load_step4_data(path):
    df = pd.read_csv(path) if path.endswith(".csv") else pd.read_excel(path)
    df.columns = [c.strip() for c in df.columns]
    df = df.dropna(subset=["Drug Name"])
    df = df[df["Drug Name"].astype(str).str.strip() != ""]

    therapeutic_map = _load_therapeutic_category_map()
    # Preserve per-row Therapeutic Category from Step 4 (supports multi-category
    # drugs like amikacin appearing in ANTIINFECTIVES, DERMATOLOGICALS, etc.).
    # Only fill in blank/NaN values using the GSI mapping as fallback.
    if "Therapeutic Category" not in df.columns:
        df["Therapeutic Category"] = ""
    _tc_missing = (
        df["Therapeutic Category"].isna()
        | df["Therapeutic Category"].astype(str).str.strip().isin(["", "nan", "None"])
    )
    if _tc_missing.any():
        df.loc[_tc_missing, "Therapeutic Category"] = (
            df.loc[_tc_missing, "Drug Name"].astype(str).str.strip().str.lower()
            .map(therapeutic_map)
        )
    if "Drug Category" not in df.columns:
        df["Drug Category"] = ""
    # Final fallback for any still-missing Therapeutic Category
    _still_missing = (
        df["Therapeutic Category"].isna()
        | df["Therapeutic Category"].astype(str).str.strip().isin(["", "nan", "None"])
    )
    if _still_missing.any():
        df.loc[_still_missing, "Therapeutic Category"] = df.loc[_still_missing, "Drug Category"]

    df["Summary Bucket"] = df.apply(assign_summary_bucket, axis=1)

    if "Evaluation Group" not in df.columns:
        df["Evaluation Group"] = ""

    group_map = df.groupby("Drug Name")["Summary Bucket"].agg(_evaluation_group_for_buckets)
    df["Evaluation Group"] = df["Drug Name"].map(group_map)

    detailed = df[df["Evaluation Group"] == "Has Guidance"]["Drug Name"].nunique()
    oem = df[df["Evaluation Group"] == "No Current Guidance"]["Drug Name"].nunique()

    # print(f"[INFO] Loaded from Step 4 -> {detailed} drugs for detailed pages | {oem} drugs for Other Evaluated Medications")

    return df

def load_specialized_genes_data(step4_path: str) -> pd.DataFrame:
    try:
        df = pd.read_excel(step4_path, sheet_name="Specialized Genes Drugs")
        df.columns = [c.strip() for c in df.columns]
        return df
    except Exception:
        return pd.DataFrame()

def load_drug_gene_catalog(step4_path: str) -> dict:
    try:
        df = pd.read_excel(step4_path, sheet_name="Drug Gene Catalog")
        df.columns = [c.strip() for c in df.columns]
        catalog = {}
        for row in df.to_dict("records"):
            drug = str(row.get("Drug Name", "")).strip().lower()
            genes_str = str(row.get("Genes", "")).strip()
            if drug and genes_str and genes_str.lower() != "nan":
                genes = [g.strip() for g in genes_str.split(",") if g.strip()]
                catalog[drug] = genes
        # print(f"[INFO] Loaded drug-gene catalog: {len(catalog)} entries")
        return catalog
    except Exception:
        # print("[WARN] Drug Gene Catalog sheet not found in Step 4 output — gene status column will be limited.")
        return {}

def load_all_evaluated_drugs(step4_path: str) -> pd.DataFrame:
    try:
        sheet_name = "Genotype and All Sources"
        try:
            df = pd.read_excel(step4_path, sheet_name=sheet_name)
        except Exception:
            sheet_name = "All Evaluated Drugs"
            df = pd.read_excel(step4_path, sheet_name=sheet_name)
        df.columns = [c.strip() for c in df.columns]
        # print(f"[INFO] Loaded {sheet_name}: {len(df)} rows | {df['Drug Name'].nunique()} unique drugs across {df['Drug Category'].nunique()} categories")
        return df
    except Exception:
        # print("[WARN] 'Genotype and All Sources' sheet not found — Executive Summary will be skipped.")
        return pd.DataFrame()

def load_qc_flags(step3_path: str) -> pd.DataFrame:
    try:
        df = pd.read_excel(step3_path, sheet_name="4_QC_FLAGS")
        df.columns = [c.strip() for c in df.columns]
        if "Status" in df.columns and "Gene" not in df.columns:
            return pd.DataFrame()
        return df
    except Exception:
        return pd.DataFrame()

def load_drug_coverage(step4_path: str) -> pd.DataFrame:
    try:
        df = pd.read_excel(step4_path, sheet_name="Drug Coverage")
        df.columns = [c.strip() for c in df.columns]
        if "In PharmCAT" not in df.columns:
            return pd.DataFrame()
        pharmcat_drugs = df[df["In PharmCAT"].astype(str).str.strip().str.lower() == "yes"].copy()
        # print(f"[INFO] Drug Coverage — {len(pharmcat_drugs)} drugs with In PharmCAT=Yes")
        return pharmcat_drugs
    except Exception:
        # print("[WARN] Drug Coverage sheet not found in Step 4 output — no-guidance list skipped.")
        return pd.DataFrame()

def _parse_pharmcat_no_guidance(sample_id: str) -> set:
    """Parse PharmCAT report HTML to get the authoritative 'Drugs With No Guidance' list."""
    import re as _re
    # Find PharmCAT report HTML
    patterns = [
        os.path.join(RESULTS_DIR, "reports", f"{sample_id}.report.html"),
        os.path.join(RESULTS_DIR, f"{sample_id}.report.html"),
    ]
    # Also search by glob
    patterns += glob.glob(os.path.join(RESULTS_DIR, "reports", "*.report.html"))
    html_path = None
    for p in patterns:
        if os.path.exists(p):
            html_path = p
            break
    if not html_path:
        # print("[WARN] PharmCAT report HTML not found — cannot parse 'Drugs With No Guidance'.")
        return set()
    try:
        with open(html_path, encoding="utf-8", errors="replace") as f:
            html = f.read()
        idx = html.find("Drugs With No Guidance")
        if idx == -1:
            # print(f"[WARN] 'Drugs With No Guidance' section not found in {os.path.basename(html_path)}")
            return set()
        end_ul = html.find("</ul>", idx)
        if end_ul == -1:
            return set()
        snippet = html[idx : end_ul + 5]
        drugs = set(d.strip().lower() for d in _re.findall(r"<li>([^<]+)</li>", snippet))
        # print(f"[INFO] PharmCAT report 'Drugs With No Guidance': {len(drugs)} drug(s)")
        return drugs
    except Exception as e:
        # print(f"[WARN] Could not parse PharmCAT report HTML: {e}")
        return set()

def load_master_data(sample_id):
    step3_files = glob.glob(os.path.join(RESULTS_DIR, "step3_*MASTER.xlsx"))
    df_genes = pd.DataFrame()
    drug_gene_map = {}
    recs_df = pd.DataFrame()

    if not step3_files:
        # print("[WARN] Could not find Step 3 MASTER file! Falling back to Step 4 data.")
        return df_genes, drug_gene_map, recs_df

    matched = [f for f in step3_files if os.path.basename(f) == f"step3_{sample_id}_MASTER.xlsx"]
    step3_path = matched[0] if matched else max(step3_files, key=os.path.getctime)
    print(f"\n[INFO] --- EXTRACTING PURE DNA DATA FROM {os.path.basename(step3_path)} ---")

    try:
        raw_df = pd.read_excel(step3_path, sheet_name="3_GENOTYPE_DETAILS")
        col_mapping = {}
        for col in raw_df.columns:
            c_lower = str(col).strip().lower()
            if c_lower in ["gene", "genes", "gene symbol"]:   col_mapping[col] = "Gene"
            elif c_lower in ["genotype", "diplotype"]:         col_mapping[col] = "Diplotype"
            elif c_lower in ["phenotype", "phenotypes"]:       col_mapping[col] = "Phenotype"

        raw_df.rename(columns=col_mapping, inplace=True)
        for req_col in ["Gene", "Diplotype", "Phenotype"]:
            if req_col not in raw_df.columns:
                raw_df[req_col] = "N/A"

        # Include Activity Score so Genotype Summary can display it alongside phenotype
        _geno_cols = ["Gene", "Diplotype", "Phenotype"]
        # Use PharmCAT Phenotype (original from report) for display if available;
        # the "Phenotype" column may contain GSI merge keys for genotype-only genes.
        if "PharmCAT Phenotype" in raw_df.columns:
            _geno_cols.append("PharmCAT Phenotype")
        if "Activity Score" in raw_df.columns:
            _geno_cols.append("Activity Score")
        df_genes = raw_df[_geno_cols].copy()
        df_genes = df_genes[df_genes["Gene"].notna() & (df_genes["Gene"] != "N/A")]

        recs_df   = pd.read_excel(step3_path, sheet_name="2_ALL_RECOMMENDATIONS")
        drug_col  = next((c for c in recs_df.columns if str(c).strip().lower() == "drug"), None)
        genes_col = next((c for c in recs_df.columns if str(c).strip().lower() in ["gene", "genes"]), None)

        if drug_col and genes_col:
            for row in recs_df.to_dict("records"):
                d     = str(row[drug_col]).strip().lower()
                g_str = str(row[genes_col]).strip()
                if d and d != "nan" and g_str and g_str != "nan":
                    genes = [x.strip() for x in g_str.replace(",", ";").split(";") if x.strip()]
                    if d not in drug_gene_map:
                        drug_gene_map[d] = set()
                    drug_gene_map[d].update(genes)
    except Exception as e:
        print(f"[ERROR] Failed to read Step 3 Master file: {e}")

    # Normalize recommendation dataframe columns for downstream matching
    try:
        if not recs_df.empty:
            recs_df.columns = [c.strip() for c in recs_df.columns]
            if "Drug" in recs_df.columns:
                recs_df["_drug_lower"] = recs_df["Drug"].astype(str).str.strip().str.lower()
            if "Gene" in recs_df.columns:
                recs_df["_gene_norm"] = recs_df["Gene"].astype(str).str.strip()
            # Ensure phenotype/diplotype/activity columns exist
            for c in ["Patient Phenotype", "Patient Diplotype", "Activity Score", "Source", "Recommendation", "Classification"]:
                if c not in recs_df.columns:
                    recs_df[c] = ""
    except Exception:
        pass

    return df_genes, drug_gene_map, recs_df

def estimate_toc_pages(df):
    category_col = "Therapeutic Category" if "Therapeutic Category" in df.columns else "Drug Category"
    categories   = df[category_col].nunique()
    unique_drugs = df["Drug Name"].nunique()
    estimated_height = (
        5 * _toc_entry_height("chapter")
        + categories * _toc_entry_height("numbered_category")
        + unique_drugs * _toc_entry_height("drug")
        + 10 * _toc_entry_height("static")
    )
    return max(1, -(-estimated_height // TOC_PAGE_BUDGET))

def _estimate_toc_pages_from_entries(toc_entries):
    used = 0
    for entry in toc_entries:
        used += _toc_entry_height(entry[0])
    return max(1, -(-used // TOC_PAGE_BUDGET))

def build_toc_pages(toc_entries, patient_name, toc_start_pg):
    page = toc_template(toc_entries, patient_name, toc_start_pg)
    return [page]

def safe_append(pages_list, template_result, current_pg):
    if isinstance(template_result, tuple):
        pages_list.extend(template_result[0])
        return template_result[1]
    else:
        pages_list.append(template_result)
        return current_pg + 1

def build_report(df, patient_name, sample_id, step4_path: str = ""):
    pages = []

    if "Evaluation Group" not in df.columns:
        df["Evaluation Group"] = "Has Guidance"

    # Use Step 4 routing derived from the 5-bucket summary, while keeping the
    # therapeutic class in Drug Category for TOC grouping.
    df_detailed = df[df["Evaluation Group"] == "Has Guidance"].copy()
    df_oem = df[df["Evaluation Group"] == "No Current Guidance"].copy()

    # print(f"[INFO] Using Step 4 routing in build_report -> {len(df_detailed['Drug Name'].unique())} detailed | "
    #       f"{len(df_oem['Drug Name'].unique())} OEM")

    df_master_genes, drug_gene_map, recs_df = load_master_data(sample_id)
    if df_master_genes.empty:
        df_master_genes = df_detailed.copy()

    drug_gene_catalog = load_drug_gene_catalog(step4_path)
    all_eval_df = load_all_evaluated_drugs(step4_path)

    step3_files = glob.glob(os.path.join(RESULTS_DIR, "step3_*MASTER.xlsx"))
    matched_s3  = [f for f in step3_files if os.path.basename(f) == f"step3_{sample_id}_MASTER.xlsx"]
    step3_path  = matched_s3[0] if matched_s3 else (max(step3_files, key=os.path.getctime) if step3_files else "")
    qc_df           = load_qc_flags(step3_path) if step3_path else pd.DataFrame()
    drug_coverage_df = load_drug_coverage(step4_path) if step4_path else pd.DataFrame()

    # ── PharmCAT-based routing ──────────────────────────────────────────
    # PRIMARY SOURCE: parse PharmCAT report HTML for the authoritative
    # "Drugs With No Guidance" list.  These drugs → No Current Guidance.
    #
    # SECONDARY: among remaining PharmCAT drugs, if ALL Step4 rows have
    # CPIC Status = "No Recommendation" AND DPWG Status in {No Recommendation,
    # Not Evaluated} — treat as No Current Guidance too.
    #
    # EVERYTHING ELSE from PharmCAT → detailed drug page (via Summary Bucket
    # routing already applied by load_step4_data).
    # ────────────────────────────────────────────────────────────────────

    # 1. Get authoritative No-Guidance set from PharmCAT report HTML
    pharmcat_no_guid_from_report = _parse_pharmcat_no_guidance(sample_id)

    # 2. Build set of all PharmCAT drug names (from Drug Coverage sheet)
    all_pharmcat_drugs = set()
    if not drug_coverage_df.empty and "Drug Name" in drug_coverage_df.columns:
        all_pharmcat_drugs = set(
            drug_coverage_df["Drug Name"].dropna().astype(str).str.strip().str.lower().unique()
        )

    # 3. Among PharmCAT drugs NOT in the No-Guidance list, check Step4 statuses:
    #    if ALL rows have only "No Recommendation"/"Not Evaluated" → No Current Guidance
    #    Also: rows where the gene has an uncertain phenotype (requires specialized
    #    testing) are treated as effectively "no recommendation" for this check.
    _no_real_rec_statuses = {"no recommendation", "not evaluated", "not annotated", ""}
    _uncertain_phenotypes = {"indeterminate", "no call", "no data available", "unknown",
                            "no result", "uncallable", "ambiguous", "nan", ""}
    _step4_demoted = set()
    pharmcat_with_guidance = all_pharmcat_drugs - pharmcat_no_guid_from_report
    for drug in pharmcat_with_guidance:
        rows = df[df["Drug Name"].astype(str).str.strip().str.lower() == drug]
        if rows.empty:
            continue

        # For each gene row, check if it's either "no recommendation" or
        # "requires specialized testing" (uncertain phenotype)
        all_rows_no_guidance = True
        for row in rows.to_dict("records"):
            pheno = str(row.get("Phenotype", "")).strip().lower()
            # If gene has uncertain phenotype → requires specialized testing, skip it
            if pheno in _uncertain_phenotypes:
                continue
            # Otherwise check if all source statuses are "no recommendation"
            cpic = str(row.get("CPIC Status", "")).strip().lower() if "CPIC Status" in rows.columns else ""
            dpwg = str(row.get("DPWG Status", "")).strip().lower() if "DPWG Status" in rows.columns else ""
            fda = str(row.get("FDA_Label Status", "")).strip().lower() if "FDA_Label Status" in rows.columns else ""
            row_statuses = {cpic, dpwg, fda}
            if not row_statuses.issubset(_no_real_rec_statuses):
                all_rows_no_guidance = False
                break

        if all_rows_no_guidance:
            _step4_demoted.add(drug)

    if _step4_demoted:
        pass
        # print(f"[INFO] {len(_step4_demoted)} PharmCAT drug(s) have only 'No Recommendation' in Step4 → No Current Guidance")
    # 4. Combine: final No Current Guidance set = PharmCAT No Guidance + Step4 demoted
    final_no_guid = pharmcat_no_guid_from_report | _step4_demoted
    
    # 4a. Add drugs with forced no_guideline status
    _forced_no_guideline = set()
    if not df.empty:
        for drug in df["Drug Name"].astype(str).str.strip().str.lower().unique():
            if drug in FORCED_NO_GUIDELINE_DRUGS:
                _forced_no_guideline.add(drug)
                final_no_guid.add(drug)
    
    if _forced_no_guideline:
        pass
        # print(f"[INFO] {len(_forced_no_guideline)} drug(s) forced to No Current Guidance:")
        for d in sorted(_forced_no_guideline):
            pass
            # print(f"       - {d}")
    
    # print(f"[INFO] Total No Current Guidance: {len(final_no_guid)} drug(s) "
    #       f"({len(pharmcat_no_guid_from_report)} from PharmCAT report + {len(_step4_demoted)} from Step4 no-rec + {len(_forced_no_guideline)} forced)")

    # 5. Apply routing: force No-Current-Guidance drugs into that group,
    #    ensure they are NOT in detailed or OEM
    if final_no_guid:
        mask = df["Drug Name"].astype(str).str.strip().str.lower().isin(final_no_guid)
        if mask.any():
            df.loc[mask, "Evaluation Group"] = "No Current Guidance"
            df.loc[mask, "Drug Category"] = "No Guideline Available"
        # Refresh partitions
        df_detailed = df[df["Evaluation Group"] == "Has Guidance"].copy()
        df_oem = df[df["Evaluation Group"] == "No Current Guidance"].copy()
        # Remove No-Guidance drugs from detailed and OEM (they go to their own section)
        df_detailed = df_detailed[~df_detailed["Drug Name"].astype(str).str.strip().str.lower().isin(final_no_guid)].copy()
        df_oem = df_oem[~df_oem["Drug Name"].astype(str).str.strip().str.lower().isin(final_no_guid)].copy()
        # print(f"[INFO] After PharmCAT routing: {df_detailed['Drug Name'].nunique()} detailed | {df_oem['Drug Name'].nunique()} OEM | {len(final_no_guid)} No Current Guidance")

    # ── Post-load demotion (Rule D) ─────────────────────────────────────
    _pharmcat_drug_names: set = set(k.lower() for k in drug_gene_map.keys())
    if not drug_coverage_df.empty and "Drug Name" in drug_coverage_df.columns:
        _pharmcat_drug_names.update(
            drug_coverage_df["Drug Name"].astype(str).str.strip().str.lower()
        )

    _rule_d_demoted = []
    for _drug_name, _grp in df.groupby("Drug Name"):
        if _grp["Drug Category"].iloc[0] == "Other Evaluated Medications":
            continue
        if _drug_name.strip().lower() not in _pharmcat_drug_names:
            df.loc[df["Drug Name"] == _drug_name, "Drug Category"] = "Other Evaluated Medications"
            _rule_d_demoted.append(_drug_name)

    if _rule_d_demoted:
        pass
        # print(f"[INFO] Rule D: {len(_rule_d_demoted)} drug(s) demoted to OEM:")
        for _d in _rule_d_demoted[:20]:
            pass
            # print(f"       - {_d}")
        if len(_rule_d_demoted) > 20:
            pass
            # print(f"       ... and {len(_rule_d_demoted) - 20} more")

        _rule_d_keys = {d.strip().lower() for d in _rule_d_demoted}
        _rule_d_mask = df_detailed["Drug Name"].astype(str).str.strip().str.lower().isin(_rule_d_keys)
        if _rule_d_mask.any():
            df_oem = pd.concat([df_oem, df_detailed[_rule_d_mask]], ignore_index=True)
            df_detailed = df_detailed[~_rule_d_mask].copy()
            df_oem.loc[
                df_oem["Drug Name"].astype(str).str.strip().str.lower().isin(_rule_d_keys),
                "Drug Category",
            ] = "Other Evaluated Medications"

    # ── Authoritative-source enforcement (strict rule requested by user)
    # Only drugs that have proper recommendations from CPIC, DPWG or FDA
    # should appear as detailed pages. Others are demoted to OEM.
    _demoted_for_no_authoritative_source = []
    try:
        # Use Step3 recommendations to precisely match patient diplotype/phenotype/activity
        def _norm_pheno(s: str) -> str:
            if s is None:
                return ""
            s = str(s).lower()
            s = re.sub(r"\(.*?\)", "", s)  # remove parenthetical
            s = s.replace("metabolizer", "").replace("function", "").strip()
            return s

        def _source_tag(src: str) -> str:
            s = str(src or "").upper()
            if "CPIC" in s:
                return "CPIC"
            if "DPWG" in s:
                return "DPWG"
            if "FDA" in s:
                return "FDA"
            return s

        recs_available = (recs_df is not None and not recs_df.empty)
        # phrases that indicate a non-meaningful placeholder recommendation
        _placeholder_phrases = [
            'does not provide a recommendation', 'no recommendation', 'no guidance',
            'no prescribing information', 'does not provide a recommendation for',
            'no cpic guidance', 'no dpwg guidance', 'no dpwg guidance.'
        ]
        # Per-drug overrides: allow targeted suppression or demotion rules
        _PER_DRUG_OVERRIDES = {
            # Remove DPWG entries for Brivaracetam (keep FDA/CPIC)
            'brivaracetam': {'suppress_sources': ['DPWG']},
            # If only placeholder recs remain for these drugs, force demotion to OEM
            'allopurinol': {'force_oem_if_only_placeholder': True},
            'elagolix': {'force_oem_if_only_placeholder': True},
            # hydralazine special ordering can be handled in templates; keep placeholder for now
            'hydralazine': {},
        }

        # Placeholder detector regex (kept in sync with HELPER/htmls_drug_wise._PLACEHOLDER_REC_RE)
        _PLACEHOLDER_REC_RE = re.compile(
            r"^(?:cpic|dpwg|fda(?:\s+label|\s+table)?)\s+"
            r"(?:has|provides)\s+no\s+"
            r"(?:recommendations?|prescribing\s+information|actionable\s+guidance)"
            r"\s+for\s",
            re.IGNORECASE,
        )

        def _is_placeholder_rec(text: str) -> bool:
            t = (text or "").strip().rstrip(".").lower()
            if not t:
                return False
            if _PLACEHOLDER_REC_RE.match(t):
                return True
            if t in ("no recommendation", "not annotated", "n/a"):
                return True
            return False

        if recs_available:
            # Pre-normalise recs drug name lower
            if "_drug_lower" not in recs_df.columns:
                recs_df["_drug_lower"] = recs_df["Drug"].astype(str).str.strip().str.lower()

            # Only protect further_testing drugs that ALSO have real analyzable
            # data (i.e. at least one non-further_testing bucket).  Drugs with
            # ONLY further_testing rows (e.g. Carvedilol → CYP2D6 unanalyzable)
            # have nothing meaningful to display and should be demoted to OEM.
            _further_testing_drugs = set()
            _ft_only_drugs = set()  # drugs with ONLY further_testing bucket (no detail page, but keep in exec summary)
            if not all_eval_df.empty and "Summary Bucket" in all_eval_df.columns:
                _buckets = all_eval_df.groupby(
                    all_eval_df["Drug Name"].astype(str).str.strip().str.lower()
                )["Summary Bucket"].apply(
                    lambda vals: set(v.strip().lower() for v in vals.astype(str) if v.strip())
                )
                _further_testing_drugs = set(
                    drug for drug, bkts in _buckets.items()
                    if "further_testing" in bkts and bkts - {"further_testing"}
                )
                _ft_only_drugs = set(
                    drug for drug, bkts in _buckets.items()
                    if bkts == {"further_testing"}
                )

            activity_score_matches = []
            demoted = _demoted_for_no_authoritative_source
            detailed_drugs = sorted(set(df_detailed["Drug Name"].astype(str).str.strip().str.lower().unique()))
            for dd in detailed_drugs:
                # Never demote further_testing drugs unless they lack mapped genes
                # in the Step 3 drug-gene recommendation map. Those drugs would
                # otherwise render an empty Genes Analyzed table.
                if dd in _further_testing_drugs:
                    if drug_gene_map.get(dd):
                        continue

                matched_recs = []
                # observed rows for this drug
                obs = all_eval_df[ all_eval_df["Drug Name"].astype(str).str.strip().str.lower() == dd ]
                if obs.empty:
                    demoted.append(dd)
                    continue

                rec_rows = recs_df[ recs_df["_drug_lower"] == dd ]
                if rec_rows.empty:
                    demoted.append(dd)
                    continue

                for r in rec_rows.to_dict("records"):
                    # Only consider recommendation rows that contain meaningful guidance
                    raw_rec_text = str(r.get("Recommendation", "") or "").strip()
                    rec_lower = raw_rec_text.lower()
                    has_rec_text = bool(raw_rec_text)
                    meaningful_rec = has_rec_text and not any(p in rec_lower for p in _placeholder_phrases)
                    has_tags = any([
                        bool(r.get("Dosing Info", False)),
                        bool(r.get("Alternative Drug", False)),
                        bool(r.get("Other Guidance", False)),
                    ])
                    clf = str(r.get("Classification", "") or "").strip().lower()
                    # Allow "no action" (Standard Use) as a valid authoritative classification
                    if not (meaningful_rec or has_tags or (clf and clf not in {"no recommendation", "not evaluated", "not annotated"})): 
                        continue
                    src_tag = _source_tag(r.get("Source", ""))
                    gene_field = str(r.get("Gene", "") or "").strip().lower()
                    rec_pheno = _norm_pheno(r.get("Patient Phenotype", "") or "")
                    rec_dip = str(r.get("Patient Diplotype", "") or "").strip().lower()
                    rec_as = r.get("Activity Score", None)

                    # iterate observed rows and attempt match
                    for o in obs.to_dict("records"):
                        ogene = str(o.get("Gene", "") or "").strip().lower()
                        odip = str(o.get("Diplotype", "") or "").strip().lower()
                        ophen = _norm_pheno(o.get("Phenotype", "") or "")
                        oas = o.get("Activity Score", None)

                        gene_match = (not gene_field) or (gene_field in ogene) or (ogene in gene_field)
                        if not gene_match:
                            continue

                        match_score = 0
                        matched = False
                        # Diplotype exact match preferred
                        if rec_dip and odip and rec_dip == odip:
                            matched = True
                            match_score = 4
                        # Activity Score match (numeric or string)
                        elif pd.notna(rec_as) and pd.notna(oas):
                            try:
                                if float(rec_as) == float(oas):
                                    matched = True
                                    match_score = 3
                            except Exception:
                                if str(rec_as).strip().lower() == str(oas).strip().lower():
                                    matched = True
                                    match_score = 3
                        # Phenotype exact/token match
                        elif rec_pheno and ophen:
                            if rec_pheno == ophen:
                                matched = True
                                match_score = 2
                            else:
                                # require token overlap rather than naive substring
                                rec_tokens = set([t for t in rec_pheno.split() if t])
                                o_tokens = set([t for t in ophen.split() if t])
                                if rec_tokens & o_tokens:
                                    matched = True
                                    match_score = 1

                        if matched:
                            matched_recs.append((src_tag, r, match_score))
                            # Record AS-match examples for later reporting
                            if match_score == 3:
                                try:
                                    activity_score_matches.append({
                                        'drug': dd,
                                        'gene_field': gene_field,
                                        'rec_diplotype': rec_dip,
                                        'obs_diplotype': odip,
                                        'rec_AS': rec_as,
                                        'obs_AS': oas,
                                        'source': src_tag,
                                    })
                                except Exception:
                                    pass
                            break

                if not matched_recs:
                    demoted.append(dd)
                    continue

                # Select sources per precedence
                srcs = set(s for s, _, _ in matched_recs)
                selected = set()
                if "CPIC" in srcs:
                    selected.add("CPIC")
                elif "DPWG" in srcs:
                    selected.add("DPWG")
                    if "FDA" in srcs:
                        selected.add("FDA")
                elif "FDA" in srcs:
                    selected.add("FDA")

                # Build combined recommendation text and classification
                chosen_recs = [ (s, r, sc) for s, r, sc in matched_recs if _source_tag(r.get("Source", "")) in selected]
                # Filter out placeholder-like recommendation rows using the
                # same detector as the HTML renderer so behaviour is consistent.
                non_placeholder = [ (s, r, sc) for s, r, sc in chosen_recs if not _is_placeholder_rec(str(r.get("Recommendation", "") or "")) ]
                # If all chosen records are placeholders, treat as no meaningful
                # recommendations and demote the drug to OEM. Honor per-drug
                # overrides which may force OEM when only placeholders remain.
                overrides = _PER_DRUG_OVERRIDES.get(dd, {})
                if not non_placeholder:
                    # If override specifically allows keeping placeholders, skip demotion
                    if overrides.get('force_oem_if_only_placeholder'):
                        demoted.append(dd)
                        continue
                    # Otherwise, default to demote placeholder-only drugs as well
                    demoted.append(dd)
                    continue
                chosen_recs = non_placeholder
                # order by match score desc to prefer diplotype/activity matches
                chosen_recs.sort(key=lambda x: x[2], reverse=True)
                rec_texts = []
                clfs = []
                for _, r, _ in chosen_recs:
                    rt = str(r.get("Recommendation", "") or "").strip()
                    if rt:
                        rec_texts.append(rt)
                    clfs.append(str(r.get("Classification", "") or "").strip().lower())

                # Join with single newline to reduce big gaps
                combined_rec = "\n".join(dict.fromkeys(rec_texts)) if rec_texts else ""
                # If override forces demotion when only placeholder recs remain,
                # check and demote accordingly.
                overrides = _PER_DRUG_OVERRIDES.get(dd, {})
                if overrides.get('force_oem_if_only_placeholder') and combined_rec:
                    # all chosen recs are placeholder-like?
                    only_placeholders = all(any(p in (rt or "").lower() for p in _placeholder_phrases) for rt in rec_texts)
                    if only_placeholders:
                        # treat as no meaningful recs → demote
                        demoted.append(dd)
                        continue
                # classification priority: strong > moderate > optional > others
                clf_priority = {"strong": 0, "moderate": 1, "optional": 2}
                clf_ranked = sorted(clfs, key=lambda x: clf_priority.get(x, 99))
                chosen_clf = clf_ranked[0] if clf_ranked else ""

                # Apply combined recommendation into df_detailed rows for this drug
                mask = df_detailed["Drug Name"].astype(str).str.strip().str.lower() == dd
                if combined_rec:
                    df_detailed.loc[mask, "Recommendation"] = combined_rec
                if chosen_clf:
                    df_detailed.loc[mask, "Classification"] = chosen_clf

            # Move demoted drugs to OEM
            for dd in demoted:
                rows_to_move = df_detailed[ df_detailed["Drug Name"].astype(str).str.strip().str.lower() == dd ]
                if not rows_to_move.empty:
                    df_oem = pd.concat([df_oem, rows_to_move], ignore_index=True)
                    df_detailed = df_detailed[ df_detailed["Drug Name"].astype(str).str.strip().str.lower() != dd ]
                    df_oem.loc[df_oem["Drug Name"].astype(str).str.strip().str.lower() == dd, "Drug Category"] = "Other Evaluated Medications"

            if demoted:
                pass
                # print(f"[INFO] {len(demoted)} drug(s) demoted to OEM (no matching CPIC/DPWG/FDA rec):")
                for d in demoted[:60]:
                    pass
                    # print(f"       - {d}")
            # Report some Activity Score matching examples to help debugging
            if activity_score_matches:
                pass
                # print(f"[INFO] Activity Score based matches captured: {len(activity_score_matches)} example(s). Showing up to 10:")
                for ex in activity_score_matches[:10]:
                    pass
                    # print(f"       - Drug: {ex.get('drug')} | Gene field: {ex.get('gene_field')} | source: {ex.get('source')} | rec_AS={ex.get('rec_AS')} obs_AS={ex.get('obs_AS')} | rec_dip={ex.get('rec_diplotype')} obs_dip={ex.get('obs_diplotype')}")
        else:
            # If no Step3 recs available, keep prior guided behavior (no-op)
            pass
    except Exception as _e:
        pass
        # print(f"[WARN] Authoritative-source enforcement failed: {_e}")

    # ── TOC & Page Generation ───────────────────────────────────────────
    df_detailed_for_toc = df_detailed.copy()
    toc_pg_count = estimate_toc_pages(df_detailed_for_toc)
    toc_start_pg = 1
    pg = 1 + toc_pg_count

    # print(f"[INFO] TOC estimate: {toc_pg_count} page(s) | content pages start at pg {pg}")

    # Introduction pages
    welcome_pg = pg
    pg = safe_append(pages, welcome_template(patient_name, pg), pg)

    how_to_read_pg = pg
    pg = safe_append(pages, how_to_read_template(patient_name, pg), pg)

    faqs_pg = pg
    pg = safe_append(pages, faqs_template(patient_name, pg), pg)

    doctor_pg = pg
    pg = safe_append(pages, doctor_page_template(df, df_master_genes, patient_name, pg), pg)

    # Executive Summary
    exec_summary_pg = pg
    from HELPER.htmls_drug_wise import executive_summary_template, coverage_statement_template
    # Executive summary shows detailed drugs + no-guideline drugs (all 5 buckets)
    try:
        detailed_drug_names = set(df_detailed['Drug Name'].astype(str).str.strip().str.lower().unique())
        # Executive summary is limited to drugs that have detailed pages plus
        # PharmCAT/no-recommendation drugs shown on the No Guideline page.
        # OEM-only drugs are intentionally excluded.
        all_summary_drugs = detailed_drug_names | (final_no_guid if final_no_guid else set())
        if not all_eval_df.empty:
            zone_df = all_eval_df[ all_eval_df['Drug Name'].astype(str).str.strip().str.lower().isin(all_summary_drugs) ].copy()
        else:
            zone_df = pd.DataFrame(columns=all_eval_df.columns if not all_eval_df.empty else ["Drug Name", "Summary Bucket", "Gene"])

        # Force final_no_guid drugs to bucket "no_guideline" in zone_df.
        # PharmCAT is authoritative: if it says "No Guidance" for a drug, that
        # overrides any synthetic "further_testing" rows from unanalyzable genes.
        if final_no_guid and not zone_df.empty:
            _ng_mask = zone_df["Drug Name"].astype(str).str.strip().str.lower().isin(final_no_guid)
            zone_df.loc[_ng_mask, "Summary Bucket"] = "no_guideline"

        # Inject placeholder rows for final_no_guid drugs NOT already in zone_df
        # so they appear in the "No Guideline Available" bucket count
        if final_no_guid:
            _existing_in_zone = set(zone_df["Drug Name"].astype(str).str.strip().str.lower().unique()) if not zone_df.empty else set()
            _missing_no_guid = final_no_guid - _existing_in_zone
            if _missing_no_guid:
                _placeholder_rows = pd.DataFrame({
                    "Drug Name": sorted(_missing_no_guid),
                    "Summary Bucket": "no_guideline",
                    "Gene": "",
                })
                # Fill missing columns
                for col in zone_df.columns:
                    if col not in _placeholder_rows.columns:
                        _placeholder_rows[col] = ""
                zone_df = pd.concat([zone_df, _placeholder_rows[zone_df.columns]], ignore_index=True)
    except Exception as _exc:
        # print(f"[WARN] Executive summary zone_df build failed: {_exc}")
        zone_df = all_eval_df
    pg = safe_append(pages, executive_summary_template(df, patient_name, pg, zone_df=zone_df), pg)

    # Coverage Statement
    coverage_statement_pg = pg
    category_col_name = "Therapeutic Category" if "Therapeutic Category" in df_detailed.columns else "Drug Category"
    present_categories_coverage = set(df_detailed[category_col_name].dropna().astype(str).unique())
    pg = safe_append(pages, coverage_statement_template(patient_name, pg, present_categories_coverage), pg)

    toc_entries = [
        ("chapter", "1. Introduction", None),
        ("static",  "Table of Contents", toc_start_pg),
        ("static",  "About Personalized Medicine", welcome_pg),
        ("static",  "About This Report", how_to_read_pg),
        ("static",  "How to Understand Your Results", faqs_pg),
        ("static",  "Note to Doctor", doctor_pg),
        ("chapter", "2. Summary of Medication Insights", None),
        ("static",  "Summary of Medication Insights", exec_summary_pg),
        ("static",  "Coverage Statement", coverage_statement_pg),
        ("chapter", "3. Drug-Specific Recommendations", None),
    ]

    # rs12777823 check
    _rs12777823_in_step1 = False
    _vcf_sample_id = sample_id.replace("step1_", "") if sample_id.startswith("step1_") else sample_id
    _step1_vcf = os.path.join(RESULTS_DIR, f"step1_{_vcf_sample_id}.vcf")
    if os.path.exists(_step1_vcf):
        try:
            with open(_step1_vcf, encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if not line.startswith("#") and "rs12777823" in line.lower():
                        _rs12777823_in_step1 = True
                        break
        except:
            pass
    # print(f"[INFO] rs12777823 in step1 VCF: {_rs12777823_in_step1}")

    # Drug Detail Pages
    category_col = "Therapeutic Category" if "Therapeutic Category" in df_detailed.columns else "Drug Category"
    present_categories = list(dict.fromkeys(df_detailed[category_col].astype(str)))
    category_order = _load_therapeutic_category_order()
    categories = [c for c in category_order if c in present_categories]
    categories.extend([c for c in present_categories if c not in categories])
    is_first_drug = True
    cat_counter = 1

    for category in categories:
        category_drugs = df_detailed[df_detailed[category_col] == category]["Drug Name"].astype(str)
        unique_drugs = list(dict.fromkeys(category_drugs))

        if not unique_drugs:
            continue

        cat_toc_idx = len(toc_entries)
        toc_entries.append(("numbered_category", f"3.{cat_counter}  {category}", None))
        cat_counter += 1

        cat_first_pg = None
        for drug in unique_drugs:
            # Render every drug under every category it belongs to so the
            # category label on the page always matches the TOC section.
            drug_rows = df_detailed[df_detailed["Drug Name"] == drug]
            if drug_rows.empty:
                continue

            drug_start_pg = pg
            if cat_first_pg is None:
                cat_first_pg = drug_start_pg

            res = drug_detail_template(
                drug_name=drug,
                category=category,
                drug_rows=drug_rows,
                patient_name=patient_name,
                curr_pg=pg,
                is_section_start=is_first_drug,
                coverage_categories=categories if is_first_drug else None,
                master_genes_df=df_master_genes,
                drug_gene_map=drug_gene_map,
                drug_gene_catalog=drug_gene_catalog,
                rs12777823_in_step1=_rs12777823_in_step1,
                per_drug_overrides=_PER_DRUG_OVERRIDES,
            )
            pg = safe_append(pages, res, pg)
            # Build category-aware anchor matching drug_detail_template's drug_id
            from HELPER.htmls_drug_wise import sanitize_id
            _toc_anchor = sanitize_id(drug) + "__" + sanitize_id(category)
            toc_entries.append(("drug", (drug, _toc_anchor), drug_start_pg))
            is_first_drug = False

        if cat_first_pg is not None:
            kind, lbl, _ = toc_entries[cat_toc_idx]
            toc_entries[cat_toc_idx] = (kind, lbl, cat_first_pg)

    # === NO GUIDELINE AVAILABLE (before OEM) ===
    toc_entries.append(("chapter", "4. Additional Information", None))

    if final_no_guid:
        # Build a DataFrame for the No Guideline Available section from the FULL
        # Drug Coverage sheet (not the PharmCAT-filtered subset, since many
        # No-Guidance drugs have In PharmCAT=No).
        _no_guid_df = pd.DataFrame()
        try:
            _dc_full = pd.read_excel(step4_path, sheet_name="Drug Coverage")
            _dc_full.columns = [c.strip() for c in _dc_full.columns]
            _no_guid_df = _dc_full[
                _dc_full["Drug Name"].astype(str).str.strip().str.lower().isin(final_no_guid)
            ].copy()
        except Exception:
            pass
        if _no_guid_df.empty:
            # Fallback: build minimal DataFrame from the drug names
            _no_guid_df = pd.DataFrame({"Drug Name": sorted(final_no_guid)})
            _no_guid_df["Drug Category"] = "No Guideline Available"
        if not _no_guid_df.empty:
            no_guid_start_pg = pg
            toc_entries.append(("section1", "No Guideline Available", no_guid_start_pg))
            res = pharmcat_no_guidance_template(_no_guid_df, patient_name, no_guid_start_pg)
            pg  = safe_append(pages, res, pg)
            # print(f"[INFO] No Guideline Available: {_no_guid_df['Drug Name'].nunique()} drugs at pg {no_guid_start_pg}")
        else:
            pass
            # print("[INFO] No drugs for 'No Guideline Available' section.")
    else:
        pass
        # print("[INFO] No drugs for 'No Guideline Available' section.")

    # === OTHER EVALUATED MEDICATIONS ===
    other_start_pg = pg
    toc_entries.append(("section1", "Other Evaluated Medications", other_start_pg))

    from HELPER.htmls_drug_wise import other_evaluated_medicines_template
    res = other_evaluated_medicines_template(df_oem, patient_name, pg, df_master_genes, drug_gene_map)
    pg = safe_append(pages, res, pg)
    # print(f"[INFO] Other Evaluated Medications: {df_oem['Drug Name'].nunique()} drugs at pg {other_start_pg}")

    # Appendix
    toc_entries.append(("chapter", "5. Appendix", None))

    spec_df = load_specialized_genes_data(step4_path)
    if not spec_df.empty:
        spec_start_pg = pg
        toc_entries.append(("section1", "Genes Requiring Specialized Testing", spec_start_pg))
        res = specialized_genes_template(spec_df, patient_name, spec_start_pg)
        pg  = safe_append(pages, res, pg)

    geno_start_pg = pg
    toc_entries.append(("section1", "Genotype Summary", geno_start_pg))
    res = genotype_summary_template(df_master_genes, patient_name, geno_start_pg, qc_df=qc_df)
    pg  = safe_append(pages, res, pg)

    from HELPER.htmls_drug_wise import disclaimer_template
    disclaimer_start_pg = pg
    toc_entries.append(("section1", "Disclaimer", disclaimer_start_pg))
    res = disclaimer_template(patient_name, disclaimer_start_pg)
    pg = safe_append(pages, res, pg)

    # Build TOC
    actual_toc_pg_count = _estimate_toc_pages_from_entries(toc_entries)
    if actual_toc_pg_count != toc_pg_count:
        delta = actual_toc_pg_count - toc_pg_count
        # print(f"[INFO] TOC adjusted: {actual_toc_pg_count} page(s) after height-aware calculation")
        toc_entries = [
            (kind, label, (pagenum + delta if isinstance(pagenum, int) and pagenum > toc_pg_count else pagenum))
            for kind, label, pagenum in toc_entries
        ]

    toc_pages = build_toc_pages(toc_entries, patient_name, toc_start_pg)
    for i, toc_page in enumerate(toc_pages):
        pages.insert(i, toc_page)

    # print(f"[INFO] Report total: {len(pages)} pages")
    return pages

def main(patient_name: str, sample_id: str = ""):
    name = patient_name.strip()

    print("\n" + "=" * 60)
    # print("  PGx Pipeline -- Step 5: PDF Report Generation")
    # print(f"  Patient: {name}")
    print("=" * 60)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(TEMP_DIR,   exist_ok=True)

    # If a specific sample_id is given (e.g. called from main.py after step4),
    # pin directly to that file so batch re-runs don't bleed into each other.
    if sample_id:
        pinned = os.path.join(RESULTS_DIR, f"step4_{sample_id}_all_recc.xlsx")
        path   = pinned if os.path.exists(pinned) else find_latest_step4()
        if not os.path.exists(pinned):
            pass
            # print(f"[WARN] Pinned step4 file not found ({pinned}), falling back to latest")
    else:
        path = find_latest_step4()
    # print(f"[DATA] {path}")

    basename  = os.path.basename(path)
    sample_id = basename.replace("step4_", "").replace("_all_recc.xlsx", "").replace("_all_recc.csv", "")

    df = load_step4_data(path)
    # print(f"[DATA] {len(df)} rows | {df['Drug Name'].nunique()} drugs | {df['Drug Category'].nunique()} categories")
    # print(f"[DATA] Categories: {sorted(df['Drug Category'].unique())}")

    pages = build_report(df, name, sample_id, step4_path=path)

    final_pdf = generate_report(
        pages=pages,
        patient_name=name,
        front_cover=FRONT_COVER,
        back_cover=BACK_COVER,
        output_folder=OUTPUT_DIR,
        temp_folder=TEMP_DIR,
        ghostscript_bin=GHOSTSCRIPT_BIN,
    )

    if final_pdf:
        print(f"\n[OK] Report ready: {final_pdf}")
    else:
        print("\n[FAIL] Generation failed -- check logs above.")

    return final_pdf

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Step 5: Generate PGx PDF report")
    ap.add_argument("--name", required=True, metavar="PATIENT_NAME")
    args = ap.parse_args()
    main(patient_name=args.name)
