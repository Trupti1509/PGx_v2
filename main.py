#!/usr/bin/env python3
import sys
import argparse
from pathlib import Path
import os

os.environ.setdefault("PGX_VERBOSE", "0")

import step1_convert_rawdna_to_vcf as step1
import step2_pharmcat as step2
import step2b_warfarin_input as warfarin_input
import step3_json_to_summary as step3
import step4_all_recc as step4
import step5_drug_wise_xcode as step5

def derive_sample_id(vcf_path: str) -> str:
    base = os.path.basename(vcf_path)
    if base.endswith(".vcf.gz"):
        return base[:-7]
    if base.endswith(".vcf"):
        return base[:-4]
    return os.path.splitext(base)[0]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("raw_input")
    parser.add_argument("--name", required=True, metavar="PATIENT_NAME", help='Patient name for the report. Example: --name "John Smith"')
    args = parser.parse_args()

    raw_path = Path(args.raw_input)
    if not raw_path.exists():
        print(f"Input file not found: {raw_path}")
        sys.exit(1)

    patient_name = args.name.strip()
    sample_basename = raw_path.stem
    step1_vcf_path = Path("results") / f"step1_{sample_basename}.vcf"
    verbose = os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}

    if verbose:
        print(f"[INFO] Patient name: {patient_name}")
    else:
        print(f"Generating report for: {patient_name}")

    print("[STEP 1/5] Raw DNA -> VCF")
    step1.run_step1(str(raw_path))

    if not step1_vcf_path.exists():
        print(f"[STEP 1/5] Failed: VCF not created at {step1_vcf_path}")
        sys.exit(1)
    print(f"  Done: {step1_vcf_path}")

    sample_id = derive_sample_id(str(step1_vcf_path))

    print("[STEP 2/5] PharmCAT analysis")
    step2.run_pharmcat(str(step1_vcf_path), sample_id)
    pharmcat_json_path = Path("results") / "reports" / f"{sample_id}.report.json"
    print(f"  Done: {pharmcat_json_path}")

    print("[STEP 2b/5] Warfarin input")
    warfarin_csv_path = warfarin_input.generate_warfarin_csv(
        vcf_path=str(step1_vcf_path),
        pharmcat_json_path=str(pharmcat_json_path),
        sample_id=sample_id,
    )
    print(f"  Done: {warfarin_csv_path}")

    print("[STEP 3/5] Summary generation")
    step3.main([])
    print("  Done: summary data prepared")

    print("[STEP 4/5] Recommendation merge")
    step4.main()
    print("  Done: recommendation catalog ready")

    print("[STEP 5/5] Final report generation")
    pdf_path = step5.main(patient_name=patient_name, sample_id=sample_id, show_summary=False)

    if pdf_path:
        timestamp = Path(pdf_path).stem.replace("report_final_", "")
        final_html = str(Path("results") / "reports_drugwise_pdf" / f"report_content_{timestamp}.html")
        print(f"PDF: {pdf_path}")
        print(f"HTML: {final_html}")
        print("STATUS: COMPLETE")
    else:
        print("[STEP 5/5] Failed: report generation failed.")

if __name__ == "__main__":
    main()