#!/usr/bin/env python3
import sys
import argparse
from pathlib import Path
import os

import step1_convert_rawdna_to_vcf as step1
import step2_pharmcat as step2
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
    # Made --name a required argument, removed the hardcoded default
    parser.add_argument("--name", required=True, metavar="PATIENT_NAME", help='Patient name for the report. Example: --name "John Smith"')
    args = parser.parse_args()

    raw_path = Path(args.raw_input)
    if not raw_path.exists():
        print(f"Input file not found: {raw_path}")
        sys.exit(1)

    patient_name = args.name.strip()
    sample_basename = raw_path.stem
    step1_vcf_path = Path("results") / f"step1_{sample_basename}.vcf"

    print(f"[INFO] Patient name: {patient_name}")

    print("\n--- Step 1 (Convert raw DNA to VCF) ---")
    step1.run_step1(str(raw_path))

    if not step1_vcf_path.exists():
        print(f" Step 1 failed: VCF not created at {step1_vcf_path}")
        sys.exit(1)

    sample_id = derive_sample_id(str(step1_vcf_path))
    print(f"[INFO] Sample ID for PharmCAT: {sample_id}")

    print("\n--- Step 2 (Run PharmCAT) ---")
    step2.run_pharmcat(str(step1_vcf_path), sample_id)

    print("\n--- Step 3 (JSON to summary Excel) ---")
    step3.main([])

    print("\n--- Step 4 (Add all sources recommendations) ---")
    step4.main()

    print("\n--- Step 5 (Final PDF report) ---")
    # Pass sample_id so step5 pins to the correct step4 file and avoids
    # picking a neighbour's output when multiple samples exist in results/.
    step5.main(patient_name=patient_name, sample_id=sample_id)

    print("\n" + "="*70)
    print(" All steps completed successfully!")
    print("Final PDF saved in: results/reports_drugwise_pdf/")
    print("="*70)

if __name__ == "__main__":
    main()