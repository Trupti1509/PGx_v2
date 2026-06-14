# STEP 2: RUN PHARMCAT ON STEP 1 VCF
#!/usr/bin/env python3
import os
import sys
import subprocess

PHARMCAT_JAR = "pharmcat_data/pharmcat-3.2.0-all.jar"
REPORT_DIR   = "results/reports"

def die(msg: str, code: int = 1) -> None:
    print(msg, file=sys.stderr)
    sys.exit(code)

def derive_sample_id(vcf_path: str) -> str:
    base = os.path.basename(vcf_path)
    if base.endswith(".vcf.gz"):
        return base[:-7]
    if base.endswith(".vcf"):
        return base[:-4]
    return os.path.splitext(base)[0]

def run_pharmcat(vcf_path: str, sample_id: str) -> None:
    os.makedirs(REPORT_DIR, exist_ok=True)
    cmd = [
        "java", "-Xmx4g", "-jar", PHARMCAT_JAR,
        "-vcf", vcf_path,
        "-o",   REPORT_DIR,
        "-bf",  sample_id,
        "-reporterHtml",
        "-reporterJson",
    ]
    print(f"[PGx] Running PharmCAT on {vcf_path}")
    print("Command:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("[ERROR] PharmCAT failed")
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
        die("PharmCAT exited with non-zero status")

    json_report = os.path.join(REPORT_DIR, f"{sample_id}.report.json")
    html_report = os.path.join(REPORT_DIR, f"{sample_id}.report.html")

    if not os.path.exists(json_report):
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
        die(f"JSON report missing: {json_report}")

    print("[PGx] PharmCAT completed successfully")
    print(f"[PGx] JSON report: {json_report}")
    print(f"[PGx] HTML report: {html_report}")

def main() -> None:
    if len(sys.argv) != 2:
        die("Usage: python3 step2_pharmcat.py <step1_output.vcf>")
    vcf_path = sys.argv[1]
    if not os.path.exists(vcf_path):
        die(f"Input VCF not found: {vcf_path}")
    if not os.path.exists(PHARMCAT_JAR):
        die(f"PharmCAT JAR not found: {PHARMCAT_JAR}\nDownload from: https://github.com/PharmGKB/PharmCAT/releases")
    sample_id = derive_sample_id(vcf_path)
    print(f"[PGx] Sample ID: {sample_id}")
    run_pharmcat(vcf_path, sample_id)

if __name__ == "__main__":
    main()