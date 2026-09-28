#!/usr/bin/env python3
"""
HELPER/generation.py
Clean orchestrator: HTML assembly -> PDF rendering -> cover merge.
"""

import os
import shutil
from datetime import datetime


def generate_final_html(pages, patient_name, output_dir):
    """Join all page HTML strings into one complete HTML document and save it."""
    from HELPER.htmls_drug_wise import styles
    from HELPER.utils import generate_b64_fonts

    # Each item in ``pages`` is a TOP-LEVEL REPORT TOPIC, not a physical page.
    # Insert an explicit paged-media break BETWEEN topics so the next topic
    # always starts on a fresh PDF page.  The content inside a topic is still
    # free to flow across as many physical pages as it needs.
    #
    # This is intentionally done here (at final HTML assembly) rather than
    # relying on selectors such as ``body > .page:not(:first-child)``.  It
    # makes the rule universal for every template and every report.
    topic_break = (
        '<div class="report-topic-break" aria-hidden="true" '
        'style="display:block; width:100%; height:0; margin:0; padding:0; '
        'border:0; break-after:page !important; '
        'page-break-after:always !important;"></div>'
    )
    content = ("\n" + topic_break + "\n").join(pages)

    fonts   = generate_b64_fonts()
    css     = styles(fonts, patient_name=patient_name)

    css += r"""
body > .page:not(:first-child) {
    break-before: auto !important;
    page-break-before: auto !important;
}

.report-topic-break {
    display: block !important;
    width: 100% !important;
    height: 0 !important;
    min-height: 0 !important;
    margin: 0 !important;
    padding: 0 !important;
    border: 0 !important;
    break-after: page !important;
    page-break-after: always !important;
}
"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pharmacogenomics Report - {patient_name}</title>
    <style>
        {css}
    </style>
</head>
<body>
    {content}
</body>
</html>"""

    os.makedirs(output_dir, exist_ok=True)
    timestamp = int(datetime.now().timestamp())
    html_path = os.path.join(output_dir, f"report_temp_{timestamp}.html")

    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html)

    if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
        print(f"[HTML] Saved temp HTML: {html_path}")
    return html_path



def _run_pagedjs(html_path, output_dir, patient_name="Patient"):
    """Render the assembled report HTML with the local Paged.js CLI."""
    from HELPER.utils import pagedjs_generator
    return pagedjs_generator(html_path, output_dir, patient_name=patient_name)


def generate_report(pages, patient_name, front_cover, back_cover,
                    output_folder, temp_folder, ghostscript_bin):
    """
    Main orchestrator called from step5.

    Flow:
        1. Join pages[] into full HTML document -> save to temp/
        2. Render HTML -> PDF via Paged.js (headless Chromium)
        3. Merge front cover + content + back cover via Ghostscript
        4. Copy HTML to output folder for debugging
        5. Return final PDF path

    Args:
        pages           : list of logical HTML report sections
        patient_name    : str
        front_cover     : path to front cover PDF
        back_cover      : path to back cover PDF
        output_folder   : final PDF + HTML copy destination
        temp_folder     : temp HTML storage
        ghostscript_bin : path to Ghostscript binary

    Returns:
        str — path to final merged PDF
    """
    from HELPER.utils import ghost

    # ── Step 1: HTML assembly ────────────────────────────────
    html_path = generate_final_html(pages, patient_name, temp_folder)

    # ── Step 2: HTML -> PDF ───────────────────────────────────
    if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
        print("[PDF] Rendering HTML -> PDF via Paged.js...")
    content_pdf = _run_pagedjs(html_path, temp_folder, patient_name=patient_name)
    if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
        print(f"[PDF] Content PDF: {content_pdf}")

    # ── Step 3: Merge covers ─────────────────────────────────
    os.makedirs(output_folder, exist_ok=True)
    timestamp = int(datetime.now().timestamp())
    final_pdf = os.path.join(output_folder, f"report_final_{timestamp}.pdf")

    if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
        print("[MERGE] Merging covers with Ghostscript...")
    success = ghost(
        fp      = front_cover,
        bp      = back_cover,
        content = content_pdf,
        output  = final_pdf,
        bin     = ghostscript_bin,
    )

    if not success:
        # Fallback: use content PDF as final if merge fails
        if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
            print("[WARN] Ghostscript merge failed — using content PDF as output.")
        shutil.copy(content_pdf, final_pdf)

    # ── Step 4: Copy HTML to output folder ───────────────────
    html_copy = os.path.join(output_folder, f"report_content_{timestamp}.html")
    shutil.copy(html_path, html_copy)
    if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
        print(f"[HTML] HTML copy saved: {html_copy}")

    if os.path.isdir(temp_folder):
        shutil.rmtree(temp_folder, ignore_errors=True)
        if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
            print(f"[CLEANUP] Removed temp folder: {temp_folder}")

    if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
        print(f"[DONE] Final PDF: {final_pdf}")
    return final_pdf, html_copy