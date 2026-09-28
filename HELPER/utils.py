#!/usr/bin/env python3
import os
import sys
import base64
import asyncio
import subprocess
import shutil
import io
import re
import html as html_lib
from html.parser import HTMLParser
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional



class _PaginationQATokenParser(HTMLParser):
    """Collect medication names explicitly marked by report templates."""

    TARGETS = {
        "pdf-qa-oem-medication": "oem",
        "pdf-qa-specialized-medication": "specialized",
        "pdf-qa-no-guideline-medication": "no_guideline",
        "pdf-qa-genotype-diplotype": "genotype_diplotype",
        "pdf-qa-genotype-phenotype": "genotype_phenotype",
    }

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self._active = []
        self.tokens = {
            "oem": [],
            "specialized": [],
            "no_guideline": [],
            "genotype_diplotype": [],
            "genotype_phenotype": [],
        }

    @staticmethod
    def _class_set(attrs):
        for key, value in attrs:
            if key == "class" and value:
                return set(str(value).split())
        return set()

    def handle_starttag(self, tag, attrs):
        for item in self._active:
            item["depth"] += 1

        classes = self._class_set(attrs)
        for class_name, kind in self.TARGETS.items():
            if class_name in classes:
                self._active.append({"kind": kind, "depth": 1, "parts": []})

    def handle_startendtag(self, tag, attrs):
        self.handle_starttag(tag, attrs)
        self.handle_endtag(tag)

    def handle_data(self, data):
        for item in self._active:
            item["parts"].append(data)

    def handle_endtag(self, tag):
        finished = []
        for item in self._active:
            item["depth"] -= 1
            if item["depth"] == 0:
                finished.append(item)

        for item in finished:
            raw = html_lib.unescape("".join(item["parts"])).strip()
            if raw:
                self.tokens[item["kind"]].append(raw)
            self._active.remove(item)


def _pagination_qa_norm(value: str) -> str:
    """Normalize layout-only differences while preserving medication wording."""
    value = html_lib.unescape(str(value or ""))
    value = (
        value.replace("\u00ad", "")
             .replace("\u200b", "")
             .replace("\u2010", "-")
             .replace("\u2011", "-")
             .replace("\u2012", "-")
             .replace("\u2013", "-")
             .replace("\u2014", "-")
             .replace("\u2212", "-")
             .replace("\u2018", "'")
             .replace("\u2019", "'")
             .replace("\u201c", '"')
             .replace("\u201d", '"')
    )
    return re.sub(r"\s+", "", value).casefold()


def _pagination_qa_alnum(value: str) -> str:
    """Second-pass normalization tolerant of PDF punctuation extraction quirks."""
    value = html_lib.unescape(str(value or "")).casefold()
    return re.sub(r"[^a-z0-9]+", "", value)


def _clean_pdf_page_for_qa(text: str) -> str:
    """
    Remove repeated page-margin text before checking content across page breaks.

    A medication name can visually continue from one physical page to the next.
    PDF text extraction may insert the footer/header text between the two halves,
    so searching the raw concatenated PDF text can falsely report a missing drug.
    """
    kept = []
    for raw_line in str(text or "").splitlines():
        line = raw_line.strip()
        if not line:
            continue

        # Paged.js margin boxes commonly extract as one combined line.
        if re.match(
            r"^Pharmacogenomics\s+Report\s*-\s*.*?\s+\d+\s+Table\s+of\s+Contents$",
            line,
            flags=re.IGNORECASE,
        ):
            continue

        # Also tolerate extractors that split the three margin-box pieces.
        if line.casefold() == "table of contents":
            continue
        if re.match(r"^Pharmacogenomics\s+Report\s*-\s*.+$", line, flags=re.IGNORECASE):
            continue
        if re.fullmatch(r"\d{1,4}", line):
            continue

        kept.append(line)

    return "\n".join(kept)


def _find_actual_section_page(page_texts, heading, signature=None, start=0):
    """Locate a real body section rather than its Table-of-Contents mention."""
    heading_n = _pagination_qa_norm(heading)
    signature_n = _pagination_qa_norm(signature) if signature else None

    for idx in range(max(0, start), len(page_texts)):
        text_n = _pagination_qa_norm(_clean_pdf_page_for_qa(page_texts[idx]))
        if heading_n not in text_n:
            continue
        if signature_n and signature_n not in text_n:
            continue
        return idx
    return None


def _qa_token_present_in_section(token: str, section_texts) -> bool:
    """
    Search a complete rendered section, not a single physical page.

    This deliberately joins all physical pages belonging to the section after
    removing repeated margin-box text.  Therefore a word or medication may
    wrap across any physical page boundary without creating a false failure.
    Multiple PDF extractors may be supplied; a match from either extractor is
    accepted because text extraction order/ligatures differ between engines.
    """
    target = _pagination_qa_norm(token)
    target_alnum = _pagination_qa_alnum(token)

    if not target:
        return True

    for section_text in section_texts:
        if section_text is None:
            continue
        normal = _pagination_qa_norm(section_text)
        if target in normal:
            return True
        if target_alnum and target_alnum in _pagination_qa_alnum(section_text):
            return True

    return False


def _extract_pdf_text_candidates(pdf_path: str):
    """Return page-text arrays from independent extractors when available."""
    candidates = []

    # PyMuPDF usually preserves Chromium/Paged.js text more faithfully and is
    # the primary validator.  PyPDF2 can occasionally crash or get interrupted on
    # some generated PDFs, so it must never be allowed to kill the QA pipeline.
    try:
        try:
            import pymupdf as fitz
        except ImportError:
            import fitz
        doc = fitz.open(pdf_path)
        candidates.append(("PyMuPDF", [(page.get_text("text") or "") for page in doc]))
        doc.close()
    except BaseException as exc:
        print(f"[PDF QA WARN] PyMuPDF text extraction unavailable: {exc}")

    # Keep PyPDF2 as an independent fallback/check, but treat parser errors and
    # interrupts as non-fatal so the report can continue using the working
    # extractor instead of aborting the build.
    try:
        from PyPDF2 import PdfReader
        reader = PdfReader(pdf_path)
        candidates.append(("PyPDF2", [(page.extract_text() or "") for page in reader.pages]))
    except BaseException as exc:
        print(f"[PDF QA WARN] PyPDF2 text extraction unavailable: {exc}")

    return candidates


def _build_section_texts(page_texts):
    """Locate dynamic report sections and return their cleaned full text."""
    no_guideline_start = _find_actual_section_page(
        page_texts,
        "No Guideline Available",
        "The medications listed below are associated with genes analysed in this report",
    )
    oem_start = _find_actual_section_page(
        page_texts,
        "Other Evaluated Medications",
        "These medications were evaluated using your genetic profile",
    )
    specialized_start = _find_actual_section_page(
        page_texts,
        "Genes Requiring Specialized Testing",
        "Cannot Be Called",
    )
    genotype_start = _find_actual_section_page(
        page_texts,
        "Genotype Summary",
        start=(specialized_start + 1) if specialized_start is not None else 0,
    )
    disclaimer_start = _find_actual_section_page(
        page_texts,
        "Disclaimer",
        "Purpose and Limitations of This Report",
        start=(genotype_start + 1) if genotype_start is not None else 0,
    )

    def joined(start_idx, end_idx):
        if start_idx is None:
            return None
        if end_idx is None or end_idx <= start_idx:
            end_idx = len(page_texts)
        return "\n".join(
            _clean_pdf_page_for_qa(page_texts[i])
            for i in range(start_idx, min(end_idx, len(page_texts)))
        )

    return {
        "no_guideline": joined(no_guideline_start, oem_start),
        "oem": joined(oem_start, specialized_start),
        "specialized": joined(specialized_start, genotype_start),
        "genotype": joined(genotype_start, disclaimer_start),
    }


def validate_dynamic_pagination(html_path: str, pdf_path: str, strict: bool = True) -> bool:
    """
    Provider/sample-agnostic HTML -> PDF content-integrity check.

    It verifies the variable-length sections that are most exposed to pagination:
      * No Guideline Available medication rows
      * Other Evaluated Medications medication names
      * Genes Requiring Specialized Testing medication names
      * Genotype Summary diplotype and phenotype strings

    The check does not know any company, gene, drug count, or page number.  It
    simply compares text explicitly emitted by the HTML template with the text
    that survived into the rendered PDF.  Section text is concatenated across
    all physical pages, so harmless page wrapping cannot cause a false failure.
    """
    if not os.path.exists(html_path) or not os.path.exists(pdf_path):
        msg = "[PDF QA] HTML or PDF is missing; cannot validate dynamic pagination."
        if strict:
            raise RuntimeError(msg)
        print("[PDF QA WARN]", msg)
        return False

    try:
        with open(html_path, "r", encoding="utf-8") as f:
            html_text = f.read()

        parser = _PaginationQATokenParser()
        parser.feed(html_text)

        extractor_candidates = _extract_pdf_text_candidates(pdf_path)
        if not extractor_candidates:
            raise RuntimeError("No PDF text extractor is available for pagination QA.")

        # Build the same logical section ranges independently for each extractor.
        section_candidates = []
        for extractor_name, page_texts in extractor_candidates:
            section_candidates.append((extractor_name, _build_section_texts(page_texts)))

        checks = {
            "no_guideline": ("No Guideline Available", "no_guideline"),
            "oem": ("Other Evaluated Medications", "oem"),
            "specialized": ("Genes Requiring Specialized Testing", "specialized"),
            "genotype_diplotype": ("Genotype Summary diplotype(s)", "genotype"),
            "genotype_phenotype": ("Genotype Summary phenotype(s)", "genotype"),
        }

        failures = []
        checked = {}

        for token_kind, (label, section_key) in checks.items():
            expected = []
            seen = set()
            for token in parser.tokens.get(token_kind, []):
                key = _pagination_qa_alnum(token)
                if key and key not in seen:
                    seen.add(key)
                    expected.append(token)

            checked[token_kind] = len(expected)
            if not expected:
                continue

            available_section_texts = [
                sections.get(section_key)
                for _, sections in section_candidates
                if sections.get(section_key) is not None
            ]

            if not available_section_texts:
                failures.append(f"{label}: section could not be located in rendered PDF.")
                continue

            missing = [
                token
                for token in expected
                if not _qa_token_present_in_section(token, available_section_texts)
            ]

            if missing:
                preview = ", ".join(missing[:10])
                suffix = f" (+{len(missing)-10} more)" if len(missing) > 10 else ""
                failures.append(
                    f"{label}: {len(missing)} expected item(s) could not be found "
                    f"in the rendered section: {preview}{suffix}"
                )

        if failures:
            message = (
                "Paged.js content-integrity QA failed. HTML content was not "
                "verified in the rendered PDF:\n- " + "\n- ".join(failures)
            )
            if strict:
                raise RuntimeError(message)
            print("[PDF QA WARN]", message)
            return False

        if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
            print(
                "[PDF QA] Dynamic content verified: "
                f"NoGuideline={checked.get('no_guideline', 0)}, "
                f"OEM={checked.get('oem', 0)}, "
                f"Specialized={checked.get('specialized', 0)}, "
                f"GenotypeDiplotypes={checked.get('genotype_diplotype', 0)}, "
                f"GenotypePhenotypes={checked.get('genotype_phenotype', 0)}."
            )
        return True

    except RuntimeError:
        raise
    except Exception as exc:
        msg = f"Dynamic pagination QA could not complete: {exc}"
        if strict:
            raise RuntimeError(msg) from exc
        print("[PDF QA WARN]", msg)
        return False


def generate_b64_fonts(font_dir: Optional[str] = None) -> Dict[str, str]:
    """
    Load DM Sans font files and encode to base64 for HTML embedding.
    Searches project root and 02_deps as fallback.
    """
    if font_dir is None:
        root = os.getcwd()
        candidates = [
            os.path.join(root, "DM_Sans", "static"),
            os.path.join(root, "02_deps", "fonts"),
            os.path.join(root, "fonts"),
        ]
        font_dir = next((c for c in candidates if os.path.exists(c)), candidates[0])

    font_files = {
        'regular': 'DMSans-Regular.ttf',
        'bold':    'DMSans-Bold.ttf',
        'italic':  'DMSans-Italic.ttf',
    }

    fonts = {}
    for key, fname in font_files.items():
        path = os.path.join(font_dir, fname)
        if os.path.exists(path):
            with open(path, 'rb') as f:
                fonts[key] = base64.b64encode(f.read()).decode('utf-8')
        else:
            print(f"[WARN] Font not found: {path}")
            fonts[key] = ''

    return fonts


def pagedjs_generator(
    html_path: str,
    temp_folder: str,
    patient_name: str = "Patient",
    timeout_ms: int = 120000,
) -> str:
    """
    Render a local HTML report to PDF with the project's local Paged.js CLI.

    Expected project setup (run once from project root):
        npm install --save-dev pagedjs-cli@0.4.3

    The function intentionally uses the local node_modules installation so
    report generation is reproducible and does not depend on a global npm
    package. Paged.js owns physical pagination; Python does not estimate page
    heights or row counts.
    """
    _ = patient_name  # Footer text is embedded in CSS during HTML assembly.

    os.makedirs(temp_folder, exist_ok=True)
    pdf_path = os.path.join(
        temp_folder,
        f"report_content_{int(datetime.now().timestamp())}.pdf"
    )

    root = Path(os.getcwd())

    # Use the exact invocation already verified in this project environment:
    #     npx pagedjs-cli ...
    # Because pagedjs-cli is installed in this project's devDependencies, npx
    # resolves that local copy first.
    npx = shutil.which("npx")
    if npx:
        cli = [npx, "pagedjs-cli"]
    else:
        local_candidates = [
            root / "node_modules" / ".bin" / "pagedjs-cli",
            root / "node_modules" / ".bin" / "pagedjs-cli.cmd",
        ]
        local_bin = next((p for p in local_candidates if p.exists()), None)
        if local_bin is None:
            raise RuntimeError(
                "Paged.js is not available. Install Node/npm and run: "
                "npm install --save-dev pagedjs-cli@0.4.3"
            )
        if str(local_bin).lower().endswith(".cmd") and sys.platform == "win32":
            cli = ["cmd.exe", "/c", str(local_bin)]
        else:
            cli = [str(local_bin)]

    # IMPORTANT FOR WSL + WINDOWS NODE/NPM:
    #
    # In this setup `npx` is launched from WSL but Node ultimately resolves
    # paths on the Windows side. Passing a WSL absolute path such as
    #   /mnt/c/Users/.../report.html
    # can therefore be misread by Node as
    #   C:\\mnt\\c\\Users\\...\\report.html
    # which does not exist.
    #
    # Use paths relative to the project root instead. This is the same form
    # that was verified manually with:
    #   npx pagedjs-cli results/.../report.html -o results/.../report.pdf
    # and it works on native Linux/WSL/Windows because subprocess runs with
    # cwd set to the project root below.
    html_arg = os.path.relpath(os.path.abspath(html_path), str(root))
    pdf_arg = os.path.relpath(os.path.abspath(pdf_path), str(root))

    cmd = cli + [
        html_arg,
        "-o", pdf_arg,
        "--page-size", "A4",
        "--media", "print",
        "--timeout", str(timeout_ms),
    ]

    if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
        print("[PDF] Rendering HTML file with Paged.js...")

    result = subprocess.run(
        cmd,
        cwd=str(root),
        capture_output=True,
        text=True,
    )

    if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
        if result.stdout.strip():
            print(result.stdout.strip())
        if result.stderr.strip():
            # pagedjs-cli sometimes writes progress/warnings to stderr even when
            # successful, so only treat a non-zero return code as failure.
            print(result.stderr.strip())

    if result.returncode != 0:
        raise RuntimeError(
            "Paged.js PDF rendering failed "
            f"(exit code {result.returncode})."
        )

    if not os.path.exists(pdf_path) or os.path.getsize(pdf_path) == 0:
        raise RuntimeError(
            f"Paged.js completed but did not create a valid PDF: {pdf_path}"
        )

    if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
        print(f"[PDF] Paged.js content PDF: {pdf_path}")

    # Verify that variable-length medication content survived pagination before
    # covers are merged.  The verifier is tolerant of a medication name that
    # spans two physical pages and of footer/header text inserted by extraction.
    validate_dynamic_pagination(
        html_path=html_path,
        pdf_path=pdf_path,
        strict=True,
    )

    return pdf_path


async def pyppeteer_generator(html_path: str, temp_folder: str, patient_name: str = "Patient") -> str:
    from pyppeteer import launch

    pdf_path = os.path.join(
        temp_folder,
        f"report_content_{int(datetime.now().timestamp())}.pdf"
    )

    launch_args = {
        'headless': True,
        'handleSIGINT':  False,
        'handleSIGTERM': False,
        'handleSIGHUP':  False,
        'args': ['--no-sandbox', '--disable-setuid-sandbox', '--disable-dev-shm-usage'],
    }

    if sys.platform == 'win32':
        edge = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"
        if os.path.exists(edge):
            launch_args['executablePath'] = edge

    browser = await launch(**launch_args)

    try:
        page = await browser.newPage()
        page.setDefaultNavigationTimeout(0)

        file_url = f'file:///{os.path.abspath(html_path).replace(os.sep, "/")}'
        print("[PDF] Rendering HTML file as PDF...")
        await page.goto(file_url, {'waitUntil': 'networkidle0'})
        await page.evaluate('document.fonts.ready')
        await asyncio.sleep(0.2)

        await page.pdf({
            'path': pdf_path,
            'format': 'A4',
            'printBackground': True,
            'displayHeaderFooter': True,
            'headerTemplate': '<div></div>',
            'footerTemplate': """
                <div style="-webkit-print-color-adjust:exact; color-adjust:exact; width:100%; padding: 0 52px; font-family:Arial, sans-serif; font-size:9px; background:white; margin-bottom:5px;">
                    <div style="height:1.5px; background-color:#00B5C8; margin-bottom:4px; width:100%;"></div>
                    <div style="display:flex; justify-content:space-between; align-items:center; width:100%;">
                        <div style="flex:1; text-align:left; color:#444444;"><span class="title"></span></div>
                        <div style="flex:1; text-align:center; color:#444444; font-weight:600;"><span class="pageNumber"></span></div>
                        <div style="flex:1; text-align:right; color:#00B5C8; font-style:italic;">Table of Contents</div>
                    </div>
                </div>
            """,
            'margin': {
                'top': '30px',
                'bottom': '64px',
                'left': '52px',
                'right': '52px'
            }
        })
    finally:
        await browser.close()

    print(f"[PDF] Content PDF: {pdf_path}")
    return pdf_path


def _add_toc_links(pdf_path: str):
    try:
        from PyPDF2 import PdfReader, PdfWriter
        from PyPDF2.generic import NameObject, DictionaryObject, ArrayObject, FloatObject, NumberObject
    except ImportError:
        print("[WARN] PyPDF2 not installed. TOC links will not be clickable. Run: pip install PyPDF2")
        return

    if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
        print("[POST-PROCESS] Reactivating TOC links and preserving internal anchors...")
    
    try:
        reader = PdfReader(pdf_path)
        writer = PdfWriter()
        
        # --- CRITICAL FIX FOR INTERNAL LINKS ---
        # writer.append() safely copies all pages AND the document's hidden Root Catalog
        # (which contains the named destinations for the body HTML links).
        if hasattr(writer, 'append'):
            writer.append(reader)
        else:
            # Failsafe for older versions of PyPDF2
            for page in reader.pages:
                writer.add_page(page)
            if '/Names' in reader.trailer['/Root']:
                writer._root_object[NameObject('/Names')] = reader.trailer['/Root']['/Names']
            if '/Dests' in reader.trailer['/Root']:
                writer._root_object[NameObject('/Dests')] = reader.trailer['/Root']['/Dests']

        # 1. Smarter TOC text detection
        toc_page_idx = None
        for i, page in enumerate(writer.pages):
            text = page.extract_text()
            if text and "Table of Contents" in text and "Introduction" in text:
                toc_page_idx = i
                break

        if toc_page_idx is not None:
            toc_page_ref = writer.pages[toc_page_idx].indirect_reference

            # 2. Generous bounding box applied to ALL pages except the TOC itself
            for i, page in enumerate(writer.pages):
                if i == toc_page_idx:
                    continue

                annotation = DictionaryObject({
                    NameObject('/Type'): NameObject('/Annot'),
                    NameObject('/Subtype'): NameObject('/Link'),
                    NameObject('/Rect'): ArrayObject([
                        FloatObject(450.0), FloatObject(10.0),
                        FloatObject(560.0), FloatObject(35.0),
                    ]),
                    NameObject('/Border'): ArrayObject([NumberObject(0), NumberObject(0), NumberObject(0)]),
                    NameObject('/H'): NameObject('/N'), 
                    NameObject('/A'): DictionaryObject({
                        NameObject('/Type'): NameObject('/Action'),
                        NameObject('/S'): NameObject('/GoTo'),
                        NameObject('/D'): ArrayObject([toc_page_ref, NameObject('/Fit')]),
                    }),
                })

                # Safely inject the footer link
                if NameObject('/Annots') in page:
                    annots = page[NameObject('/Annots')]
                    if hasattr(annots, 'get_object'):
                        annots = annots.get_object()
                    annots.append(annotation)
                else:
                    page[NameObject('/Annots')] = ArrayObject([annotation])
        else:
            print("[WARN] Could not locate TOC page. Footer links not added.")

        # Write the final fixed PDF
        with open(pdf_path, "wb") as f:
            writer.write(f)
        if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
            print("[POST-PROCESS] Internal links preserved and TOC activated successfully!")
        
    except Exception as e:
        print(f"[WARN] Failed to process links: {e}")

# ============================================================================
# INTRO PAGE VISUAL STYLE
# ============================================================================

def _style_intro_pdf_pages(
    pdf_path: str,
    fill_hex: str = "#F5FAFD",
    border_hex: str = "#A9D7E5",
    border_width: float = 1.0,
    border_inset: float = 12.0,
):
    """
    Apply a true full-page tint + border to the complete introductory section.

    Styled range:
        About Personalized Medicine
        ...
        About This Report
        ...
        How to Understand Your Results
        ...
        Note to Doctor

    Styling stops immediately BEFORE the real:
        Summary of Medication Insights
    """

    try:
        try:
            import pymupdf as fitz  # current PyMuPDF API
        except ImportError:
            import fitz  # backward-compatible fallback

    except ImportError:
        print(
            "[WARN] PyMuPDF is not installed. "
            "Run: pip install pymupdf"
        )
        return


    def _hex_to_rgb(hex_color):
        h = hex_color.lstrip("#")

        return (
            int(h[0:2], 16) / 255.0,
            int(h[2:4], 16) / 255.0,
            int(h[4:6], 16) / 255.0,
        )


    border_rgb = _hex_to_rgb(
        border_hex
    )

    fill_rgb = (
        _hex_to_rgb(fill_hex)
        if fill_hex
        else None
    )


    doc = fitz.open(
        pdf_path
    )


    # ============================================================
    # PASS 1:
    # Find the ACTUAL intro start and summary start.
    # ============================================================

    intro_start_idx = None
    summary_start_idx = None


    for page_index in range(
        len(doc)
    ):

        text = (
            doc[page_index]
            .get_text("text")
            or ""
        )


        # --------------------------------------------------------
        # ACTUAL "About Personalized Medicine" page.
        #
        # Requiring the subheading prevents the Table of Contents
        # from being mistaken for the start page.
        # --------------------------------------------------------

        if (
            intro_start_idx is None
            and "About Personalized Medicine" in text
            and "What Is Personalized Medicine and Why Does It Matter?"
            in text
        ):

            intro_start_idx = (
                page_index
            )


        # --------------------------------------------------------
        # ACTUAL "Summary of Medication Insights" page.
        #
        # Note to Doctor mentions the summary by name, so simply
        # searching for the title is NOT enough.
        #
        # "Action Required" + "Standard Use" are unique to the
        # actual summary page.
        # --------------------------------------------------------

        if (
            summary_start_idx is None
            and "Summary of Medication Insights" in text
            and "Action Required" in text
            and "Standard Use" in text
        ):

            summary_start_idx = (
                page_index
            )


    # ============================================================
    # SAFETY CHECKS
    # ============================================================

    if intro_start_idx is None:

        print(
            "[WARN] Could not find the actual "
            "'About Personalized Medicine' page. "
            "Intro page styling was skipped."
        )

        doc.close()
        return


    if summary_start_idx is None:

        print(
            "[WARN] Could not find the actual "
            "'Summary of Medication Insights' page. "
            "Intro page styling was skipped."
        )

        doc.close()
        return


    if (
        summary_start_idx
        <= intro_start_idx
    ):

        print(
            "[WARN] Intro styling boundaries are invalid. "
            f"start={intro_start_idx + 1}, "
            f"summary={summary_start_idx + 1}. "
            "Styling was skipped."
        )

        doc.close()
        return


    if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
        print(
            "[PDF STYLE] Intro range detected: "
            f"PDF pages {intro_start_idx + 1} "
            f"through {summary_start_idx}"
        )


    # ============================================================
    # PASS 2:
    # Style EVERY page between intro start and summary.
    #
    # This guarantees that Note to Doctor and any continuation
    # pages are included.
    # ============================================================

    styled_pages = []


    for page_index in range(
        intro_start_idx,
        summary_start_idx
    ):

        page = doc[
            page_index
        ]

        page_rect = (
            page.rect
        )


        # --------------------------------------------------------
        # 1. FULL PHYSICAL-PAGE PALE BLUE TINT
        # --------------------------------------------------------

        if fill_rgb is not None:

            page.draw_rect(
                page_rect,

                color=None,
                fill=fill_rgb,

                # Put tint underneath existing text/tables.
                overlay=False,
            )


        # --------------------------------------------------------
        # 2. THIN FULL-PAGE BORDER
        # --------------------------------------------------------

        frame_rect = fitz.Rect(
            page_rect.x0
            + border_inset,

            page_rect.y0
            + border_inset,

            page_rect.x1
            - border_inset,

            page_rect.y1
            - border_inset,
        )


        page.draw_rect(
            frame_rect,

            color=border_rgb,

            fill=None,

            width=border_width,

            # Keep border visible above page content.
            overlay=True,
        )


        styled_pages.append(
            page_index + 1
        )


    # ============================================================
    # SAVE SAFELY
    # ============================================================

    temp_path = (
        pdf_path
        + ".intro_styled.pdf"
    )


    doc.save(
        temp_path,
        garbage=4,
        deflate=True,
    )


    doc.close()


    os.replace(
        temp_path,
        pdf_path
    )


    if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
        print(
            "[PDF STYLE] Intro pages styled: "
            + ", ".join(
                str(p)
                for p in styled_pages
            )
        )

def ghost(fp: str, bp: str, content: str, output: str, bin: str) -> bool:
    if not os.path.exists(content):
        print(f"[ERROR] Content PDF not found: {content}")
        return False

    if not os.path.exists(fp) or not os.path.exists(bp):
        print("[WARN] Cover PDFs not found - outputting content PDF only.")
        shutil.copy(content, output)
        return True

    cmd = [
        bin,
        '-q', '-dNOPAUSE', '-dBATCH',
        '-sDEVICE=pdfwrite',
        f'-sOutputFile={output}',
        fp, content, bp
    ]

    if os.environ.get("PGX_VERBOSE", "0").lower() not in {"0", "false", "no", "off"}:
        print(f"[GHOST] Merging: front + content + back -> {os.path.basename(output)}")

    try:
        if sys.platform == 'win32':
            startupinfo = subprocess.STARTUPINFO()
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            subprocess.run(cmd, startupinfo=startupinfo,
                           check=True, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
        else:
            subprocess.run(cmd, check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                           
        # Run PyPDF2 immediately after Ghostscript merge completes
        if os.path.exists(output):
            # ------------------------------------------------------------
            # Give the introductory/instruction pages their own visual
            # treatment before rebuilding the clickable TOC/footer links.
            # ------------------------------------------------------------

            _style_intro_pdf_pages(
                output,

                fill_hex="#F5FAFD",
                border_hex="#A9D7E5",
                border_width=1.0,
                border_inset=12.0,
            )

            _add_toc_links(output)

        return True

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Ghostscript failed: {e}")
        return False
    except FileNotFoundError:
        print(f"[ERROR] Ghostscript not found at: {bin}")
        return False