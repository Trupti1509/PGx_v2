#!/usr/bin/env python3
import os
import sys
import base64
import asyncio
import subprocess
import shutil
import io
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional


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
                'bottom': '48px',
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
        print("[POST-PROCESS] Internal links preserved and TOC activated successfully!")
        
    except Exception as e:
        print(f"[WARN] Failed to process links: {e}")


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
            _add_toc_links(output)
            
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Ghostscript failed: {e}")
        return False
    except FileNotFoundError:
        print(f"[ERROR] Ghostscript not found at: {bin}")
        return False