#!/usr/bin/env python3
"""
HELPER/htmls_drug_wise.py

"""

import re
import os
import html as html_lib
import pandas as pd

root = os.getcwd()
INTRO_ICONS_DIR = os.path.join(root, '02_deps', 'intro_icons')



def _comma_safe_chunks(items, max_chars=220):
    """
    Split a dynamic medication list into short, pagination-safe comma groups.

    This is intentionally provider/gene/drug agnostic.  It does NOT decide
    physical pages.  It only prevents Paged.js from having to fragment one
    very large inline text box across a page boundary.

    Each returned chunk is small enough to be kept together with
    ``break-inside: avoid`` while the overall list can still flow across any
    number of physical pages.
    """
    chunks = []
    current = []
    current_chars = 0

    for value in items or []:
        text = str(value or '').strip()
        if not text:
            continue

        extra = len(text) + (2 if current else 0)  # comma + space

        if current and current_chars + extra > max_chars:
            chunks.append(current)
            current = [text]
            current_chars = len(text)
        else:
            current.append(text)
            current_chars += extra

    if current:
        chunks.append(current)

    return chunks

# ============================================================================
# FDA DDI REFERENCE
# "Other Medications That May Compound This Result"
# ============================================================================

DDI_EXCEL_PATH = os.path.join(
    root,
    '02_deps',
    'Xcode_DDI_Reference.xlsx'
)

_DDI_REFERENCE_CACHE = None


def _load_ddi_reference():
    """
    Load the FDA-derived DDI workbook once.

    Returns:
        CYP dataframe
        Transporter dataframe
        gene -> (FDA target, reference sheet) mapping
    """
    global _DDI_REFERENCE_CACHE

    if _DDI_REFERENCE_CACHE is not None:
        return _DDI_REFERENCE_CACHE

    empty = (
        pd.DataFrame(),
        pd.DataFrame(),
        {}
    )

    if not os.path.exists(DDI_EXCEL_PATH):
        print(
            f"[WARN] DDI reference not found: "
            f"{DDI_EXCEL_PATH}"
        )
        _DDI_REFERENCE_CACHE = empty
        return _DDI_REFERENCE_CACHE

    try:
        cyp_df = pd.read_excel(
            DDI_EXCEL_PATH,
            sheet_name='CYP_Reference'
        )

        transporter_df = pd.read_excel(
            DDI_EXCEL_PATH,
            sheet_name='Transporter_Reference'
        )

        mapping_df = pd.read_excel(
            DDI_EXCEL_PATH,
            sheet_name='Gene_Mapping'
        )

        gene_map = {}

        for row in mapping_df.to_dict('records'):

            gene = str(
                row.get(
                    'Gene in Xcode Life Report',
                    ''
                )
            ).strip().upper()

            target_raw = str(
                row.get(
                    'Matches This Table As',
                    ''
                )
            ).strip()

            covered = str(
                row.get(
                    'Covered?',
                    ''
                )
            ).strip().lower()

            if (
                not gene
                or covered not in {'yes', 'partial'}
            ):
                continue

            match = re.match(
                r'^\s*(.*?)\s*'
                r'\((CYP_Reference|Transporter_Reference)\)'
                r'\s*$',
                target_raw
            )

            if match:
                gene_map[gene] = (
                    match.group(1).strip(),
                    match.group(2)
                )

        # These two aliases are relevant to the report but are not
        # currently separate rows in your workbook's Gene_Mapping tab.
        #
        # FDA groups CYP3A4 / CYP3A5 as CYP3A.
        # ABCG2 encodes the BCRP transporter.
        gene_map.setdefault(
            'CYP3A5',
            ('CYP3A', 'CYP_Reference')
        )

        gene_map.setdefault(
            'ABCG2',
            ('BCRP', 'Transporter_Reference')
        )

        _DDI_REFERENCE_CACHE = (
            cyp_df,
            transporter_df,
            gene_map
        )

    except Exception as exc:

        print(
            f"[WARN] Could not load DDI reference: "
            f"{exc}"
        )

        _DDI_REFERENCE_CACHE = empty

    return _DDI_REFERENCE_CACHE


def _ddi_direction(phenotype: str):
    """
    Determine which kind of interaction would compound
    the patient's genetic result.

    Reduced activity -> inhibitors
    Increased activity -> inducers
    Normal / unknown -> no DDI alert
    """

    p = str(
        phenotype or ''
    ).strip().lower()

    reduced_terms = (
        'intermediate',
        'poor',
        'decreased',
        'reduced',
        'slow',
        'low function',
        'no function',
    )

    increased_terms = (
        'ultrarapid',
        'rapid',
        'increased',
        'gain of function',
    )

    if any(
        term in p
        for term in reduced_terms
    ):
        return 'reduced'

    if any(
        term in p
        for term in increased_terms
    ):
        return 'increased'

    return None


def generate_compounding_alert_html(
    drug_name: str,
    gene_pheno_map: dict,
    gene_status_map: dict,
    fs_header: str = '14px',
    fs_body: str = '11.5px',
) -> str:
    """
    Create the:
        Other Medications That May Compound This Result

    box.

    The box is shown only when:

    1. The gene was actually analyzed.
    2. The patient has a non-normal phenotype.
    3. That gene/enzyme/transporter exists in the FDA DDI reference.
    4. At least one relevant inhibitor/inducer exists.
    """

    (
        cyp_df,
        transporter_df,
        gene_map
    ) = _load_ddi_reference()

    if not gene_map:
        return ""

    sections = []

    # Highest-significance interactions first.
    strength_rank = {
        'strong': 0,
        'moderate': 1,
        'weak': 2,
    }

    for gene in sorted(gene_pheno_map):

        # Only use a real patient genetic result.
        if (
            str(
                gene_status_map.get(
                    gene,
                    ''
                )
            ).strip().lower()
            != 'analyzed'
        ):
            continue

        gene_key = (
            str(gene)
            .strip()
            .upper()
        )

        mapping = gene_map.get(
            gene_key
        )

        if not mapping:
            continue

        target, table_name = mapping

        # A gene can occasionally have multiple possible
        # phenotypes because of phasing uncertainty.
        directions = {
            'reduced': [],
            'increased': [],
        }

        for (
            _,
            phenotype_text
        ) in gene_pheno_map.get(
            gene,
            []
        ):

            # Your existing code joins multiple phenotypes using " / "
            for phenotype in str(
                phenotype_text
            ).split(' / '):

                phenotype = (
                    phenotype.strip()
                )

                direction = (
                    _ddi_direction(
                        phenotype
                    )
                )

                if (
                    direction
                    and phenotype
                    not in directions[
                        direction
                    ]
                ):
                    directions[
                        direction
                    ].append(
                        phenotype
                    )

        for (
            direction,
            phenotypes
        ) in directions.items():

            if not phenotypes:
                continue

            # IMPORTANT:
            #
            # Reduced / slow metabolism:
            #     inhibitors can compound the result.
            #
            # Rapid / increased metabolism:
            #     inducers can compound the result.
            role = (
                'Inhibitor'
                if direction == 'reduced'
                else 'Inducer'
            )

            if (
                table_name
                == 'CYP_Reference'
            ):

                source_df = cyp_df
                target_col = 'Enzyme'

            else:

                source_df = (
                    transporter_df
                )
                target_col = (
                    'Transporter'
                )

            if source_df.empty:
                continue

            matches = source_df[
                source_df[
                    target_col
                ]
                .astype(str)
                .str.strip()
                .str.upper()
                .eq(
                    target.upper()
                )
                &
                source_df[
                    'Role'
                ]
                .astype(str)
                .str.strip()
                .str.lower()
                .eq(
                    role.lower()
                )
            ].copy()

            if matches.empty:
                continue

            # Strong -> Moderate -> Weak
            matches[
                '_strength_rank'
            ] = (
                matches[
                    'Strength'
                ]
                .astype(str)
                .str.strip()
                .str.lower()
                .map(
                    lambda x:
                    strength_rank.get(
                        x,
                        99
                    )
                )
            )

            matches = (
                matches.sort_values(
                    '_strength_rank',
                    kind='stable'
                )
            )

            medications = []

            current_drug = (
                str(drug_name)
                .strip()
                .casefold()
            )

            for raw_med in (
                matches['Drug']
                .tolist()
            ):

                med = str(
                    raw_med
                ).strip()

                if (
                    not med
                    or med.casefold()
                    == current_drug
                ):
                    continue

                display_med = med.title()

                if (
                    display_med
                    not in medications
                ):
                    medications.append(
                        display_med
                    )

                # Match the mock:
                # maximum 4 medications
                if len(
                    medications
                ) == 4:
                    break

            if not medications:
                continue

            phenotype_display = (
                ' / '.join(
                    phenotypes
                )
            )

            pills_html = ''.join(

                f'<span style="'
                f'background:#dbeafe; '
                f'color:#075985; '
                f'border:1px solid #93c5fd; '
                f'border-radius:12px; '
                f'padding:2px 8px; '
                f'font-size:11px; '
                f'font-weight:600; '
                f'display:inline-block; '
                f'margin:2px 4px 2px 0;'
                f'">{med}</span>'

                for med
                in medications
            )

            # Use the FDA table target name for display.
            # This matters because some report genes map to a different
            # FDA interaction-system label:
            #
            # CYP3A5  -> CYP3A
            # SLCO1B1 -> OATP1B1
            # ABCG2   -> BCRP
            system_label = target

            interaction_text = (
                f'Your <strong>{gene_key}</strong> result shows '
                f'<strong>{phenotype_display}</strong> status and may influence '
                f'how your body handles <strong>{str(drug_name).title()}</strong>. '
                f'The following are examples of FDA-listed medications or substances '
                f'that affect <strong>{system_label}</strong> and may further modify '
                f'drug exposure or response if taken together:'
            ) 

            sections.append(
                f'<div style="'
                f'font-size:{fs_body}; '
                f'color:#1e293b; '
                f'line-height:1.5; '
                f'margin-bottom:4px;'
                f'">'
                f'{interaction_text}'
                f'</div>'

                f'<div style="'
                f'line-height:1.5; '
                f'margin-bottom:5px;'
                f'">'
                f'{pills_html}'
                f'</div>'
            )

    if not sections:
        return ""

    # If more than one altered gene applies to the drug,
    # keep everything inside one DDI box instead of creating
    # several separate boxes.
    separator = (
        '<div style="'
        'height:1px; '
        'background:#bfdbfe; '
        'margin:6px 0;'
        '"></div>'
    )

    return (

        f'<div style="'
        f'background:#eff6ff; '
        f'border:1px solid #bfdbfe; '
        f'border-left:4px solid #0284c7; '
        f'border-radius:6px; '
        f'padding:7px 11px; '
        f'margin:0 0 10px 0; '
        f'page-break-inside:avoid; '
        f'break-inside:avoid;'
        f'">'

        f'<div style="'
        f'font-size:{fs_header}; '
        f'font-weight:700; '
        f'color:#075985; '
        f'margin-bottom:4px;'
        f'">'

        f'Other Medications That May '
        f'Further Influence This Result '

        f'</div>'

        f'{separator.join(sections)}'

        f'<div style="'
        f'font-size:10px; '
        f'color:#64748b; '
        f'font-style:italic; '
        f'line-height:1.4;'
        f'">'

        f'If you are currently taking '
        f'any of the above, mention this '
        f'to your prescriber. '

        f'This is general reference '
        f'information based on your '
        f'genetic result — it is not '
        f'a personalized check of your '
        f'current medications.'

        f'</div>'

        f'</div>'
    )

def sanitize_id(name):
    return re.sub(r'[^a-z0-9_-]', '_', str(name).lower().strip())

def icon_badge(bg_from, bg_to):
    return f"""
    <div style="
        width: 36px; height: 36px; border-radius: 50%;
        background: linear-gradient(135deg, {bg_from}, {bg_to});
        flex-shrink: 0;">
    </div>
    """


# ============================================================================
# GENOTYPE-LABEL DISPLAY TRANSLATION
# ============================================================================
# PharmCAT emits VKORC1 / IFNL3 diplotypes in raw rs-ID form, e.g.:
#   "rs9923231 reference (C)/rs9923231 reference (C)"
# This is the lookup key the GSI matches against, but it's unreadable for both
# the patient and the doctor.  These tables translate the raw diplotype into
# (pretty_diplotype, pretty_phenotype) pairs:
#   • Diplotype column: shows a compact rs-ID format like "rs9923231 C/C"
#   • Phenotype column: shows the haplotype interpretation, e.g. "-1639 G/G"
#     (VKORC1) or simply "C/C" (IFNL3).
# Heterozygous variants get a single canonical orientation regardless of the
# order PharmCAT produced.

_PRETTY_GENOTYPE_LABELS = {
    "VKORC1": {
        "rs9923231 reference (c)/rs9923231 reference (c)": ("rs9923231 C/C", "-1639 G/G"),
        "rs9923231 reference (c)/rs9923231 variant (t)":   ("rs9923231 C/T", "-1639 G/A"),
        "rs9923231 variant (t)/rs9923231 reference (c)":   ("rs9923231 C/T", "-1639 G/A"),
        "rs9923231 variant (t)/rs9923231 variant (t)":     ("rs9923231 T/T", "-1639 A/A"),
    },
    "IFNL3": {
        "rs12979860 reference (c)/rs12979860 reference (c)": ("rs12979860 C/C", "C/C"),
        "rs12979860 reference (c)/rs12979860 variant (t)":   ("rs12979860 C/T", "C/T"),
        "rs12979860 variant (t)/rs12979860 reference (c)":   ("rs12979860 C/T", "C/T"),
        "rs12979860 variant (t)/rs12979860 variant (t)":     ("rs12979860 T/T", "T/T"),
    },
    # rs12777823 — warfarin dosing in African Americans (CYP2B6-related locus)
    # PharmCAT / step4 emit chromosomal HGVS notation as the phenotype key.
    # Translate to a compact, human-readable form.
    "RS12777823": {
        "nc_000010.11:g.94645745g=/nc_000010.11:g.94645745g=":    ("G/G", "Reference (G/G) — no variant"),
        "nc_000010.11:g.94645745g=/nc_000010.11:g.94645745g>a":   ("G/A", "Heterozygous (G/A)"),
        "nc_000010.11:g.94645745g>a/nc_000010.11:g.94645745g=":   ("G/A", "Heterozygous (G/A)"),
        "nc_000010.11:g.94645745g>a/nc_000010.11:g.94645745g>a":  ("A/A", "Homozygous variant (A/A)"),
    },
}


def pretty_genotype_pair(gene: str, raw_diplotype: str):
    """
    For VKORC1 / IFNL3 / rs12777823, translate PharmCAT's raw diplotype /
    HGVS phenotype string to (pretty_diplotype, pretty_phenotype).
    Returns (raw, raw) for any other gene or unrecognised label.
    """
    g = str(gene).strip().upper()
    table = _PRETTY_GENOTYPE_LABELS.get(g)
    if not table:
        return raw_diplotype, raw_diplotype
    raw = str(raw_diplotype).strip().lower()
    if raw in table:
        return table[raw]
    return raw_diplotype, raw_diplotype


# ============================================================================
# STYLES
# ============================================================================

def styles(fonts, patient_name="Patient"):
    regular = fonts.get('regular', '')
    bold    = fonts.get('bold', '')
    italic  = fonts.get('italic', '')

    # Paged.js renders the physical footer through @page margin boxes.
    # Escape only the characters that can break a CSS quoted string.
    footer_label = f"Pharmacogenomics Report - {patient_name}"
    footer_label = footer_label.replace('\\', '\\\\').replace('\"', '\\"')

    return f"""
    @font-face {{
        font-family: 'DM Sans';
        src: url(data:font/truetype;charset=utf-8;base64,{regular}) format('truetype');
        font-weight: 400; font-style: normal;
    }}
    @font-face {{
        font-family: 'DM Sans';
        src: url(data:font/truetype;charset=utf-8;base64,{bold}) format('truetype');
        font-weight: 700; font-style: normal;
    }}
    @font-face {{
        font-family: 'DM Sans';
        src: url(data:font/truetype;charset=utf-8;base64,{italic}) format('truetype');
        font-weight: 400; font-style: italic;
    }}

    * {{ 
        margin: 0; 
        padding: 0; 
        box-sizing: border-box; 
        font-variant-ligatures: none;
        -webkit-font-variant-ligatures: none;
        font-feature-settings: "liga" 0, "clig" 0;
        text-rendering: optimizeSpeed;
    }}

html, body {{
        font-family: 'DM Sans', Arial, sans-serif;
        font-size: 11.5px;
        background: transparent;
        color: #1a1a1a;
        width: 100%;
        margin: 0;
        padding: 0;
    }}


    /* ================================================================
       PAGED.JS PHYSICAL PAGE MODEL
       ================================================================
       HTML wrappers represent report TOPICS, not fixed-height paper pages.
       Paged.js decides the real physical page breaks after layout.
    */
    @page {{
        size: A4;
        margin: 30px 52px 64px 52px;

        @bottom-left {{
            content: "{footer_label}";
            font-family: 'DM Sans', Arial, sans-serif;
            font-size: 9px;
            color: #444444;
            border-top: 1.5px solid #00B5C8;
            padding-top: 4px;
            vertical-align: top;
            text-align: left;
        }}

        @bottom-center {{
            content: counter(page);
            font-family: 'DM Sans', Arial, sans-serif;
            font-size: 9px;
            font-weight: 700;
            color: #444444;
            border-top: 1.5px solid #00B5C8;
            padding-top: 4px;
            vertical-align: top;
            text-align: center;
        }}

        @bottom-right {{
            content: "Table of Contents";
            font-family: 'DM Sans', Arial, sans-serif;
            font-size: 9px;
            font-style: italic;
            color: #00B5C8;
            border-top: 1.5px solid #00B5C8;
            padding-top: 4px;
            vertical-align: top;
            text-align: right;
        }}
    }}

    /* Page numbers are resolved AFTER layout from each link target. */
    .toc-page-number::after {{
        content: target-counter(attr(href), page);
    }}

    .toc-page-number {{
        color: #374151;
        font-weight: 700;
        text-decoration: none;
    }}

    p, li {{
        orphans: 3;
        widows: 3;
    }}

    /* ── A4 Page Setup ──
       NO fixed height — let Chromium paginate naturally so content that is
       too long for one page flows to the next instead of being clipped.
       The visible footer (cyan rule + report title + page number + TOC text)
       is rendered by Paged.js @page margin boxes. The bottom page margin is
       reserved by the paged-media layout itself, so topic wrappers do not need
       artificial footer padding.
    */
    /*
       IMPORTANT PAGED.JS RULE:
       A .page is a LOGICAL REPORT SECTION, not a fixed-height sheet.
       It may fragment across 2, 3, or more physical PDF pages.

       Force the START of every top-level report section onto a fresh page.
       We deliberately use break-before instead of break-after: when a section
       itself fragments across multiple pages, Paged.js can ignore/resolve an
       end break differently than expected. break-before on the NEXT section is
       deterministic and prevents the next drug/topic from being pulled into
       unused space on the previous section's last page.
    */
    .page {{
        width: 100%;
        min-height: 0;
        box-sizing: border-box;
        background: white;
        position: relative;
        padding-bottom: 0;

        page-break-inside: auto !important;
        break-inside: auto !important;

        page-break-after: auto !important;
        break-after: auto !important;
    }}

    /* First logical section (TOC) starts normally. Every later topic starts fresh. */
    body > .page:not(:first-child) {{
        page-break-before: always !important;
        break-before: page !important;
    }}

    .intro-page {{
        background: transparent !important;
        -webkit-print-color-adjust: exact;
        print-color-adjust: exact;
        padding-bottom: 0 !important;
    }}

    /* Keep an intro title attached to the first content block below it. */
    .intro-page .page-title {{
        page-break-after: avoid !important;
        break-after: avoid !important;
    }}

   /* ── Pagination Rules (Let Chromium flow content naturally) ── */
    .welcome-box, .faq-card-wrap, .drug-header {{
        page-break-inside: avoid;
        break-inside: avoid;
    }}

    .drug-header {{
        page-break-after: avoid !important;
        break-after: avoid !important;
    }}

    /*
       Keep each phenotype-specific recommendation entry together when it fits.
       Paged.js can otherwise fragment a nested bordered recommendation card at
       the page boundary and fail to paint the first fragment, which makes the
       beginning of the recommendation appear missing in the PDF.
       The outer CPIC/DPWG/FDA source block remains splittable.
    */
    .recommendation-entry {{
        page-break-inside: avoid !important;
        break-inside: avoid !important;
    }}

    .recommendation-source-label {{
        page-break-after: avoid !important;
        break-after: avoid !important;
    }}

    /* Keep the large content blocks fluid; do not force every row or heading onto
       a fresh paper page. Chromium should decide the physical break according to
       the actual rendered height. */
    /* Keep table rows intact across physical pages.  Tables themselves may
       flow for as many pages as needed, and thead is repeated by Paged.js. */
    tr {{
        page-break-inside: avoid !important;
        break-inside: avoid !important;
    }}

    td, th {{
        page-break-inside: auto;
        break-inside: auto;
    }}
    
    table {{
        page-break-inside: auto;
        break-inside: auto;
        table-layout: fixed;
        word-wrap: break-word;
        word-break: break-word;
    }}
    
    thead {{
        display: table-header-group;
        page-break-inside: auto !important;
        break-inside: auto !important;
    }}
    
    td, th {{
        overflow-wrap: break-word;
        word-wrap: break-word;
        word-break: break-word;
        hyphens: auto;
    }}
    tbody {{
        page-break-inside: auto;
        break-inside: auto;
    }}
    /* .info-box intentionally excluded — allow large content boxes to split across pages */

    /* Pagination-safe grid rows for static introductory tables.
       Paged.js can occasionally consume part of a native <tr> at a page
       boundary without painting it.  These block/grid rows are atomic, so a
       whole row moves to the next physical page instead of losing text. */
    .intro-safe-grid {{
        width: 100%;
        box-sizing: border-box;
    }}

    .intro-safe-grid-header,
    .intro-safe-row {{
        page-break-inside: avoid !important;
        break-inside: avoid !important;
        box-sizing: border-box;
    }}

    .intro-safe-grid-header {{
        page-break-after: avoid !important;
        break-after: avoid !important;
    }}

    .crosslisted-category-block {{
        page-break-inside: auto !important;
        break-inside: auto !important;
        margin-top: 10px;
    }}

    /* Keep the category heading attached to the repeated drug title. */
    .crosslisted-category-block .crosslisted-category-heading,
    .crosslisted-category-block .drug-header {{
        page-break-after: avoid !important;
        break-after: avoid !important;
    }}

    /* Keep each SMALL secondary-category box intact when it fits, but allow
       page breaks between boxes.  Do not make the whole cross-listed block
       atomic, because long category-specific About/WITM text may need to flow. */
    .crosslisted-about-box,
    .crosslisted-impact-box,
    .crosslisted-reference-box {{
        page-break-inside: avoid !important;
        break-inside: avoid !important;
    }}

    .intro-section {{
        page-break-inside: auto !important;
        break-inside: auto !important;
    }}

    .intro-heading {{
        page-break-after: avoid !important;
        break-after: avoid !important;
        break-inside: avoid !important;
    }}

    h1, h2, h3, .page-title, .section-block-header {{
        page-break-after: avoid !important;
        break-after: avoid !important;
        break-inside: avoid !important;
    }}

    /* ── Page Title ── */
    .page-title {{
        text-align: center;
        font-size: 16px;
        font-weight: 700;
        color: #023D79;
        margin-bottom: 6px;
        letter-spacing: -0.3px;
        padding-bottom: 10px;
        border-bottom: 2px solid #4DB7D0;
        width: 95%;
        margin-left: auto;
        margin-right: auto;
    }}
    .page-subtitle {{
        text-align: center;
        font-size: 11.5px;
        color: #555;
        margin-bottom: 28px;
    }}

    /* ── Welcome info boxes ── */
    .welcome-box {{
        display: flex;
        border: 1px solid #dde2ea;
        border-radius: 10px;
        overflow: hidden;
        margin-bottom: 18px;
        background: white;
        box-shadow: 0 1px 4px rgba(0,0,0,0.05);
    }}
    .welcome-box-icon {{
        width: 52px;
        min-height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
        font-size: 16px;
        flex-shrink: 0;
    }}
    .welcome-box-body {{
        padding: 16px 18px;
        flex: 1;
    }}
    .welcome-box-title {{
        font-size: 11.5px;
        font-weight: 700;
        color: #1a1a1a;
        margin-bottom: 8px;
    }}
    .welcome-box-text {{
        font-size: 11.5px;
        line-height: 1.65;
        color: #333;
    }}
    .welcome-box-text ul {{
        margin: 8px 0;
        padding-left: 18px;
    }}
    .welcome-box-text li {{ margin-bottom: 3px; }}

    /* ── FAQ Card Grid ── */
    .faq-grid {{
        display: grid;
        grid-template-columns: 1fr 1fr;
        gap: 20px 18px;
        margin-bottom: 18px;
    }}

    /* Card wrapper gives space for the outside number */
    .faq-card-wrap {{
        position: relative;
        padding-top: 16px;
        padding-left: 16px;
    }}

    /* The numbered circle — sits OUTSIDE the card, overlapping top-left */
    .faq-number {{
        position: absolute;
        top: 0;
        left: 0;
        width: 32px;
        height: 32px;
        border-radius: 50%;
        background: white;
        border: 2px solid #E91E8C;
        display: flex;
        align-items: center;
        justify-content: center;
        font-size: 12px;
        font-weight: 700;
        color: #E91E8C;
        z-index: 2;
    }}

    /* The actual card box */
    .faq-card {{
        background: white;
        border: 1.5px solid #dde2ea;
        border-radius: 10px;
        padding: 22px 16px 16px 22px;
        box-shadow: 0 1px 4px rgba(0,0,0,0.05);
        min-height: 108px;
    }}

    .faq-card-title {{
        font-size: 11.5px;
        font-weight: 700;
        color: #1a1a1a;
        margin-bottom: 7px;
    }}
    .faq-card-text {{
        font-size: 9.5px;
        line-height: 1.6;
        color: #444;
    }}

    .info-box {{
        border: 1px solid #d0e4f7;
        border-left: 4px solid #1a73e8;
        border-radius: 8px;
        padding: 14px 18px;
        margin-bottom: 14px;
        background: #f4f9ff;
    }}
    .info-box-title {{
        font-size: 11.5px;
        font-weight: 700;
        color: #1a1a1a;
        margin-bottom: 8px;
    }}
    .info-box-content {{
        font-size: 11.5px;
        line-height: 1.65;
        color: #333;
    }}
    .info-box-content ul {{
        margin: 6px 0;
        padding-left: 18px;
    }}
    .info-box-content li {{ margin-bottom: 4px; }}

    /* ── Drug page ── */
    .drug-header {{
        border: 1.5px solid #4DB7D0;                  /* teal outline */
        background: linear-gradient(90deg,#E7F7FB,#F5FCFF);  /* light teal inside */
        color: #000000;
        padding: 7px 12px;
        border-radius: 8px;
        margin-bottom: 6px;
        display: flex;
        justify-content: flex-start;
        align-items: center;
    }}
    .drug-title {{
        font-size: 12px;
        font-weight: 700;
        margin: 0;
        font-family: 'DM Sans', sans-serif;
        text-transform: lowercase;
    }}
    .drug-title::first-letter {{
        text-transform: uppercase;
    }}
    .drug-cat-tag {{ display: none; }}


    /* ── Gene table ── */
    .section-block {{
        border: 1px solid #dde2ea;
        border-radius: 8px;
        margin-bottom: 12px;
        /* overflow: hidden REMOVED to allow perfect pagination */
    }}
    .section-block-header {{
        background: #f3f5f8;
        padding: 7px 12px;
        font-size: 11.5px;
        font-weight: 700;
        color: #1a1a1a;
        border-bottom: 1px solid #dde2ea;
        border-top-left-radius: 7px;
        border-top-right-radius: 7px;
    }}

    table {{ 
            width: 100%; 
            border-collapse: collapse; 
            font-size: 11.5px; 
            table-layout: fixed; /* Forces strict column widths */
        }}
    th {{
        background: #0D3B7A;
        color: white;
        padding: 7px 10px;
        text-align: left;
        font-weight: 700;
    }}
    
    /* Apply rounding directly to headers so we don't need overflow:hidden */
    thead tr:first-child th:first-child {{ border-top-left-radius: 7px; }}
    thead tr:first-child th:last-child {{ border-top-right-radius: 7px; }}

    td {{ padding: 5px 8px; border-bottom: 1px solid #e9ecef; }}
    tbody tr:last-child td {{ border-bottom: none; }}


    /* ── TOC ── */
    .toc-category {{
        font-weight: 700;
        color: #0D3B7A;
        font-size: 11.5px;
        margin: 12px 0 5px 0;
        padding-left: 8px;
        border-left: 3px solid #0D3B7A;
    }}
    /* Dynamic page reference used by compact cross-listed drug entries.
       Paged.js resolves the physical page of the primary/full monograph. */
    a.cross-ref-page {{
        color: #0D3B7A;
        font-weight: 700;
        text-decoration: none;
    }}
    a.cross-ref-page::after {{
        content: target-counter(attr(href), page);
    }}

    .toc-row {{
        display: flex;
        justify-content: space-between;
        padding: 4px 8px;
        margin-bottom: 1px;
        border-bottom: 1px dotted #ddd;
        font-size: 11.5px;
    }}
    .toc-pg {{ color: #888; font-weight: 600; }}

    @media print {{
        .page {{ box-shadow: none; margin: 0; }}
    }}
    
    /* Keep the blue header rows borderless so only the dark-blue band shows. */
    table tr[style*="background:#4DB7D0"] th,
    table tr[style*="background:#023D79"] th,
    table th[style*="background:#4DB7D0"],
    table th[style*="background:#023D79"] {{
        border-color: transparent !important;
    }}


/* ================================================================
   UNIVERSAL REPORT PAGINATION
   ================================================================ */

/*
   Large sections are allowed to continue naturally onto another page.
   Chromium decides the break based on the ACTUAL rendered height.
*/
.flow-section {{
    page-break-inside: auto !important;
    break-inside: auto !important;
}}


/*
   Small logical units that should remain together.
*/
.keep-together {{
    page-break-inside: avoid !important;
    break-inside: avoid !important;
}}


/*
   Prevent a heading from being left alone at the bottom of a page.
*/
.keep-with-next {{
    page-break-after: avoid !important;
    break-after: avoid !important;
}}


/*
   Large cards can continue naturally onto another physical PDF page.
   Keep the box itself simple so Paged.js only has to fragment normal flow.
*/
.flow-card {{
    page-break-inside: auto !important;
    break-inside: auto !important;
    overflow: visible !important;
    max-height: none !important;
    height: auto !important;
}}


/*
   Keep the beginning of a card together:
   header + short introductory information.
*/
.flow-card-top {{
    page-break-inside: avoid !important;
    break-inside: avoid !important;
}}


/*
   Dynamic pagination safety for variable-length medication content.
   Nothing here depends on provider, gene, category, drug count, or page count.
*/
.oem-card {{
    /* Keep ordinary cards together. If a future card is taller than a page,
       its block-level medication rows below provide safe fallback break points. */
    page-break-inside: avoid !important;
    break-inside: avoid !important;
    overflow: visible !important;
    max-height: none !important;
    height: auto !important;
}}

.oem-card-head {{
    page-break-after: avoid !important;
    break-after: avoid !important;
}}

.oem-medication-list {{
    display: block;
    padding: 0;
    page-break-inside: auto !important;
    break-inside: auto !important;
    overflow: visible !important;
    orphans: 2;
    widows: 2;
}}

.oem-medication-chunk {{
    display: block;
    page-break-inside: avoid !important;
    break-inside: avoid !important;
    overflow: visible !important;
    margin: 0;
    padding: 0;
}}

/* Compact comma-separated OEM display.
   Each medication remains individually marked for PDF QA, but visually the
   list reads exactly like normal prose instead of a bullet list. */
.oem-medication-item {{
    display: inline;
    white-space: nowrap;
}}

.oem-card-count {{
    page-break-before: avoid !important;
    break-before: avoid !important;
}}

.specialized-gene-card {{
    page-break-inside: auto !important;
    break-inside: auto !important;
    overflow: visible !important;
    max-height: none !important;
    height: auto !important;
}}

.specialized-gene-top {{
    page-break-inside: avoid !important;
    break-inside: avoid !important;
}}

.specialized-affected-section {{
    page-break-inside: auto !important;
    break-inside: auto !important;
}}

.specialized-med-group {{
    page-break-inside: auto !important;
    break-inside: auto !important;
    margin: 0 0 7px 0;
}}

.specialized-med-category {{
    page-break-after: avoid !important;
    break-after: avoid !important;
    font-weight: 700;
    color: #0D3B7A;
    font-size: 11.5px;
    line-height: 1.35;
    margin-bottom: 2px;
}}

/* Compact comma-separated Specialized Testing display.
   IMPORTANT: the complete medication list is NOT one large fragmenting box.
   Python breaks the dynamic list into short invisible-layout chunks.  Each
   chunk is kept together, while the overall list may span unlimited pages.
   Visually this remains ordinary comma-separated text — no bullets. */
.specialized-med-list {{
    display: block;
    padding: 0;
    page-break-inside: auto !important;
    break-inside: auto !important;
    overflow: visible !important;
    color: #374151;
    font-size: 11.5px;
    line-height: 1.4;
}}

.specialized-med-first-fragment,
.specialized-med-chunk {{
    display: block;
    page-break-inside: avoid !important;
    break-inside: avoid !important;
    overflow: visible !important;
    margin: 0;
    padding: 0;
}}

.specialized-med-item {{
    display: inline;
    white-space: normal;
}}

/* ================================================================
   PAGINATION-SAFE VARIABLE DATA GRIDS
   ================================================================
   These are intentionally NOT native <table>/<tr> structures.  Chromium /
   Paged.js can occasionally consume the first fragment of a long table row at
   a physical page boundary.  Each data row below is its own normal block, so a
   row that fits on one page moves intact to the next page.  The outer section
   remains fully fluid and can span any number of physical pages.
*/
.genotype-grid-header,
.genotype-grid-row {{
    display: grid;
    grid-template-columns: 20% 30% 50%;
    width: 100%;
    box-sizing: border-box;
}}

.genotype-grid-header {{
    background: #174A86;
    color: #ffffff;
    border-radius: 7px 7px 0 0;
    page-break-after: avoid !important;
    break-after: avoid !important;
    page-break-inside: avoid !important;
    break-inside: avoid !important;
}}

.genotype-grid-row {{
    border-bottom: 1px solid #f3f4f6;
    page-break-inside: avoid !important;
    break-inside: avoid !important;
}}

.genotype-grid-cell {{
    padding: 6px 9px;
    min-width: 0;
    overflow-wrap: anywhere;
    word-break: normal;
    box-sizing: border-box;
}}

.analytical-grid-header,
.analytical-grid-row {{
    display: grid;
    grid-template-columns: 18% 82%;
    width: 100%;
    box-sizing: border-box;
}}

.analytical-grid-header {{
    background: #174A86;
    color: #ffffff;
    page-break-after: avoid !important;
    break-after: avoid !important;
    page-break-inside: avoid !important;
    break-inside: avoid !important;
}}

.analytical-grid-row {{
    border-bottom: 1px solid #fde68a;
    page-break-inside: avoid !important;
    break-inside: avoid !important;
}}

.analytical-grid-cell {{
    padding: 5px 8px;
    min-width: 0;
    overflow-wrap: anywhere;
    word-break: normal;
    box-sizing: border-box;
}}

.no-guideline-grid {{
    width: 100%;
    box-sizing: border-box;
}}

/*
   IMPORTANT: do not use CSS Grid for this multi-page list.
   Paged.js 0.4.3 can lose the leading fragment of a grid/table-like row
   when the parent is fragmented at a physical page boundary.

   Each row below is therefore a simple atomic block.  The two cells are
   inline-block children purely for the visual two-column layout.  Paged.js
   may break BETWEEN rows, but never needs to fragment a row itself.
*/
.no-guideline-grid-header,
.no-guideline-grid-row {{
    display: block;
    width: 100%;
    box-sizing: border-box;
    font-size: 0;
    page-break-inside: avoid !important;
    break-inside: avoid !important;
    -webkit-column-break-inside: avoid;
    overflow: hidden;
}}

.no-guideline-grid-header {{
    background: #023D79;
    color: #ffffff;
    border: 1px solid #023D79;
    page-break-after: avoid !important;
    break-after: avoid !important;
}}

.no-guideline-grid-row {{
    display: flex;
    flex-wrap: nowrap;
    border-left: 1px solid #dde2ea;
    border-right: 1px solid #dde2ea;
    border-bottom: 1px solid #e5e7eb;
    min-height: 34px;
}}

.no-guideline-grid-cell {{
    display: block;
    flex: 0 0 auto;
    vertical-align: top;
    padding: 5px 8px;
    min-width: 0;
    overflow-wrap: anywhere;
    word-break: normal;
    box-sizing: border-box;
    font-size: 12px;
}}

.no-guideline-grid-cell:first-child {{
    width: 40%;
    flex-basis: 40%;
}}

.no-guideline-grid-cell:last-child {{
    width: 60%;
    flex-basis: 60%;
}}

/*
   Generic multi-page table.
*/
.flow-table {{
    width: 100%;
    border-collapse: collapse;

    page-break-inside: auto !important;
    break-inside: auto !important;
}}

.flow-table thead {{
    display: table-header-group;
}}

.flow-table tbody {{
    page-break-inside: auto !important;
    break-inside: auto !important;
}}

.flow-table tr {{
    page-break-inside: avoid !important;
    break-inside: avoid !important;
}}

/* Universal heading behavior: never strand a heading at page bottom. */
.page-title,
.section-block-header,
.keep-with-next {{
    page-break-after: avoid !important;
    break-after: avoid !important;
}}

/* Genotype Summary + Analytical Notes are ONE continuous topic. */
.genotype-summary-section {{
    page-break-inside: auto !important;
    break-inside: auto !important;
}}

.analytical-notes-box {{
    page-break-inside: auto !important;
    break-inside: auto !important;
    -webkit-box-decoration-break: clone;
    box-decoration-break: clone;
}}

.analytical-notes-heading {{
    page-break-after: avoid !important;
    break-after: avoid !important;
}}
    """

def _wrap_page(
    content,
    patient_name,
    page_num,
    page_id="",
    background_color=None
):
    """
    Wrap report content in a page.

    background_color is now used only as a flag identifying
    introductory pages. The actual full-page tint is applied
    after PDF generation in utils.py.
    """

    _ = page_num

    pid = (
        f' id="{page_id}"'
        if page_id
        else ''
    )

    # Existing intro templates already pass background_color.
    # We use that only to assign the intro-page class.
    page_class = (
        ' intro-page'
        if background_color
        else ''
    )

    return f"""
    <div class="page{page_class}"{pid}>
        {content}
    </div>
    """


def _wrap_page_OLD_DEPRECATED(content, patient_name, page_num, page_id=""):
    """Kept for reference only \u2014 see _wrap_page above for current behavior."""
    pid = f' id="{page_id}"' if page_id else ''

    footer_label = f"{patient_name}\u2019s Personalized Medicine Report"

    toc_cell = '<div style="flex:1;"></div>'

    # Three-cell footer: patient label | centered page number | TOC link.
    # IMPORTANT: position must be 'absolute' (per .page wrapper), NOT 'fixed'.
    # 'position: fixed' causes Ghostscript / Chromium HTML-to-PDF to emit ONE
    # footer for the whole document — every page then shows the same number
    # (e.g. "53" on every page) because fixed elements are positioned relative
    # to the viewport, not the printed page.  'position: absolute; bottom:0'
    # within each .page (which itself has page-break-after:always) gives a
    # per-page footer with the correct page number.
    page_num_cell = (
        f'<div style="flex:1; text-align:center;">'
        f'<span style="font-size:12px; color:#444444; background:white; '
        f'padding:0 10px; font-weight:600;">{page_num}</span>'
        f'</div>'
    )
    page_footer = f"""
        <div style="position:absolute; bottom:0; left:0; right:0;
                    font-family: 'DM Sans', sans-serif; z-index:900;
                    background:white; padding:5px 0 6px 0;">
            <div style="height:1.5px; background-color:#00B5C8;
                        margin-bottom:5px; width:100%;"></div>
            <div style="display:flex; justify-content:space-between;
                        align-items:center;">
                <div style="flex:1; text-align:left;">
                    <span style="font-size:12px; color:#444444; background:white;
                                 padding-right:5px;">{footer_label}</span>
                </div>
                {page_num_cell}
                {toc_cell}
            </div>
        </div>
    """

    # padding-bottom raised from 30px → 60px so trailing content (citation
    # footers, "continued" rows) stops bleeding onto the footer line.
    # overflow:hidden clips any tail that still tries to escape.
    return f"""
    <div class="page" style="position:relative; min-height:889px;
                              box-sizing:border-box; padding-bottom:60px;
                              overflow:hidden;"{pid}>
        {content}
        {page_footer}
    </div>
    """

# ============================================================================
# 1. WELCOME PAGE
# ============================================================================

def welcome_template(name, pg):

    def _h2(text):
        return (
            f'<div style="'
            f'font-size:12px; '
            f'font-weight:700; '
            f'color:#023D79; '
            f'margin-top:10px; '
            f'margin-bottom:4px; '
            f'padding-bottom:3px; '
            f'border-bottom:1.5px solid #e5e7eb;'
            f'">'
            f'{text}'
            f'</div>'
        )


    def _para(text):
        return (
            f'<p style="'
            f'font-size:11.5px; '
            f'color:#374151; '
            f'line-height:1.6; '
            f'margin-bottom:6px;'
            f'">'
            f'{text}'
            f'</p>'
        )


    # ========================================================================
    # TABLE USED ON PAGE 2
    # ========================================================================

    pgx_rows = [

        (
            "Clopidogrel (antiplatelet)",

            "Alternative antiplatelet therapy may be recommended "
            "for CYP2C19 intermediate or poor metabolizers, "
            "especially in ACS/PCI settings."
        ),

        (
            "Warfarin (blood thinner)",

            "PGx-guided dosing can help estimate the warfarin "
            "starting dose when CYP2C9, VKORC1, CYP4F2, and "
            "clinical factors are available."
        ),

        (
            "SSRIs (antidepressants)",

            "CYP2C19 and CYP2D6 variants affect how drugs like "
            "escitalopram and sertraline are metabolized. Rapid "
            "or ultrarapid metabolism may reduce exposure to some "
            "antidepressants, which can affect response."
        ),

        (
            "Statins (cholesterol drugs)",

            "Variants in SLCO1B1, ABCG2, and CYP2C9 can affect "
            "statin exposure and the risk of statin-associated "
            "muscle symptoms; the level of concern differs by "
            "statin and dose."
        ),

        (
            "Thiopurines (chemotherapy / immunosuppression)",

            "TPMT and NUDT15 variants affect tolerance to active "
            "thiopurine metabolites and can increase the risk of "
            "severe myelosuppression."
        ),
    ]


    table_body = ""


    for i, (drug, guidance) in enumerate(
        pgx_rows
    ):

        bg = (
            "#f9fafb"
            if i % 2 == 0
            else "#ffffff"
        )


        table_body += (
            f'<div class="intro-safe-row" style="display:grid; grid-template-columns:30% 70%; background:{bg};">'
            f'<div style="padding:6px 10px; font-size:12px; font-weight:600; border-left:1px solid #e5e7eb; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#1a1a1a; box-sizing:border-box;">{drug}</div>'
            f'<div style="padding:6px 10px; font-size:12px; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#374151; line-height:1.5; box-sizing:border-box;">{guidance}</div>'
            f'</div>'
        )


    # ========================================================================
    # PAGE 1
    #
    # General PGx explanation + types of pharmacogenes.
    #
    # We deliberately STOP before "PGx in Everyday Clinical Practice".
    # This prevents the table from spilling silently onto another physical page.
    # ========================================================================

    page1_content = f"""
    <div style="padding-top:10px;">

        <div
            class="page-title"
            id="welcome"
        >
            About Personalized Medicine
        </div>


        {_h2(
            "What Is Personalized Medicine and Why Does It Matter?"
        )}

        {_para(
            "Medicine has traditionally been prescribed based on "
            "population averages: a standard drug at a standard dose "
            "for everyone with the same diagnosis. But we are not all "
            "the same. Your genes, lifestyle, age, and other biological "
            "factors shape how your body responds to a medication. "
            "Personalized medicine, sometimes called precision medicine, "
            "is an approach that uses these individual differences to "
            "guide treatment choices."
        )}

        {_para(
            "One of the most important drivers of this variability is "
            "pharmacogenomics (PGx): the study of how your genetic "
            "makeup influences how your body processes and responds "
            "to drugs. Even small differences in specific genes can "
            "change how quickly you break down a drug, how much of it "
            "reaches your bloodstream, and how strongly you respond to it."
        )}

        {_para(
            "Knowing this in advance helps your doctor make smarter "
            "choices &mdash; potentially avoiding drugs that won't work "
            "well for you, choosing safer alternatives, or adjusting "
            "doses before problems arise."
        )}


        {_h2(
            "Where the Science Stands"
        )}

        {_para(
            "Pharmacogenomics is not new; it has been studied for decades, "
            "but clinical adoption has accelerated dramatically. Today, "
            "major health institutions around the world use PGx testing "
            "to guide prescribing in oncology, psychiatry, cardiology, "
            "infectious disease, and more."
        )}

        {_para(
            "Authoritative bodies like CPIC (Clinical Pharmacogenetics "
            "Implementation Consortium), ClinPGx (formerly PharmGKB), "
            "and the FDA have published evidence-based guidelines that "
            "translate genetic findings directly into prescribing "
            "recommendations."
        )}

        {_para(
            "At the same time, the field is still growing. Not every "
            "drug has robust PGx evidence yet, and not every genetic "
            "variant has been fully characterized. This is why this "
            "report is structured carefully: only gene-drug combinations "
            "with established clinical guidelines are presented as "
            "actionable recommendations. Everything else is shown "
            "transparently, with clear explanations of what is and "
            "isn't yet known."
        )}


        {_h2(
            "The Different Types of Pharmacogenes"
        )}

        {_para(
            "Not all pharmacogenomic genes work the same way. This report "
            "includes several types of genes, each of which influences "
            "drug response through a different mechanism. Understanding "
            "this helps you interpret the phenotype labels you will encounter:"
        )}


        <ul style="
            font-size:11.5px;
            color:#374151;
            line-height:1.6;
            margin-top:4px;
            margin-bottom:8px;
            padding-left:18px;
        ">

            <li style="margin-bottom:4px;">
                <strong>
                    1. Drug-Metabolizing Enzymes
                    (CYP genes, DPYD, TPMT, NUDT15, UGT1A1)
                </strong>:
                These genes encode enzymes that chemically break down drugs
                in the liver or gut. Variants affect how quickly a drug is
                cleared from the body. These are described using the
                metabolizer phenotype framework
                (Normal, Intermediate, Poor, Rapid/Ultrarapid Metabolizer).
            </li>

            <li style="margin-bottom:4px;">
                <strong>
                    2. Drug Transporter Genes
                    (SLCO1B1, ABCG2)
                </strong>:
                These genes encode proteins that physically move drugs across
                cell membranes. Variants are described as Decreased Function
                or Poor Function because the drug is not being metabolized
                differently &mdash; it is being transported differently.
            </li>

            <li style="margin-bottom:4px;">
                <strong>
                    3. Drug Target / Mechanism Genes
                    (VKORC1, RYR1, CACNA1S)
                </strong>:
                These encode the molecular targets or mechanisms relevant to
                drug response. Results may therefore use labels other than
                metabolizer phenotypes.
            </li>

            <li style="margin-bottom:4px;">
                <strong>
                    4. Immune and Hypersensitivity Genes
                    (HLA-A, HLA-B)
                </strong>:
                These genes are involved in immune recognition. Certain
                alleles are associated with severe immune reactions to
                specific medications and may require specialized
                high-resolution typing.
            </li>

            <li>
                <strong>
                    5. Response Prediction Genes
                    (IFNL3)
                </strong>:
                Some genes help predict the likelihood of treatment response
                rather than primarily predicting toxicity.
            </li>

        </ul>

    </div>
    """


    # ========================================================================
    # PAGE 2
    #
    # Clinical examples + safety reminder.
    # ========================================================================

    page2_content = f"""
    <div style="padding-top:10px;">


        {_h2(
            "PGx in Everyday Clinical Practice"
        )}

        {_para(
            "Pharmacogenomics-guided prescribing is increasingly used "
            "in clinical settings. The examples below illustrate how "
            "genetic information can contribute to medication decisions."
        )}


        <div class="intro-safe-grid" style="margin-bottom:14px;">
            <div class="intro-safe-grid-header" style="display:grid; grid-template-columns:30% 70%; background:#4DB7D0;">
                <div style="padding:5px 8px; font-size:12px; font-weight:700; text-align:left; color:white; box-sizing:border-box;">Drug / Drug Class</div>
                <div style="padding:5px 8px; font-size:12px; font-weight:700; text-align:left; color:white; box-sizing:border-box;">How PGx Can Guide Prescribing</div>
            </div>
            {table_body}
        </div>


        <div style="
            border:1px solid #d97706;
            border-left:4px solid #d97706;
            border-radius:6px;
            background:#fffbeb;
            padding:10px 16px;
            page-break-inside:avoid;
            break-inside:avoid;
        ">

            <div style="
                font-size:11.5px;
                font-weight:700;
                color:#92400e;
                margin-bottom:4px;
            ">
                Important to Know
            </div>

            <div style="
                font-size:11.5px;
                line-height:1.6;
                color:#374151;
            ">
                PGx information is one input into a prescribing decision,
                not a prescription in itself. Your doctor combines this
                genetic information with your medical history, current
                medications, organ function, and clinical judgment.
                This report is designed to support that conversation.
            </div>

        </div>

    </div>
    """


    # Let Chromium fill physical pages naturally inside this section.
    # We keep one logical wrapper so there are no artificial half-empty
    # continuation pages.
    combined_content = page1_content + page2_content

    page = _wrap_page(
        combined_content,
        name,
        pg,
        page_id="welcome",
        background_color="#f3f4f6"
    )

    return [page], pg + 1

# ============================================================================
# 2. HOW TO READ PAGE
# ============================================================================

def how_to_read_template(name, pg):

    # ========================================================================
    # Helpers
    # ========================================================================

    def _h2(text):
        return (
            f'<div style="'
            f'font-size:12px; '
            f'font-weight:700; '
            f'color:#023D79; '
            f'margin-top:10px; '
            f'margin-bottom:4px; '
            f'padding-bottom:3px; '
            f'border-bottom:1.5px solid #e5e7eb;'
            f'">'
            f'{text}'
            f'</div>'
        )


    def _para(text):
        return (
            f'<p style="'
            f'font-size:11.5px; '
            f'color:#374151; '
            f'line-height:1.6; '
            f'margin-top:0; '
            f'margin-bottom:6px;'
            f'">'
            f'{text}'
            f'</p>'
        )


    def _continued_title():
        # Continuation headings are intentionally omitted. Chromium paginates
        # this logical section using the actual rendered height.
        return ""


    def _table2(
        h1,
        h2,
        w1,
        rows
    ):

        body = ""


        for i, (a, b) in enumerate(
            rows
        ):

            bg = (
                "#f9fafb"
                if i % 2 == 0
                else "#ffffff"
            )


            body += (
                f'<div class="intro-safe-row" style="display:grid; grid-template-columns:{w1} 1fr; background:{bg};">'
                f'<div style="padding:5px 9px; font-size:11.5px; font-weight:600; border-left:1px solid #e5e7eb; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#1a1a1a; box-sizing:border-box;">{a}</div>'
                f'<div style="padding:5px 9px; font-size:11.5px; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#374151; line-height:1.45; box-sizing:border-box;">{b}</div>'
                f'</div>'
            )


        return f"""
        <div class="intro-safe-grid" style="margin-bottom:10px;">
            <div class="intro-safe-grid-header" style="display:grid; grid-template-columns:{w1} 1fr; background:#4DB7D0;">
                <div style="padding:5px 8px; font-size:11.5px; font-weight:700; color:white; box-sizing:border-box;">{h1}</div>
                <div style="padding:5px 8px; font-size:11.5px; font-weight:700; color:white; box-sizing:border-box;">{h2}</div>
            </div>
            {body}
        </div>
        """


    def _table3(
        h1,
        h2,
        h3,
        w1,
        w2,
        rows
    ):

        body = ""


        for i, (a, b, c) in enumerate(
            rows
        ):

            bg = (
                "#f9fafb"
                if i % 2 == 0
                else "#ffffff"
            )


            body += (
                f'<div class="intro-safe-row" style="display:grid; grid-template-columns:{w1} {w2} 1fr; background:{bg};">'
                f'<div style="padding:5px 8px; font-size:11.5px; font-weight:600; border-left:1px solid #e5e7eb; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#1a1a1a; box-sizing:border-box;">{a}</div>'
                f'<div style="padding:5px 8px; font-size:11.5px; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#374151; line-height:1.4; box-sizing:border-box;">{b}</div>'
                f'<div style="padding:5px 8px; font-size:11.5px; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#374151; line-height:1.4; box-sizing:border-box;">{c}</div>'
                f'</div>'
            )


        return f"""
        <div class="intro-safe-grid" style="margin-bottom:10px;">
            <div class="intro-safe-grid-header" style="display:grid; grid-template-columns:{w1} {w2} 1fr; background:#4DB7D0;">
                <div style="padding:5px 8px; font-size:11.5px; font-weight:700; color:white; box-sizing:border-box;">{h1}</div>
                <div style="padding:5px 8px; font-size:11.5px; font-weight:700; color:white; box-sizing:border-box;">{h2}</div>
                <div style="padding:5px 8px; font-size:11.5px; font-weight:700; color:white; box-sizing:border-box;">{h3}</div>
            </div>
            {body}
        </div>
        """


    # ========================================================================
    # Content tables
    # ========================================================================

    content_rows = [

        (
            "Detailed Drug Reports",

            "For each drug with an established PGx guideline, you will find "
            "the genes analyzed, your genotype and phenotype, a plain-language "
            "explanation of what it means for you, and the clinical "
            "recommendation from recognized guidelines."
        ),

        (
            "Other Evaluated Medications",

            "Drugs that were evaluated but do not have a specific clinical "
            "action or recommendation in this report. This may occur because "
            "the relevant genetic result was unavailable or because the "
            "available evidence does not produce a patient-specific clinical "
            "recommendation."
        ),

        (
            "No Guideline Available",

            "Drugs for which a relevant gene-drug relationship was identified "
            "but no formal prescribing recommendation is available for the "
            "specific result represented in this report."
        ),

        (
            "Genes Requiring Specialized Testing",

            "Genes that could not be reliably assessed from the submitted DNA "
            "data, together with the reason and the medications for which the "
            "missing information may be relevant."
        ),

        (
            "Genotype Summary",

            "A consolidated reference of the genes evaluated, reported "
            "diplotypes, and resulting phenotypes for this sample."
        ),
    ]


    key_terms = [

        (
            "Gene",

            "The specific gene being analyzed. Genes provide instructions "
            "for producing enzymes, transporters, receptors, and other "
            "proteins that can influence medication response."
        ),

        (
            "Diplotype",

            "Your reported combination of genetic alleles for that gene "
            "(for example, *1/*2). The exact notation depends on the gene "
            "and the allele-definition system used."
        ),

        (
            "Phenotype",

            "The predicted functional interpretation of the genetic result. "
            "Examples include Normal Metabolizer, Intermediate Metabolizer, "
            "Poor Metabolizer, Normal Function, or other gene-specific "
            "functional labels."
        ),
    ]


    sources = [

        (
            "CPIC",

            "Peer-reviewed pharmacogenomic guidelines developed to help "
            "clinicians use available genetic test results in prescribing "
            "decisions.",

            "Recommendation strength and evidence are shown where available."
        ),

        (
            "DPWG",

            "Evidence-based pharmacogenomic recommendations developed by "
            "the Dutch Pharmacogenetics Working Group.",

            "Used where relevant recommendations are available."
        ),

        (
            "ClinPGx",

            "A curated resource of pharmacogenomic gene-drug associations, "
            "guideline information, and supporting evidence.",

            "Evidence annotations may be displayed where available."
        ),

        (
            "FDA",

            "Drug-label and pharmacogenomic information published by the "
            "U.S. Food and Drug Administration.",

            "FDA testing or actionable-PGx information may be displayed "
            "when applicable."
        ),
    ]


    # ========================================================================
    # PAGE 1
    #
    # Report structure + distinction between the two non-detailed sections.
    # ========================================================================

    page1_content = f"""
    <div style="padding-top:10px;">

        <div
            class="page-title"
            id="how_to_read"
        >
            About This Report
        </div>


        {_h2(
            "What This Report Contains"
        )}

        {_para(
            "This is a pharmacogenomics (PGx) report that examines how "
            "genetic variation may influence response to a range of "
            "medications. The report is organized into the following sections:"
        )}


        {_table2(
            "Section",
            "What It Tells You",
            "32%",
            content_rows
        )}



        {_h2(
            "The Difference Between "
            "'No Guideline Available' and "
            "'Other Evaluated Medications'"
        )}

        {_para(
            "These two sections serve different purposes and should not "
            "be interpreted as equivalent."
        )}

        {_para(
            "<strong>No Guideline Available:</strong> "
            "A relevant gene-drug relationship was identified, but the "
            "available guideline information does not provide a specific "
            "patient-management recommendation for the result shown in "
            "this report."
        )}

        {_para(
            "<strong>Other Evaluated Medications:</strong> "
            "The medication was evaluated against the available genetic "
            "information but is not presented as a detailed actionable "
            "drug report. This may occur because the relevant genetic "
            "result was unavailable or because the evaluated result does "
            "not generate a patient-specific clinical recommendation."
        )}


        <div style="
            border:1px solid #4DB7D0;
            border-left:4px solid #4DB7D0;
            border-radius:6px;
            background:#f0fffe;
            padding:9px 14px;
            margin-top:10px;
            page-break-inside:avoid;
            break-inside:avoid;
        ">

            <div style="
                font-size:11.5px;
                font-weight:700;
                color:#023D79;
                margin-bottom:3px;
            ">
                Important
            </div>

            <div style="
                font-size:11.5px;
                line-height:1.55;
                color:#374151;
            ">
                A medication appearing in either of these sections should
                not be interpreted as safe, unsafe, effective, or
                ineffective based on genetics alone. Medication decisions
                should always incorporate the full clinical picture.
            </div>

        </div>

    </div>
    """


    # ========================================================================
    # PAGE 2
    #
    # Key terminology + *1/reference explanation + evidence sources.
    # ========================================================================

    page2_content = f"""
    <div style="padding-top:10px;">

        {_continued_title()}


        {_h2(
            "What Do the Key Terms Mean?"
        )}

        {_para(
            "Drug-specific entries use genetic terminology to describe "
            "what was evaluated and how the result was interpreted."
        )}


        {_table2(
            "Term",
            "Plain-Language Explanation",
            "20%",
            key_terms
        )}


        {_h2(
            "About *1 and Reference Results"
        )}

        {_para(
            "A *1 or Reference result means that no defining alternative "
            "allele was identified at the relevant positions available "
            "for analysis. Because raw DNA files may not contain every "
            "position used to define an allele, a *1 or Reference result "
            "does not mean that the entire gene was sequenced or that all "
            "possible variants were excluded. Results should be interpreted "
            "within the marker coverage of the submitted DNA file."
        )}


        {_h2(
            "Where Do the Recommendations Come From?"
        )}

        {_para(
            "Clinical recommendations and evidence annotations in this "
            "report may be derived from the following recognized "
            "pharmacogenomic resources:"
        )}


        {_table3(
            "Source",
            "What It Provides",
            "How It Is Used in This Report",
            "18%",
            "47%",
            sources
        )}

    </div>
    """


    # ========================================================================
    # PAGE 3
    #
    # GSI + analysis method + currency + limitations.
    # ========================================================================

    page3_content = f"""
    <div style="padding-top:10px;">

        {_continued_title()}


        {_h2(
            "A Note on the GSI Catalog"
        )}

        {_para(
            "The GSI catalog is an additional curated reference used by "
            "this report-generation pipeline to identify genes that have "
            "reported associations with a medication beyond the gene-drug "
            "pairs for which a formal prescribing recommendation is shown. "
            "These associations may relate to metabolism, transport, "
            "toxicity, treatment response, or another pharmacogenomic "
            "mechanism."
        )}

        {_para(
            "Genes listed as additional GSI-associated genes are included "
            "for informational completeness. Their presence does not mean "
            "that PharmCAT analyzed that gene-drug pair or that a clinical "
            "prescribing recommendation is available."
        )}


        {_h2(
            "How Were Your Results Generated?"
        )}

        {_para(
            "Your results were generated from the raw DNA data file "
            "submitted for analysis. Available pharmacogenomic markers "
            "were processed through the report pipeline, including "
            "<strong>PharmCAT v3.2.0</strong>, to derive supported genetic "
            "calls and phenotype interpretations for applicable genes."
        )}

        {_para(
            "The resulting genetic information was then matched against "
            "available pharmacogenomic recommendations and supporting "
            "reference data to determine which medication findings could "
            "be presented in this report."
        )}

        {_h2(
                    "How This Report Is Structured"
                )}
        
                {_para(
                    "The report begins with guidance on interpreting "
                    "pharmacogenomic results, followed by a summary of medication "
                    "insights and detailed recommendations organized by therapeutic "
                    "category. Medications relevant to multiple categories are "
                    "grouped together, with a single shared recommendation, while "
                    "the Table of Contents lists them under each applicable "
                    "category and page number. Additional sections include other "
                    "evaluated medications, genes requiring specialized testing, "
                    "and a consolidated genotype summary."
                )}
        


        {_h2(
            "Guideline Currency"
        )}

        {_para(
            "Pharmacogenomic guidelines and drug-label information can "
            "change as new evidence becomes available. The recommendations "
            "shown in this report reflect the guideline and reference data "
            "used by the analysis pipeline at the time the report was "
            "generated. For high-stakes prescribing decisions, clinicians "
            "should confirm that they are using the most current applicable "
            "guideline or prescribing information."
        )}


        <div style="
            border:1px solid #dc2626;
            border-left:4px solid #dc2626;
            border-radius:6px;
            background:#fef2f2;
            padding:10px 16px;
            margin-top:12px;
            page-break-inside:avoid;
            break-inside:avoid;
        ">

            <div style="
                font-size:11.5px;
                font-weight:700;
                color:#991b1b;
                margin-bottom:4px;
            ">
                Technical Limitations
            </div>

            <div style="
                font-size:11.5px;
                line-height:1.6;
                color:#374151;
            ">

                Coverage varies by gene and by the raw DNA platform used.
                Some pharmacogenes can be assessed when the required
                markers are present, while others cannot be fully or
                reliably characterized from standard SNP-array data.

                <br><br>

                Structurally complex genes, copy-number changes,
                high-resolution HLA alleles, certain mitochondrial
                variants, and rare variants may require dedicated
                clinical testing or sequencing methods.

                <br><br>

                Any gene that could not be reliably assessed for this
                sample is identified separately in the
                <strong>Genes Requiring Specialized Testing</strong>
                section. A missing or indeterminate result should not be
                interpreted as a normal genetic result.

            </div>

        </div>

    </div>
    """


    # Keep the complete section in one logical wrapper. Chromium may use
    # as many physical PDF pages as needed and will use the remaining space
    # on each page before continuing.
    combined_content = (
        page1_content
        + page2_content
        + page3_content
    )

    page = _wrap_page(
        combined_content,
        name,
        pg,
        page_id="how_to_read",
        background_color="#f3f4f6"
    )

    return [page], pg + 1

# ============================================================================
# 3. FAQs PAGE
# ============================================================================

def faqs_template(name, pg):
    def _h2(text):
        return (f'<div style="font-size:12px; font-weight:700; color:#023D79; '
                f'margin-top:10px; margin-bottom:4px; padding-bottom:3px; '
                f'border-bottom:1.5px solid #e5e7eb;">{text}</div>')

    def _para(text):
        return f'<p style="font-size:11.5px; color:#374151; line-height:1.6; margin-bottom:6px;">{text}</p>'

    def _continued_title():
        # No forced continuation title. The section flows naturally across
        # physical PDF pages.
        return ""

    pheno_rows = [
        ("Normal (Extensive) Metabolizer",
         "Your enzyme works as expected for the general population. A 'Normal' result does not mean you are immune to side effects &mdash; it means genetics are unlikely to cause an unexpected response at standard doses.",
         "Standard doses are generally appropriate."),
        ("Intermediate Metabolizer",
         "Your enzyme works at a reduced rate, typically because you carry one fully functional copy and one reduced-function copy of the gene.",
         "Drug levels may be somewhat higher than average; dose adjustments or slower titration may be considered for some drugs."),
        ("Poor Metabolizer",
         "Your enzyme has very little or no activity, often because both copies of the gene carry loss-of-function variants.",
         "Drug may accumulate to high levels (increasing side-effect risk) or, if it requires activation by the enzyme, may have reduced effect."),
        ("Rapid or Ultrarapid Metabolizer",
         "Your enzyme works faster than normal, clearing the drug more quickly.",
         "Drug may be less effective at standard doses; higher doses or alternative drugs may be needed. For prodrugs, rapid conversion may lead to unexpectedly strong effects."),
        ("Indeterminate",
         "The available data are insufficient or ambiguous to assign a clear functional phenotype.",
         "No recommendation is generated from an indeterminate result. Further testing may be appropriate if the drug is clinically important."),
        ("<b>Possible</b> [Phenotype]",
         "The diplotype contains a variant with uncertain or incompletely characterized functional impact, preventing a fully definitive assignment (e.g., CYP3A5 *1/*9).",
         "The recommendation is still provided based on the most likely interpretation, but clinicians should be aware of reduced certainty and may wish to use additional clinical monitoring or seek confirmatory testing.")
    ]
    pheno_body = ""
    for i, (pheno, meaning, implication) in enumerate(pheno_rows):
        bg = "#f9fafb" if i % 2 == 0 else "#ffffff"
        pheno_body += (f'<div class="intro-safe-row" style="display:grid; grid-template-columns:24% 36% 40%; background:{bg};">'
                       f'<div style="padding:5px 8px; font-size:12px; font-weight:600; border-left:1px solid #e5e7eb; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#1a1a1a; box-sizing:border-box;">{pheno}</div>'
                       f'<div style="padding:5px 8px; font-size:12px; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#374151; line-height:1.5; box-sizing:border-box;">{meaning}</div>'
                       f'<div style="padding:5px 8px; font-size:12px; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#374151; line-height:1.5; box-sizing:border-box;">{implication}</div>'
                       f'</div>')

    transporter_rows = [
        ("Normal Function",
         "Transporter moves drug into target tissue at the expected rate.",
         "Standard drug exposure; no transporter-based dose concern."),
        ("Decreased Function",
         "Transporter moves drug less efficiently. This phenotype reflects reduced transporter activity, which may alter drug distribution and increase systemic drug exposure.",
         "Moderately increased drug plasma concentrations; some drugs may require dose adjustment or monitoring."),
        ("Poor Function",
         "Transporter has very little or no activity (two loss-of-function copies).",
         "Substantially higher drug exposure; dose reduction or alternative drug is often recommended.")
    ]
    transporter_body = ""
    for i, (pheno, meaning, implication) in enumerate(transporter_rows):
        bg = "#f9fafb" if i % 2 == 0 else "#ffffff"
        transporter_body += (f'<div class="intro-safe-row" style="display:grid; grid-template-columns:24% 36% 40%; background:{bg};">'
                       f'<div style="padding:5px 8px; font-size:12px; font-weight:600; border-left:1px solid #e5e7eb; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#1a1a1a; box-sizing:border-box;">{pheno}</div>'
                       f'<div style="padding:5px 8px; font-size:12px; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#374151; line-height:1.5; box-sizing:border-box;">{meaning}</div>'
                       f'<div style="padding:5px 8px; font-size:12px; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#374151; line-height:1.5; box-sizing:border-box;">{implication}</div>'
                       f'</div>')

    susc_rows = [
        ("Malignant Hyperthermia Susceptibility", "RYR1, CACNA1S",
         "Variants in the ryanodine receptor or calcium channel gene create abnormally unstable calcium channels in skeletal muscle. When exposed to triggering agents, these can release massive amounts of calcium, causing a catastrophic hypermetabolic crisis.",
         "Halogenated volatile anesthetics and succinylcholine are relatively contraindicated. This finding must be communicated to any anesthesiologist before any procedure."),
        ("Uncertain Susceptibility", "CACNA1S",
         "A report may show Uncertain Susceptibility when available data cannot fully exclude clinically relevant CACNA1S variation. This is <strong>different</strong> from a confirmed susceptibility.",
         "If an RYR1 susceptibility is already confirmed, the same anesthetic precautions already apply. The CACNA1S 'uncertain' result means the standard SNP array cannot rule out structural or rare CACNA1S variants. If independent confirmation is needed, in vitro contracture testing (IVCT) can be requested."),
        ("Adverse Reaction Risk", "HLA-B",
         "Certain HLA-B alleles are strongly associated with severe cutaneous adverse reactions (Stevens-Johnson syndrome / toxic epidermal necrolysis) with specific drugs.",
         "High-resolution HLA allele status generally cannot be established reliably from a standard consumer SNP-array file. When an HLA-B result is clinically required for prescribing, dedicated validated HLA testing may be necessary.")
    ]
    susc_body = ""
    for i, (pheno, genes, meaning, implication) in enumerate(susc_rows):
        bg = "#f9fafb" if i % 2 == 0 else "#ffffff"
        susc_body += (f'<div class="intro-safe-row" style="display:grid; grid-template-columns:20% 10% 40% 30%; background:{bg};">'
                       f'<div style="padding:5px 8px; font-size:12px; font-weight:600; border-left:1px solid #e5e7eb; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#1a1a1a; box-sizing:border-box;">{pheno}</div>'
                       f'<div style="padding:5px 8px; font-size:12px; font-weight:600; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#1e40af; box-sizing:border-box;">{genes}</div>'
                       f'<div style="padding:5px 8px; font-size:12px; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#374151; line-height:1.5; box-sizing:border-box;">{meaning}</div>'
                       f'<div style="padding:5px 8px; font-size:12px; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#374151; line-height:1.5; box-sizing:border-box;">{implication}</div>'
                       f'</div>')

    variant_rows = [
        ("VKORC1", "Reported as: -1639 G/A (rs9923231 C/T)",
         "VKORC1 is the direct molecular target of warfarin. The -1639 variant reduces VKORC1 expression, making patients more sensitive to vitamin K antagonists (warfarin, acenocoumarol, phenprocoumon). The result is a variant call, not a metabolizer phenotype."),
        ("CYP4F2", "Reported as diplotype (e.g., *1/*5)",
         "CYP4F2 metabolizes vitamin K1 (phylloquinone). Reduced CYP4F2 activity means less vitamin K is broken down, resulting in higher vitamin K levels and a slightly higher warfarin dose requirement."),
        ("IFNL3", "Reported as rs12979860 genotype (C/C, C/T, or T/T)",
         "IFNL3 (also called IL28B) does not metabolize drugs. It encodes interferon lambda, and its variant status predicts likelihood of response to interferon-based hepatitis C therapies. C/C = most favorable response; T/T = least favorable; C/T (heterozygous) = intermediate.")
    ]
    variant_body = ""
    for i, (gene, pfmt, meaning) in enumerate(variant_rows):
        bg = "#f9fafb" if i % 2 == 0 else "#ffffff"
        variant_body += (
            f'<div class="intro-safe-row" style="display:grid; grid-template-columns:15% 35% 50%; background:{bg};">'
            f'<div style="padding:5px 8px; font-size:12px; font-weight:600; border-left:1px solid #e5e7eb; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#1e40af; box-sizing:border-box;">{gene}</div>'
            f'<div style="padding:5px 8px; font-size:12px; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#374151; line-height:1.5; box-sizing:border-box;">{pfmt}</div>'
            f'<div style="padding:5px 8px; font-size:12px; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#374151; line-height:1.5; box-sizing:border-box;">{meaning}</div>'
            f'</div>'
        )

    dpyd_rows = [
        ("2.0", "Normal Metabolizer", "Standard dosing is appropriate."),
        ("1.5", "Intermediate Metabolizer", "Reduce starting dose by 50%; titrate based on clinical response and therapeutic drug monitoring."),
        ("1.0", "Intermediate Metabolizer", "Reduce starting dose by 50%; titrate based on clinical response and therapeutic drug monitoring."),
        ("0.5", "Poor Metabolizer", "Avoid fluoropyrimidines. If no alternative exists, use &ge;75% dose reduction with intensive monitoring."),
        ("0.0", "Poor Metabolizer", "Fluoropyrimidines are contraindicated. Avoid entirely or use only under exceptional circumstances with substantial dose reduction and intensive monitoring.")
    ]
    dpyd_body = ""
    for i, (asc, pheno, implication) in enumerate(dpyd_rows):
        bg = "#f9fafb" if i % 2 == 0 else "#ffffff"
        dpyd_body += (f'<div class="intro-safe-row" style="display:grid; grid-template-columns:15% 35% 50%; background:{bg};">'
                       f'<div style="padding:5px 8px; font-size:12px; font-weight:600; text-align:center; border-left:1px solid #e5e7eb; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#1a1a1a; box-sizing:border-box;">{asc}</div>'
                       f'<div style="padding:5px 8px; font-size:12px; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#374151; line-height:1.5; box-sizing:border-box;">{pheno}</div>'
                       f'<div style="padding:5px 8px; font-size:12px; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#374151; line-height:1.5; box-sizing:border-box;">{implication}</div>'
                       f'</div>')

    rec_rows = [
        ("Dose adjustment",
         "The evidence suggests starting at a lower or higher dose, titrating more slowly, or setting a maximum daily dose based on your phenotype."),
        ("Alternative drug",
         "A different medication in the same class is preferred for safety or efficacy reasons, given your genetic result."),
        ("Monitor / Other guidance",
         "No immediate dose change is required, but additional monitoring, specific laboratory tests, or counselling about risks is recommended."),
        ("No action required",
         "Your genetic result does not indicate a need to change the drug or dose, though your doctor still considers all other clinical factors."),
    ]
    rec_body = ""
    for i, (rtype, meaning) in enumerate(rec_rows):
        bg = "#f9fafb" if i % 2 == 0 else "#ffffff"
        rec_body += (f'<div class="intro-safe-row" style="display:grid; grid-template-columns:28% 72%; background:{bg};">'
                     f'<div style="padding:5px 8px; font-size:12px; font-weight:600; border-left:1px solid #e5e7eb; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#1a1a1a; box-sizing:border-box;">{rtype}</div>'
                     f'<div style="padding:5px 8px; font-size:12px; border-right:1px solid #e5e7eb; border-bottom:1px solid #e5e7eb; color:#374151; line-height:1.5; box-sizing:border-box;">{meaning}</div>'
                     f'</div>')

    # ========================================================================
    # PAGE 1
    # Metabolizer + transporter phenotypes only.
    # ========================================================================

    page_content1 = f"""
    <div style="padding-top:10px;">

        <div
            class="page-title"
            id="faqs"
        >
            How to Understand Your Results
        </div>

        {_h2("Understanding Metabolizer Phenotypes")}

        {_para(
            "Most pharmacogenomics results for drug-metabolizing enzyme "
            "genes are expressed as a metabolizer status. This describes "
            "how efficiently your enzyme breaks down a particular drug."
        )}

        <div class="intro-safe-grid" style="margin-bottom:10px;">
            <div class="intro-safe-grid-header" style="display:grid; grid-template-columns:24% 36% 40%; background:#4DB7D0;">
                <div style="padding:4px 6px; font-size:12px; font-weight:700; color:white;">Phenotype</div><div style="padding:4px 6px; font-size:12px; font-weight:700; color:white;">What It Means</div><div style="padding:4px 6px; font-size:12px; font-weight:700; color:white;">Common Clinical Implication</div>
            </div>{pheno_body}
        </div>

        {_h2("Understanding Transporter Function Phenotypes")}

        {_para(
            "Transporter genes such as SLCO1B1 and ABCG2 use a different "
            "phenotype vocabulary because they move drugs into and out "
            "of cells rather than metabolizing them."
        )}

        <div class="intro-safe-grid" style="margin-bottom:6px;">
            <div class="intro-safe-grid-header" style="display:grid; grid-template-columns:24% 36% 40%; background:#4DB7D0;">
                <div style="padding:4px 6px; font-size:12px; font-weight:700; color:white;">Phenotype</div><div style="padding:4px 6px; font-size:12px; font-weight:700; color:white;">What It Means</div><div style="padding:4px 6px; font-size:12px; font-weight:700; color:white;">Common Clinical Implication</div>
            </div>{transporter_body}
        </div>

    </div>
    """

    # ========================================================================
    # PAGE 2
    # SLCO1B1 example + susceptibility phenotypes.
    # ========================================================================

    page_content2 = f"""
    <div style="padding-top:10px;">

        {_continued_title()}

        {_h2(
            "Why Can the Same Genetic Result Lead to Different "
            "Recommendations for Different Drugs?"
        )}

        {_para(
            "<strong>Example: SLCO1B1 and statins.</strong> "
            "An SLCO1B1 Decreased Function result can increase plasma "
            "exposure to several statins. However, the clinical action "
            "is not identical for every statin because each medication "
            "has a different baseline exposure-response relationship "
            "and risk profile. For example, simvastatin and lovastatin "
            "may carry greater concern for muscle toxicity at increased "
            "exposure than some other statins. This is why the same "
            "genetic result can lead to different recommendations for "
            "different medications."
        )}

        {_h2("Understanding Susceptibility Phenotypes")}

        {_para(
            "Some genes in this report do not primarily affect how a drug "
            "is metabolized. Instead, certain variants can indicate "
            "susceptibility to a physiological or immune-mediated adverse "
            "reaction when specific medications are used."
        )}

        <div class="intro-safe-grid" style="margin-bottom:10px;">
            <div class="intro-safe-grid-header" style="display:grid; grid-template-columns:20% 10% 40% 30%; background:#4DB7D0;">
                <div style="padding:4px 6px; font-size:12px; font-weight:700; color:white;">Phenotype</div><div style="padding:4px 6px; font-size:12px; font-weight:700; color:white;">Gene(s)</div><div style="padding:4px 6px; font-size:12px; font-weight:700; color:white;">What It Means</div><div style="padding:4px 6px; font-size:12px; font-weight:700; color:white;">Clinical Implication</div>
            </div>{susc_body}
        </div>

    </div>
    """

    page_content3 = f"""
    <div style="padding-top:10px;">
        {_continued_title()}
        {_h2("Variant-Classified Genes (Not Metabolizer-Based)")}
        {_para("The following genes do not fit the metabolizer or transporter framework. Their results are reported as specific variant names or genotype calls rather than phenotype categories:")}
        <div class="intro-safe-grid" style="margin-bottom:10px;">
            <div class="intro-safe-grid-header" style="display:grid; grid-template-columns:15% 35% 50%; background:#4DB7D0;">
                <div style="padding:4px 6px; font-size:12px; font-weight:700; text-align:left; color:white; box-sizing:border-box;">Gene</div>
                <div style="padding:4px 6px; font-size:12px; font-weight:700; text-align:left; color:white; box-sizing:border-box;">Phenotype Format</div>
                <div style="padding:4px 6px; font-size:12px; font-weight:700; text-align:left; color:white; box-sizing:border-box;">What It Means</div>
            </div>
            {variant_body}
        </div>

        {_h2("Understanding the DPYD Activity Score")}
        {_para("DPYD is the primary enzyme responsible for breaking down fluoropyrimidine drugs (capecitabine, fluorouracil, tegafur, flucytosine). Unlike other genes in this report, DPYD phenotype is expressed using an <strong>Activity Score (AS)</strong> system that reflects the combined functional impact of both gene copies.")}
        <div class="intro-safe-grid" style="margin-bottom:10px;">
            <div class="intro-safe-grid-header" style="display:grid; grid-template-columns:15% 35% 50%; background:#4DB7D0;">
                <div style="padding:4px 6px; font-size:12px; font-weight:700; text-align:center; color:white;">Activity Score</div><div style="padding:4px 6px; font-size:12px; font-weight:700; color:white;">Phenotype</div><div style="padding:4px 6px; font-size:12px; font-weight:700; color:white;">Clinical Implication</div>
            </div>{dpyd_body}
        </div>

        {_h2("DPYD Diplotype Notation")}
        {_para("Unlike most genes in this report, which use <strong>star allele notation</strong> (e.g., *1/*5, *1/*9), DPYD diplotypes are reported in <strong>cDNA nucleotide notation</strong> (e.g., c.1057C>T/c.1484A>G). This is the established international convention for DPYD because the gene's variant landscape does not lend itself to the star allele system.")}
        {_para("Each 'c.' entry represents a specific nucleotide change at a defined position in the DPYD gene sequence. The diplotype shown lists one variant on each chromosome copy. You do not need to interpret the raw notation &mdash; the phenotype (e.g., Poor Metabolizer) and Activity Score (e.g., 0.0) are the clinically actionable result.")}
    </div>
    """

    page_content4 = f"""
    <div style="padding-top:10px;">
        {_continued_title()}
        <div style="page-break-inside: avoid; break-inside: avoid; margin-bottom: 12px;">
            {_h2("Understanding the PRODRUG Label")}
            {_para("Some drugs in this report are labeled <strong>PRODRUG</strong>. This label is clinically relevant for pharmacogenomics because the relationship between a metabolizer phenotype and clinical outcome differs:")}
            <ul style="font-size:11.5px; color:#374151; line-height:1.6; margin-top:4px; margin-bottom:8px; padding-left:18px;">
                <li>For a <strong>standard active drug</strong>, a Poor Metabolizer of the breakdown enzyme means the drug <strong>accumulates</strong> &rarr; higher toxicity risk.</li>
                <li>For a <strong>prodrug</strong> that must be converted to an active form, a Poor Metabolizer of the activating enzyme may not generate enough active drug &rarr; reduced efficacy.</li>
            </ul>
            {_para("In this report:")}
            <ul style="font-size:11.5px; color:#374151; line-height:1.6; margin-top:4px; margin-bottom:8px; padding-left:18px;">
                <li><strong>Capecitabine and Tegafur</strong> are prodrugs of 5-fluorouracil. They are converted to 5-FU in the body. The DPYD enzyme then breaks down 5-FU. DPYD Poor Metabolizer status means 5-FU accumulates to toxic levels.</li>
                <li><strong>5-Fluorouracil</strong> is the active cytotoxic agent itself (not a prodrug). It is the active metabolite that capecitabine and tegafur are converted into.</li>
                <li><strong>Azathioprine, Mercaptopurine, Thioguanine</strong> are converted to active thiopurine metabolites. TPMT and NUDT15 inactivate the toxic metabolites. A Poor Metabolizer here means the toxic metabolites are not cleared quickly enough.</li>
            </ul>
        </div>

        <div style="page-break-inside: avoid; break-inside: avoid; margin-bottom: 12px;">
            {_h2("What If a Gene Shows More Than One Diplotype or Phenotype?")}
            {_para("For some genes (such as SLCO1B1, NAT2, CYP4F2), your report may show multiple possible diplotypes in a single row. This happens when the analysis cannot definitively phase your variants &mdash; that is, it can identify the variants you carry but cannot always determine with certainty which variants sit on the same chromosome copy. Rather than exclude ambiguous results, the report shows all clinically plausible diplotype-phenotype combinations so your clinician has the full picture.")}
            {_para(
                "For <strong>CYP4F2 and warfarin specifically:</strong> "
                "When more than one plausible diplotype is shown, the "
                "possible interpretations may have similar or overlapping "
                "clinical implications. Warfarin dosing should not be "
                "interpreted from CYP4F2 alone; CYP2C9, VKORC1, available "
                "clinical factors, and the completeness of the submitted "
                "genetic data should all be considered together."
            )}
        </div>

        <div style="page-break-inside: avoid; break-inside: avoid; margin-bottom: 12px;">
            {_h2("What If Different Genes Point to Different Outcomes for the Same Drug?")}
            {_para("Some drugs are affected by more than one gene. This report presents the recommendation for each gene separately. Your clinician should consider all findings together. The combined effect of multiple gene variants can be additive &mdash; so the overall clinical picture may be more nuanced than any single result alone.")}
        </div>

        <div style="page-break-inside: avoid; break-inside: avoid; margin-bottom: 12px;">
            {_h2("How to Read the Clinical Recommendations")}
            {_para("Each drug section ends with a Clinical Recommendations block, sourced from CPIC, DPWG, ClinPGx, or the FDA. Here is how to read them:")}
            <div class="intro-safe-grid" style="margin-bottom:8px;">
                <div class="intro-safe-grid-header" style="display:grid; grid-template-columns:28% 72%; background:#4DB7D0;">
                    <div style="padding:4px 6px; font-size:12px; font-weight:700; color:white;">Recommendation Type</div><div style="padding:4px 6px; font-size:12px; font-weight:700; color:white;">What It Means in Practice</div>
                </div>{rec_body}
            </div>
        </div>
    </div>
    """

    page_content5 = f"""
    <div style="padding-top:10px;">
        {_continued_title()}
        {_h2("What Is the 'Other Evaluated Medications' Section?")}
        {_para("This section lists drugs that were assessed as part of the analysis, but for which there is no recognized guideline-level clinical recommendation. These are included because:")}
        <ul style="font-size:11.5px; color:#374151; line-height:1.6; margin-top:4px; margin-bottom:8px; padding-left:18px;">
            <li>Your result may become clinically relevant as evidence evolves.</li>
            <li>Your doctor may find the genotype or phenotype information useful even without a formal guideline recommendation.</li>
            <li>Transparency matters &mdash; knowing what was tested (and what was not) is important context.</li>
        </ul>
        {_para("These entries should not be used as the basis for prescribing changes on their own, but can be a useful reference for your healthcare team.")}

        {_h2("What If a Drug You Take Is Not Listed?")}
        {_para("If a medication you currently take does not appear anywhere in this report, it most likely means that either (a) there is no established pharmacogenomic guideline linking any of your tested genes to that drug, or (b) the relevant gene markers were not captured in your submitted data file.")}
        {_para("The absence of a listing does not mean genetics are irrelevant &mdash; it means the evidence has not yet reached the threshold for a guideline recommendation, or the drug falls outside the scope of this analysis.")}

        {_h2("What 'Not Evaluated by Platform' Means")}
        {_para("In some drug entries, you will see a gene listed with the notation <strong>'Not evaluated by platform'</strong> instead of a diplotype and phenotype. This means the gene's relevant variants were not captured in the DNA dataset submitted for analysis. This is a data limitation of consumer-grade SNP array genotyping, not a problem with the gene itself.")}
        {_para("Common reasons include: the gene's variants are located in regions not covered by standard SNP arrays (e.g., CYP2D6 structural variants), the gene requires copy-number variant detection, or the variant is rare. 'Not evaluated by platform' means <strong>no conclusion can be drawn</strong>. If the drug is relevant, a dedicated pharmacogenomic test from a clinical laboratory may be ordered.")}

        <div style="border:1px solid #4DB7D0; border-left:4px solid #4DB7D0; border-radius:6px; background:#f0fffe; padding:8px 14px; margin-bottom:8px; page-break-inside: avoid; break-inside: avoid;">
            <div style="font-size:11.5px; font-weight:700; color:#023D79; margin-bottom:3px;">Reminder</div>
            <div style="font-size:12px; line-height:1.6; color:#374151;">
                A strong evidence level (e.g., CPIC A, ClinPGx 1A) means the gene-drug relationship is well-established; it does not automatically mean you must change your medication. All prescribing decisions remain at the discretion of your qualified healthcare provider.
            </div>
        </div>
    </div>
    """

    combined_content = (
        page_content1
        + page_content2
        + page_content3
        + page_content4
        + page_content5
    )

    page = _wrap_page(
        combined_content,
        name,
        pg,
        page_id="faqs",
        background_color="#f3f4f6"
    )

    return [page], pg + 1

# ============================================================================
# 4. FOR YOUR DOCTOR PAGE
# ============================================================================

def doctor_page_template(df, master_genes_df, name, pg):
    _BAD = {"indeterminate", "no call", "no data available", "unknown", "no result",
            "nan", ""}

    def _pheno_significance(pheno_lower):
        if "ultrarapid" in pheno_lower or "ultra-rapid" in pheno_lower:
            return "Ultrarapid metabolism &mdash; standard dose may be sub-therapeutic"
        if "poor" in pheno_lower:
            return "Poor metabolism &mdash; drug accumulation risk; alternatives may be needed"
        if "intermediate" in pheno_lower:
            return "Intermediate metabolism &mdash; dose adjustment or monitoring may be needed"
        if "normal" in pheno_lower or "extensive" in pheno_lower:
            return "Normal metabolism &mdash; standard dosing applies"
        if "decreased function" in pheno_lower or "poor function" in pheno_lower:
            return "Reduced transporter/enzyme function &mdash; increased drug exposure possible"
        if "increased function" in pheno_lower:
            return "Increased transporter function &mdash; reduced drug exposure possible"
        if "malignant hyperthermia" in pheno_lower or "susceptibility" in pheno_lower:
            return "Susceptibility to adverse reaction &mdash; avoid triggering agents"
        if "deficient" in pheno_lower:
            return "Enzyme deficiency &mdash; risk of toxicity with substrates"
        return "Variant phenotype &mdash; clinical review recommended"

    def _section_hdr(title, color="#0097A7"):
        return (f'<div style="font-size:12px; font-weight:700; color:{color}; '
                f'margin-top:10px; margin-bottom:4px; padding-bottom:3px; '
                f'border-bottom:1.5px solid #e5e7eb;">{title}</div>')

    page_content1 = f"""
    <div style="padding-top: 10px;">
        <div class="page-title" id="doctor">Note to Doctor</div>

        {_section_hdr("Report Scope")}
        <div style="font-size:11.5px; line-height:1.6; color:#374151; margin-bottom:8px;">
            This pharmacogenomics report was generated using SNP array-based genotyping for <strong>{name}</strong>.
            Results reflect pharmacogenomically relevant genetic variants for the genes listed in the Genotype Summary.
            Recommendations are derived from CPIC, DPWG, and FDA prescribing information, as applicable. PharmCAT v3.2.0 was used for variant calling and phenotype assignment.
            Only gene-drug pairs with a defined phenotype and a recognized guideline-linked recommendation are surfaced as clinical recommendations.
        </div>

        {_section_hdr("Evidence Sources")}
        <div style="font-size:11.5px; line-height:1.6; color:#374151; margin-bottom:8px;">
            CPIC (cpicpgx.org) and DPWG guidelines are peer-reviewed and regularly updated by expert working groups. DPWG guidelines are produced by Dutch clinical pharmacogenomics experts and are widely used when CPIC has not yet issued a guideline for a specific gene-drug pair. ClinPGx (formerly PharmGKB, clinpgx.org) curates gene-drug association evidence that underpins guideline development. FDA label annotations reflect current approved prescribing information. Guideline classifications (Strong, Moderate, Optional) are included per source where available.
        </div>

        {_section_hdr("How to Identify Key Findings")}
        <div style="font-size:11.5px; line-height:1.6; color:#374151; margin-bottom:8px;">
            To identify your most critical findings, please refer to the <strong>Summary of Medication Insights</strong> on the following page. That summary provides a prioritized snapshot of all medications evaluated, categorized by the level of clinical action required (e.g., Action Required, Use With Caution).
        </div>

        {_section_hdr("Guidance on Missing Genes for Key Drug Decisions")}
        <div style="font-size:11.5px; line-height:1.6; color:#374151; margin-bottom:8px;">
            <p style="margin-bottom:4px;"><strong>CYP2C9 (not callable from this dataset)</strong><br>CYP2C9 affects multiple clinically important drugs including warfarin, several NSAIDs, glipizide, glyburide, and losartan. For warfarin in particular, CYP2C9 is the dominant metabolizing enzyme and the primary variable in all validated pharmacogenetic dosing algorithms (including the CPIC-recommended WarfarinDosing.org calculator). This report\'s warfarin findings are incomplete without CYP2C9 data. A dedicated CYP2C9 pharmacogenomic assay from a clinical laboratory is strongly recommended before attempting genotype-guided warfarin dosing.</p>
            <p style="margin-bottom:4px;"><strong>TPMT (not callable from this dataset)</strong><br>TPMT and NUDT15 both independently affect thiopurine toxicity (azathioprine, mercaptopurine, thioguanine). Without TPMT data, the risk assessment is incomplete. In patients of European descent, TPMT has historically been the primary predictor of thiopurine myelosuppression; in East Asian populations, NUDT15 is generally more informative. Before initiating thiopurine therapy, consider ordering a dedicated TPMT enzyme activity assay (available from most clinical biochemistry laboratories) or a clinical-grade TPMT genotyping panel. Monitor CBC closely during dose escalation regardless.</p>
            <p style="margin-bottom:4px;"><strong>CYP2D6 (not callable from this dataset)</strong><br>CYP2D6 affects over 25% of commonly prescribed drugs, including codeine, tramadol, tamoxifen, many antidepressants (paroxetine, fluoxetine, nortriptyline, venlafaxine), and many antipsychotics. For patients requiring any of these medications &mdash; particularly codeine or tramadol (where CYP2D6 poor or ultrarapid metabolizer status carries serious safety implications), or tamoxifen (where CYP2D6 status affects endoxifen generation and treatment efficacy) &mdash; a dedicated clinical CYP2D6 genotyping panel using long-read sequencing or a validated clinical assay is recommended.</p>
            <p style="margin-bottom:4px;"><strong>HLA-B (not callable from this dataset)</strong><br>HLA-B allele-level typing is required before prescribing abacavir (HLA-B*57:01), allopurinol (HLA-B*58:01), and certain anticonvulsants (HLA-B*15:02) to prevent severe cutaneous adverse reactions. The risk of HLA-B*58:01-associated allopurinol SJS is markedly higher in individuals of Han Chinese, Korean, Thai, and Vietnamese descent. Dedicated high-resolution HLA typing from a clinical immunogenetics laboratory is required.</p>
        </div>

        {_section_hdr("Technical Interpretation Notes")}
        <div style="font-size:11.5px; line-height:1.6; color:#374151; margin-bottom:8px;">
            <ul style="margin:0; padding-left:16px;">
                <li>Diplotype calling performed using PharmCAT v3.2.0 with standard reference data.</li>
                <li>SNP array coverage may not detect all star alleles; rare or novel variants may be absent.</li>
                <li>Copy number variants (e.g., CYP2D6 gene duplications) require dedicated confirmatory testing.</li>
                <li>HLA alleles (HLA-A, HLA-B) require specialized high-resolution typing and are not included here.</li>
                <li>Unphased results (e.g., CYP4F2) indicate phase-ambiguous variants; all clinically plausible interpretations are reported for completeness.</li>
                <li>DPYD diplotype is reported in cDNA notation (c.XXXXN&gt;N) per international convention, not star allele notation.</li>
                <li>The genotypic phenotypes in this report represent predicted drug metabolism based on germline DNA and do not account for phenoconversion due to concomitant medications, renal or hepatic impairment, or other clinical factors.</li>
            </ul>
        </div>

        {_section_hdr("A Note on Phenoconversion")}
        <div style="font-size:11.5px; line-height:1.6; color:#374151; margin-bottom:8px;">
            <p style="margin-bottom:4px;">A patient\'s genetically predicted phenotype may not reflect their actual functional phenotype if they are concurrently taking medications that inhibit or induce the same enzyme. <strong>Example:</strong> A patient who is a genotypic Normal Metabolizer for CYP2D6 may function as a Poor Metabolizer if they are also taking a strong CYP2D6 inhibitor such as fluoxetine, paroxetine, or bupropion.</p>
            <p style="margin-bottom:4px;">Clinicians should always cross-check the patient\'s medication list for known inhibitors and inducers of any gene flagged in this report. Common phenoconverting agents include:</p>
            <div class="intro-safe-grid" style="margin-bottom:8px;">
                <div class="intro-safe-grid-header" style="display:grid; grid-template-columns:15% 45% 40%; background:#4DB7D0;">
                    <div style="padding:4px 6px; font-size:11px; font-weight:700; color:white;">Gene</div>
                    <div style="padding:4px 6px; font-size:11px; font-weight:700; color:white;">Common Strong Inhibitors (may convert to Poor Metabolizer)</div>
                    <div style="padding:4px 6px; font-size:11px; font-weight:700; color:white;">Common Inducers (may increase metabolism)</div>
                </div>
                <div class="intro-safe-row" style="display:grid; grid-template-columns:15% 45% 40%; background:#f9fafb;"><div style="padding:4px 6px; font-size:11px; font-weight:600; border:1px solid #e5e7eb;">CYP2D6</div><div style="padding:4px 6px; font-size:11px; border:1px solid #e5e7eb;">Fluoxetine, paroxetine, bupropion, terbinafine, quinidine</div><div style="padding:4px 6px; font-size:11px; border:1px solid #e5e7eb;">None</div></div>
                <div class="intro-safe-row" style="display:grid; grid-template-columns:15% 45% 40%; background:#ffffff;"><div style="padding:4px 6px; font-size:11px; font-weight:600; border:1px solid #e5e7eb;">CYP2C19</div><div style="padding:4px 6px; font-size:11px; border:1px solid #e5e7eb;">Fluvoxamine, fluconazole, ticlopidine</div><div style="padding:4px 6px; font-size:11px; border:1px solid #e5e7eb;">Rifampin, carbamazepine</div></div>
                <div class="intro-safe-row" style="display:grid; grid-template-columns:15% 45% 40%; background:#f9fafb;"><div style="padding:4px 6px; font-size:11px; font-weight:600; border:1px solid #e5e7eb;">CYP3A4/5</div><div style="padding:4px 6px; font-size:11px; border:1px solid #e5e7eb;">Ketoconazole, itraconazole, clarithromycin, grapefruit juice</div><div style="padding:4px 6px; font-size:11px; border:1px solid #e5e7eb;">Rifampin, carbamazepine, phenytoin, St. John's wort</div></div>
                <div class="intro-safe-row" style="display:grid; grid-template-columns:15% 45% 40%; background:#ffffff;"><div style="padding:4px 6px; font-size:11px; font-weight:600; border:1px solid #e5e7eb;">CYP2C9</div><div style="padding:4px 6px; font-size:11px; border:1px solid #e5e7eb;">Fluconazole, amiodarone, miconazole</div><div style="padding:4px 6px; font-size:11px; border:1px solid #e5e7eb;">Rifampin</div></div>
            </div>
            <p style="margin-bottom:4px;">This report does not perform phenoconversion analysis. That assessment must be performed by the prescribing clinician.</p>
        </div>

        {_section_hdr("Appropriate Use")}
        <div style="font-size:11.5px; line-height:1.6; color:#374151;">
            <strong>This report is intended to support, not replace, clinical judgment.</strong>
            Decisions about prescribing or discontinuing medications should integrate this report with the patient\'s
            full clinical picture &mdash; including comorbidities, concomitant medications, organ function, and patient preferences. The pharmacogenomic findings presented here are one input among many in the prescribing decision. This report does not constitute a prescription, a contraindication, or a guarantee of drug response.
        </div>
    </div>
    """

    pg1 = _wrap_page(page_content1, name, pg, page_id="doctor", background_color="#f3f4f6")
    
    return [pg1], pg + 1
# # ============================================================================
# 5. EXECUTIVE SUMMARY PAGE
# ============================================================================

def executive_summary_template(df, patient_name, pg, zone_df=None):
    """
    5-Bucket Executive Summary ONLY (Action Required → No Guideline)
    """
    if zone_df is None or zone_df.empty:
        zone_df = df.copy()

    bucket_map = {
        "action_required":   ("Action Required",                    "#fee2e2", "#dc2626"),
        "monitoring":        ("Use With Caution / Monitoring",      "#fef3c7", "#d97706"),
        "standard_use":      ("Standard Use",                       "#ecfdf5", "#16a34a"),
        "further_testing":   ("Further Testing Required",           "#f3e8ff", "#8b5cf6"),
        "no_guideline":      ("No Guideline Available",              "#f3f4f6", "#6b7280"),
    }

    bucket_priority = {
        "further_testing": 1,
        "action_required": 2,
        "monitoring": 3,
        "standard_use": 4,
        "no_guideline": 5,
    }

    def _normalize_bucket(value: str) -> str:
        bucket = str(value or "").strip().lower().replace(" / ", "_")
        bucket = bucket.replace(" ", "_").replace("-", "_")
        if bucket in bucket_priority:
            return bucket
        if bucket in {"use_with_caution", "use_with_caution_monitoring", "monitoring"}:
            return "monitoring"
        if bucket in {"action_required", "actionrequired"}:
            return "action_required"
        if bucket in {"further_testing", "further_testing_required", "furthertesting"}:
            return "further_testing"
        if bucket in {"standard_use", "standarduse"}:
            return "standard_use"
        if bucket in {"no_guideline", "noguideline", "no guidance"}:
            return "no_guideline"
        return "no_guideline"

    zone = zone_df.copy()
    if "Summary Bucket" not in zone.columns:
        if "Best Status" in zone.columns:
            zone["Summary Bucket"] = zone["Best Status"].astype(str).map(lambda v: _normalize_bucket(v))
        else:
            zone["Summary Bucket"] = "no_guideline"

    zone["Drug Name"] = zone["Drug Name"].astype(str).str.strip()
    zone["Summary Bucket"] = zone["Summary Bucket"].astype(str).map(_normalize_bucket)
    if "Gene" in zone.columns:
        zone["Gene"] = zone["Gene"].astype(str).str.strip()
    else:
        zone["Gene"] = ""

    def _genes_for_group(group_df):
        genes = [g for g in group_df["Gene"].tolist() if g and g.lower() not in ("nan", "none")]
        return ", ".join(dict.fromkeys(genes))

    cat_col_name = "Therapeutic Category" if "Therapeutic Category" in zone.columns else "Drug Category"
    
    drug_records = []
    for drug_name, grp in zone.groupby("Drug Name", sort=False):
        buckets = [b for b in grp["Summary Bucket"].tolist() if b]
        worst_bucket = min(buckets, key=lambda b: bucket_priority.get(b, 99)) if buckets else "no_guideline"
        
        # Get the first category for the drug to build the anchor
        cat_val = ""
        if cat_col_name in grp.columns:
            valid_cats = grp[cat_col_name].dropna().astype(str).tolist()
            if valid_cats:
                cat_val = valid_cats[0]

        drug_records.append({
            "drug": drug_name,
            "drug_display": drug_name.title(),
            "genes": _genes_for_group(grp),
            "bucket": worst_bucket,
            "category": cat_val,
        })

    summary = {}
    for record in drug_records:
        summary[record["bucket"]] = summary.get(record["bucket"], 0) + 1

    html = f"""
    <div style="padding-top: 10px;">
        <div class="page-title">Summary of Medication Insights</div>
        <p style="font-size:12px; color:#374151; margin-bottom:16px;">
            Based on your genetic profile, <strong>{len(drug_records)} medications</strong> were evaluated.
        </p>
    """

    for bucket_key in ["action_required", "monitoring", "further_testing", "standard_use", "no_guideline"]:
        title, bg_color, accent_color = bucket_map[bucket_key]
        count = summary.get(bucket_key, 0)
        if count == 0:
            continue

        drugs_html = ""
        for row in sorted((r for r in drug_records if r["bucket"] == bucket_key), key=lambda x: x["drug_display"]):
            drug = row["drug_display"]
            genes = row["genes"]
            cat_val = row.get("category", "")
            
            # "No Guideline Available" drugs link to the no_guideline_available page
            if bucket_key == "no_guideline":
                link_target = "no_guideline_available"
            else:
                drug_id = sanitize_id(drug)
                link_target = drug_id + "__" + sanitize_id(cat_val) if cat_val else drug_id
            drugs_html += f"""
            <div style="padding:4px 0; border-bottom:1px solid #e5e7eb; font-size:11.5px; width: calc(50% - 6px); display: inline-block; box-sizing: border-box;">
                <a href="#{link_target}" style="color:#1a1a1a; text-decoration:none;">
                    <strong>{drug}</strong>
                </a>
                <span style="color:#6b7280; font-size:12px; margin-left:8px;">{genes}</span>
            </div>"""

        html += f"""
        <div style="margin-bottom:14px; border:none; border-left:4px solid {accent_color};
                    border-radius:8px; background:{bg_color}; padding:12px 14px; page-break-inside: auto;">
            <div style="font-size:11.5px; font-weight:700; color:{accent_color}; margin-bottom:8px;">
                {title} ({count} medications)
            </div>
            <div style="display: flex; flex-wrap: wrap; gap: 0 12px;">{drugs_html}</div>
        </div>"""

    html += """
        <div style="font-size:12px; color:#6b7280; margin-top:8px;">
            * All prescribing decisions must be made by your healthcare provider.
        </div>
    </div>
    """
    return _wrap_page(html, patient_name, pg, page_id="executive_summary")

    # ── Option A (legacy): keyword-based 3-bucket summary ──────────────────────
    CLF_ORDER = {"strong": 0, "moderate": 1, "optional": 2, "unspecified": 3, "no recommendation": 4}

    _SECTION_SPLIT_RE = re.compile(
        r'\b(Other\s+Considerations|Implications?|Affected\s+Subgroup|Submitted\s+Genotype)\b',
        re.IGNORECASE,
    )
    _RECOMMENDATION_LEAD_RE = re.compile(
        r'^\s*(?:.*?\bRecommendation\s*:?\s+)?(.+)$',
        re.DOTALL,
    )

    def _patient_rec_only(rec_text: str) -> str:
        if not rec_text:
            return ""
        head = _SECTION_SPLIT_RE.split(rec_text, maxsplit=1)[0]
        m = re.match(
            r'^\s*Recommendation\s*:?\s+(.+)$',
            head, flags=re.DOTALL,
        )
        return (m.group(1) if m else head).strip()

    # ── No-action detection ────────────────────────────────────────────────────
    # Ordered longest-first to avoid partial-match shadowing.
    NO_ACTION_PHRASES = [
        # Explicit no-action statements
        "no specific action", "no action required", "no action is needed",
        "no change needed", "no change is needed",
        "no dose adjustment", "no dose adjustment is anticipated",
        "no need to avoid", "no indication to", "no reason to avoid",
        "no alternative needed", "no recommendation", "no adjustment",
        "no clinical action", "no clinically significant",
        # Standard / recommended starting dose — these are PGx "go ahead" phrases
        "prescribe desired starting dose",
        "start with normal starting dose",    # azathioprine/MP/thioguanine NM
        "normal starting dose",               # same; substring match catches it
        "standard recommended dose", "standard recommended starting",
        "standard starting daily dose", "standard starting dose",
        "recommended starting dose", "recommended starting dosing",
        "initiate therapy with recommended", "initiate therapy with standard",
        "initiate standard dose", "initiate standard",
        "initiate therapy at standard dose",
        "use standard dose", "prescribe standard",
        "standard dose recommended",
        # Routine titration phrases — "adjust dose based on disease guidelines"
        # is NOT a PGx-specific dose change; it means follow normal clinical practice
        "adjust doses based on disease",
        "adjust doses based on clinical",
        "adjust doses based on specific population",
        "based on disease-specific guidelines",
        "based on disease-specific and specific population",
        # Standard-care TDM phrases — TDM is routine for many drugs regardless of PGx
        "use therapeutic drug monitoring to guide",
        "therapeutic drug monitoring to guide dose",
        "subsequent doses should be adjusted according to therapeutic",
        # Other standard-use confirmations
        "standard prescribing", "standard dosing",
        "standard of care dosing",
        "routine dosing algorithms",
        "no dose adjustment is anticipated",
    ]

    # ── Action-required keywords — specific phrases only ───────────────────────
    # Single words like "alternative" or "adjust" are too broad and match routine
    # clinical guidance.  Use multi-word or highly specific phrases.
    ACTION_REQUIRED_KEYWORDS = [
        "contraindicated",
        "do not use",
        "do not prescribe",
        "avoid this",  # "avoid this drug" / "avoid this agent" — not bare "avoid"
        "avoid using",
        "consider an alternative",
        "use an alternative",
        "use alternative",
        "select an alternative",
        "switch to an alternative",
        "alternative agent",
        "alternative drug",
        "dose reduction",
        "reduce dose", "reduce the dose",
        "decrease dose", "decreased dose",
        "decreased starting dose",
        "lower starting dose", "reduce starting dose",
        "not recommended",
        # Ultrarapid metabolizer — needs higher dose
        "increase the dose",
        "inadequate drug efficacy",
    ]

    # ── Use-with-caution keywords — specific phrases to avoid false positives ──
    # Bare "monitor" / "monitoring" match routine TDM instructions that appear in
    # NO_ACTION phrases above.  Bare "may increase" matches dose-escalation notes.
    CAUTION_KEYWORDS = [
        "closely monitor", "close monitoring",
        "more frequent monitoring", "additional monitoring",
        "enhanced monitoring", "increased monitoring",
        "monitor for adverse", "monitor for toxicity",
        "monitor for side effect",
        "use with caution", "prescribe with caution", "exercise caution",
        "proceed with caution",
        "increased risk of adverse", "increased risk of toxicity",
        "increased risk of side effect",
        "at increased risk of",
        "greater risk of", "higher risk of",
        "may increase the risk", "may increase toxicity",
        "may increase plasma", "may increase serum",
        "may decrease efficacy", "may decrease plasma",
        "lower dose", "start at a lower dose",
        "slower titration", "slower dose titration",
    ]

    # Trim "Other Considerations" and similar sections so keyword checks only
    # see the patient-specific guidance, not disease/population caveats.
    _EXEC_SPLIT_RE = re.compile(
        r'\s*\b(?:Other Considerations|Implications|Affected Subgroup|'
        r'Submitted Genotype|Recommendation:)\b',
        re.IGNORECASE,
    )
    def _exec_trim(text: str) -> str:
        head = _EXEC_SPLIT_RE.split(text, maxsplit=1)[0]
        return re.sub(r'^\s*Recommendation\s*:?\s+', '', head, flags=re.IGNORECASE).strip()

    def classify_action(rec_text, clf_text):
        rec = _exec_trim(rec_text).lower()
        clf = clf_text.lower()
        if any(phrase in rec for phrase in NO_ACTION_PHRASES):
            return "Standard Use"
        if any(k in rec for k in ACTION_REQUIRED_KEYWORDS):
            return "Action Required"
        # "avoid" alone is checked here (after no-action so "no need to avoid" is safe)
        if "avoid" in rec and "no need to avoid" not in rec and "need not avoid" not in rec:
            return "Action Required"
        if clf == "strong":
            return "Action Required"
        if any(k in rec for k in CAUTION_KEYWORDS):
            return "Use With Caution"
        if clf == "moderate":
            return "Use With Caution"
        return "Standard Use"

    def _action_label(rec_text, clf_text, pheno_text=""):
        rec = _exec_trim(rec_text).lower()
        clf = clf_text.lower()
        pheno = pheno_text.lower()

        # ── No action ──────────────────────────────────────────────────────────
        if any(phrase in rec for phrase in NO_ACTION_PHRASES):
            if "normal metabolizer" in pheno or "normal function" in pheno or pheno in ("normal",):
                return "Normal metabolism — standard dosing"
            if "intermediate metabolizer" in pheno:
                return "Intermediate metabolism — standard dosing"
            if "poor function" in pheno and "no action" in rec:
                return "Reduced function — standard dosing"
            if "indeterminate" in pheno or "no result" in pheno or "no data" in pheno:
                return "Genotype undetermined"
            if "no recommendation" in clf:
                return "No PGx guidance for your genotype"
            return "Standard prescribing applies"

        # ── Action required ───────────────────────────────────────────────────
        if "contraindicated" in rec:
            return "Contraindicated"
        if "do not use" in rec or ("avoid" in rec and "no need to avoid" not in rec):
            return "Avoid — use alternative"
        if "consider an alternative" in rec or "use an alternative" in rec:
            return "Consider alternative drug"
        if "not recommended" in rec:
            return "Not recommended"
        if "dose reduction" in rec or "reduce dose" in rec or "reduce the dose" in rec:
            return "Dose reduction recommended"
        if "alternative" in rec:
            return "Alternative drug recommended"

        # ── Use with caution ──────────────────────────────────────────────────
        if "closely monitor" in rec or "close monitoring" in rec:
            return "Monitor closely"
        if "increased monitoring" in rec or "additional monitoring" in rec:
            return "Increased monitoring recommended"
        if "caution" in rec:
            return "Use with caution"
        if "lower dose" in rec or "start at a lower dose" in rec:
            return "Start with lower dose"
        if "increased risk" in rec or "greater risk" in rec:
            return "Increased risk — exercise caution"

        # ── Fallback by classification ─────────────────────────────────────────
        if clf == "strong":
            return "Strong clinical guidance — see detail page"
        if clf == "moderate":
            return "Moderate clinical guidance — see detail page"
        return "Standard prescribing applies"

    work = df.copy()
    if "Classification" not in work.columns:
        work["Classification"] = "Unspecified"
    if "Recommendation" not in work.columns:
        work["Recommendation"] = ""

    has_rec = work["Classification"].astype(str).str.lower() != "no recommendation"
    has_text = work["Recommendation"].astype(str).str.strip().str.len() > 0
    work = work[has_rec | has_text].copy()

    work["_clf_rank"] = work["Classification"].astype(str).str.lower().map(
        lambda x: CLF_ORDER.get(x, 99)
    )
    gene_col = "Gene" if "Gene" in work.columns else None

    if "Drug Category" in df.columns:
        _detail_drug_names = set(
            df[df["Drug Category"].astype(str).str.strip() != "Other Evaluated Medications"]
            ["Drug Name"].astype(str).str.strip().str.lower()
        )
    else:
        _detail_drug_names = set()

    drug_records = []
    for drug_name, grp in work.groupby("Drug Name"):
        best_row = grp.sort_values("_clf_rank").iloc[0]
        rec_text  = str(best_row.get("Recommendation", "")).strip()
        clf_text  = str(best_row.get("Classification", "")).strip()
        pheno_text= str(best_row.get("Phenotype", "")).strip()
        gene_txt  = str(best_row.get(gene_col, "")) if gene_col else ""
        bucket    = classify_action(rec_text, clf_text)
        label     = _action_label(rec_text, clf_text, pheno_text)
        has_page  = str(drug_name).strip().lower() in _detail_drug_names

        drug_records.append({
            "drug":     drug_name.title(),
            "gene":     gene_txt,
            "label":    label,
            "bucket":   bucket,
            "has_page": has_page,
        })

    buckets = {
        "Action Required": [],
        "Use With Caution": [],
        "Standard Use": [],
    }
    for rec in drug_records:
        buckets[rec["bucket"]].append(rec)

    total_drugs = len(drug_records)

    def _drug_link(record):
        drug = record["drug"]
        anchor = sanitize_id(drug)
        target = f"#{anchor}" if record.get("has_page") else "#other_evaluated"
        return (
            f'<a href="{target}" style="color:#1a1a1a; text-decoration:none; '
            f'border-bottom:1px dotted #94a3b8;">'
            f'<strong>{drug}</strong></a>'
        )

    def _bucket_html(title, icon, border_color, bg_color, title_color, records):
        n = len(records)

        if n == 0:
            body = (
                '<div style="font-size:11.5px; color:#9ca3af; font-style:italic;">'
                'No medications in this category for your genetic profile.'
                '</div>'
            )

        elif title == "Standard Use":
            lines = ""
            for r in sorted(records, key=lambda x: x["drug"]):
                gene_badge = (
                    f'<span style="background:#e0f2fe; color:#0369a1; '
                    f'font-size:11.5px; font-weight:700; padding:1px 4px; '
                    f'border-radius:3px; margin-left:4px;">'
                    f'{r["gene"]}</span>'
                ) if r["gene"] else ""

                lines += (
                    f'<div style="padding:3px 0; '
                    f'border-bottom:1px solid #e5e7eb; '
                    f'font-size:12px; color:#374151; '
                    f'display:flex; justify-content:space-between;">'
                    f'<span>{_drug_link(r)}{gene_badge}</span>'
                    f'<span style="color:#6b7280; '
                    f'font-style:italic;">{r["label"]}</span>'
                    f'</div>'
                )
            body = f'<div style="margin-top:4px;">{lines}</div>'

        else:
            lines = ""
            for r in sorted(records, key=lambda x: x["drug"]):
                gene_badge = (
                    f'<span style="background:#e0f2fe; color:#0369a1; '
                    f'font-size:11.5px; font-weight:700; '
                    f'padding:1px 4px; border-radius:3px; '
                    f'margin-left:4px;">'
                    f'{r["gene"]}</span>'
                ) if r["gene"] else ""

                lines += (
                    f'<div style="padding:3px 0; '
                    f'border-bottom:1px solid rgba(0,0,0,0.06); '
                    f'font-size:12px; display:flex; '
                    f'justify-content:space-between; '
                    f'align-items:center;">'
                    f'<span>{_drug_link(r)}{gene_badge}</span>'
                    f'<span style="color:{title_color}; '
                    f'font-weight:600; font-size:12px; '
                    f'font-style:italic;">{r["label"]}</span>'
                    f'</div>'
                )
            body = f'<div style="margin-top:4px;">{lines}</div>'

        # Changed page-break-inside to 'auto' so long lists safely split across pages
        return f"""
        <div style=""
            margin-bottom:12px;
            border:none;
            border-left:4px solid {border_color};
            border-radius:8px;
            background:{bg_color};
            padding:10px 14px;
            page-break-inside: auto; 
            break-inside: auto;
        ">
            <div style=""
                font-size:12px;
                font-weight:700;
                color:{title_color};
                margin-bottom:6px;
                page-break-after: avoid;
                break-after: avoid;
            ">
                {title}
                <span style="
                    font-size:11.5px;
                    font-weight:500;
                    color:#6b7280;
                ">
                    ({n} medication{'s' if n != 1 else ''})
                </span>
            </div>

            {body}
        </div>
        """

    # --- THIS IS THE SECTION THAT WAS UN-INDENTED. IT IS NOW FIXED ---
    ar_html  = _bucket_html("Action Required",  "●", "#dc2626", "#fef2f2", "#dc2626", buckets["Action Required"])
    uwc_html = _bucket_html("Use With Caution", "●", "#d97706", "#fffbeb", "#92400e", buckets["Use With Caution"])
    su_html  = _bucket_html("Standard Use",     "●", "#16a34a", "#f0fdf4", "#15803d", buckets["Standard Use"])

    content = f"""
    <div style="padding-top: 10px;">
        <div class="page-title" id="executive_summary">Summary of Medication Insights</div>
        <p style="font-size:12px; color:#444; line-height:1.6; margin-bottom:14px;">
            Based on your genetic profile, <strong>{total_drugs} medication{'s' if total_drugs != 1 else ''}</strong>
            with pharmacogenomic guidance were identified. The findings below are grouped by the level of
            action they may require. Refer to the Drug-Specific Recommendations section for full details.
        </p>
        {ar_html}
        {uwc_html}
        {su_html}
        <div style="font-size:12px; color:#9ca3af; margin-top:8px; line-height:1.5;">
            * Classification is based on recommendation text and guideline evidence level. All decisions should be
            reviewed by your healthcare provider in the context of your full clinical history.
        </div>
    </div>
    """
    return _wrap_page(content, patient_name, pg, page_id="executive_summary")

# ============================================================================
# COVERAGE STATEMENT PAGE
# ============================================================================

def coverage_statement_template(patient_name, pg, present_categories):
    all_categories = [
        "MUSCULO-SKELETAL SYSTEM DRUGS",
        "CARDIOVASCULAR SYSTEM DRUGS",
        "ALIMENTARY TRACT AND METABOLISM DRUGS",
        "NERVOUS SYSTEM DRUGS",
        "RESPIRATORY SYSTEM DRUGS",
        "ANTIINFECTIVES FOR SYSTEMIC USE",
        "GENITO URINARY SYSTEM AND SEX HORMONES",
        "BLOOD AND BLOOD FORMING ORGAN DRUGS",
        "DERMATOLOGICALS",
        "ANTINEOPLASTIC AND IMMUNOMODULATING AGENTS",
        "SENSORY ORGAN DRUGS",
        "VARIOUS DRUG CLASSES IN ATC",
        "ANTIPARASITIC PRODUCTS, INSECTICIDES AND REPELLENTS",
        "SYSTEMIC HORMONAL PREPARATIONS, EXCL. SEX HORMONES AND INSULINS",
        "NO GROUP ASSIGNED"
    ]

    # Clean the present_categories for robust case-insensitive matching
    covered_set = {c.strip().lower() for c in present_categories}

    rows_html = ""
    for cat in all_categories:
        is_covered = cat.lower() in covered_set
        
        # Format the category name in Title Case and apply manual corrections
        display_cat = cat.title()
        display_cat = display_cat.replace("In Atc", "In ATC")
        display_cat = display_cat.replace("Genito Urinary", "Genitourinary")
        display_cat = display_cat.replace("Antiinfectives", "Anti-infectives")

        if is_covered:
            # Green checkmark and bold text
            icon = '<span style="color:#10b981; font-weight:bold; font-size:11.5px;">✓</span>'
            text_style = 'color:#111827; font-weight:700;'
        else:
            # Subdued grey text
            icon = '<span style="color:#d1d5db; font-size:11.5px;">—</span>'
            text_style = 'color:#9ca3af; font-weight:400;'

        rows_html += f"""
        <tr style="border-bottom:1px solid #f3f4f6;">
            <td style="padding:10px 14px; text-align:center; width:80px;">{icon}</td>
            <td style="padding:10px 14px; {text_style} font-size:12px; letter-spacing:0.3px;">{display_cat}</td>
        </tr>
        """

    content = f"""
    <div style="padding-top:10px;">
        <div class="page-title" id="coverage_statement" style="margin-bottom: 12px;">Coverage Statement</div>
        <p style="font-size:11.5px; color:#4b5563; line-height:1.6; margin-bottom:18px;">
            This report evaluates medications across the following therapeutic categories.
            Categories with a <span style="color:#10b981; font-weight:bold;">✓</span> indicate that relevant medications with pharmacogenomic guidance were identified in your profile.
        </p>

        <table style="width:100%; border-collapse:collapse; border:1px solid #dde2ea; background:#ffffff;">
            <thead style="display: table-header-group; background:#0D3B7A; border-bottom:1px solid #dde2ea;">
                <tr>
                    <th style="padding:8px 14px; text-align:center; font-size:12px; font-weight:700; color:#ffffff; width:80px;">Status</th>
                    <th style="padding:8px 14px; text-align:left; font-size:12px; font-weight:700; color:#ffffff;">Therapeutic Category</th>
                </tr>
            </thead>
            <tbody>
                {rows_html}
            </tbody>
        </table>

        <div style="font-size:11.5px; color:#6b7280; margin-top:14px; line-height:1.5;">
            * All medications with pharmacogenomic relevance are shown. Actions are indicated where applicable in the subsequent pages.<br><br>
            A dash (&mdash;) indicates that no medications in this therapeutic category had an established pharmacogenomic guideline result based on your tested genes. It does not indicate that the category was outside the scope of analysis &mdash; all 15 categories listed above are evaluated for every report.
        </div>
    </div>
    """
    return _wrap_page(content, patient_name, pg, page_id="coverage_statement")

# ============================================================================
# 6. DISCLAIMER PAGE (unlisted — appended after Genotype Summary)
# ============================================================================

def disclaimer_template(patient_name, pg):
    def _h2(text):
        return (f'<div style="font-size:12px; font-weight:700; color:#023D79; '
                f'margin-top:10px; margin-bottom:4px; padding-bottom:3px; '
                f'border-bottom:1.5px solid #e5e7eb;">{text}</div>')

    def _para(text):
        return f'<p style="font-size:11.5px; color:#374151; line-height:1.6; margin-bottom:6px;">{text}</p>'

    page_content = f"""
    <div style="padding-top: 10px;">
        <div class="page-title" id="disclaimer">Disclaimer</div>

        {_h2("Purpose and Limitations of This Report")}
        {_para("This pharmacogenomics report is provided for informational purposes only. It is intended to be used as a supplementary tool to support, not replace, clinical decision-making by qualified healthcare professionals. It does not constitute medical advice, a diagnosis, or a treatment recommendation.")}

        {_h2("Not a Substitute for Medical Care")}
        {_para("The findings in this report must be interpreted in the full context of your individual health situation by a licensed physician, pharmacist, or other qualified healthcare provider. Genetic information is one factor among many that influences medication safety and efficacy. Other critical factors include your age, weight, kidney and liver function, other medications you take (including over-the-counter drugs and supplements), and underlying medical conditions.")}
        {_para("Never start, stop, or change any medication based solely on the results of this report. Always consult your doctor or pharmacist before making any change to your treatment.")}

        {_h2("Scope of Analysis")}
        {_para("This report covers only the gene-drug combinations for which established pharmacogenomic guidelines exist at the time of analysis. It does not capture all possible genetic influences on drug response, and it does not assess genes or variants outside the scope of the submitted data file and analysis pipeline. The absence of a finding does not mean a drug is safe or appropriate for you &mdash; it means genetic guidance is not yet available for that specific combination.")}

        {_h2("Evolving Evidence")}
        {_para("Pharmacogenomics is a rapidly developing field. Guidelines and evidence levels change as new research emerges. The recommendations in this report reflect the state of published evidence at the time of analysis. Xcode Life does not guarantee that the information in this report remains current indefinitely, and we recommend periodic review with your healthcare provider, particularly if you are starting a new medication for which PGx guidance has recently been updated.")}

        {_h2("Analytical Limitations")}
        {_para("Genotype calls are derived from the raw DNA file submitted by the user, analyzed using the PharmCAT v3.2.0 pipeline. Accuracy depends on the quality of the submitted data. Certain genes cannot be reliably typed using standard VCF-based analysis and are therefore excluded from clinical recommendations. Results for genes with complex structural variation may include multiple possible diplotype assignments where phasing is uncertain.")}

        {_h2("No Guarantee of Outcomes")}
        {_para("Genetic information provides probabilistic guidance, not certainty. Having a particular metabolizer phenotype does not guarantee that you will or will not experience a specific response or side effect. Individual drug responses are influenced by many non-genetic factors that this report cannot account for.")}

        {_h2("Privacy")}
        {_para("Your genetic data is sensitive personal health information. Xcode Life handles all data in accordance with applicable privacy laws and our privacy policy. This report is intended for your personal use and, where relevant, for sharing with your healthcare providers. Your raw DNA data file is processed solely to generate this report and is securely deleted from our active servers within 24 hours of report generation.")}

        <div style="border:1px solid #4DB7D0; border-left:4px solid #4DB7D0; border-radius:6px; background:#f0fffe; padding:10px 16px; margin-top:8px; page-break-inside: avoid; break-inside: avoid;">
            <div style="font-size:11.5px; font-weight:700; color:#023D79; margin-bottom:4px;">In Summary</div>
            <div style="font-size:11.5px; line-height:1.6; color:#374151;">
                This report is a tool to help you and your doctor have a more informed conversation about your medications. Use it alongside professional medical guidance, not instead of it.
            </div>
        </div>
    </div>
    """
    return _wrap_page(page_content, patient_name, pg, page_id="disclaimer")
# ============================================================================
# TABLE OF CONTENTS
# ============================================================================

def toc_template(tocitems, name: str, pg: int) -> str:
    cat_colors = ['#e3f2fd', '#fce4ec', '#e8f5e9', '#fff8e1', '#f3e5f5', '#e0f7fa']
    cat_index  = 0

    static_anchors = {
        'About Personalized Medicine':    ('welcome',           'About Personalized Medicine'),
        'About This Report':              ('how_to_read',       'About This Report'),
        'How to Understand Your Results': ('faqs',              'How to Understand Your Results'),
        'Note to Doctor':                 ('doctor',            'Note to Doctor'),
        'Table of Contents':              ('tocpage',           'Table of Contents'),
        'Summary of Medication Insights': ('executive_summary', 'Summary of Medication Insights'),
        'Coverage Statement':             ('coverage_statement', 'Coverage Statement'),
        'Other Evaluated Medications':    ('other_evaluated',   'Other Evaluated Medications'),
        'Genes Requiring Specialized Testing': ('specialized_genes', 'Genes Requiring Specialized Testing'),
        'Genotype Summary':               ('genotype_summary',  'Genotype Summary'),
        'Medications Evaluated — No Current Guidance': ('no_guideline_available', 'No Guideline Available'),
        'No Guideline Available':         ('no_guideline_available', 'No Guideline Available'),
        'Disclaimer':                     ('disclaimer',        'Disclaimer'),
    }

    rows_html = ""
    current_category = None

    for kind, label, pagenum in tocitems:
        if kind == 'chapter':
            rows_html += f"""
            <tr style="background:#023D79;">
                <td colspan="2" style="padding:8px 14px; font-size:11.5px; font-weight:700; color:white; letter-spacing:0.3px;">
                    {label}
                </td>
            </tr>
            """

        elif kind in ('category', 'numbered_category'):
            bg = cat_colors[cat_index % len(cat_colors)]
            cat_index += 1
            current_category = label
            rows_html += f"""
            <tr>
                <td colspan="2" style="padding:7px 14px 7px 20px; background:{bg}; border-left:3px solid #4DB7D0; font-size:11.5px; font-weight:700; color:#023D79; letter-spacing:0.5px; text-transform:uppercase;">
                    {label}
                </td>
            </tr>
            """

        elif kind == 'drug':
            # Support optional 4th element as custom anchor (for multi-category drugs)
            if isinstance(label, tuple):
                display_label, anchor = label[0].title(), label[1]
            else:
                anchor = sanitize_id(label)
                display_label = label.title()
            rows_html += f"""
            <tr style="border-bottom:1px solid #f5f7fa;">
                <td style="padding:5px 14px 5px 24px; font-size:12px; color:#1a1a1a;">
                    <a href="#{anchor}" style="color:#1a1a1a; text-decoration:none;">{display_label}</a>
                </td>
                <td style="padding:5px 8px; font-size:12px; color:#374151; font-weight:700; text-align:right; width:55px; white-space:nowrap;">
                    <a href="#{anchor}" class="toc-page-number" aria-label="Page for {display_label}"></a>
                </td>
            </tr>
            """

        elif kind == 'static' or kind == 'section1':
            anchor_data = static_anchors.get(label)
            if anchor_data:
                anchor, display = anchor_data
            else:
                anchor  = label.lower().replace(' ', '_')
                display = label

            rows_html += f"""
            <tr style="border-bottom:1px solid #f0f4f8;">
                <td style="padding:6px 14px 6px 20px; font-size:11.5px; font-weight:500; color:#1a1a1a;">
                    <a href="#{anchor}" style="color:#1a1a1a; text-decoration:none;">{display}</a>
                </td>
                <td style="padding:6px 8px; font-size:11.5px; color:#374151; font-weight:700; text-align:right; width:55px; white-space:nowrap;">
                    <a href="#{anchor}" class="toc-page-number" aria-label="Page for {display}"></a>
                </td>
            </tr>
            """

    inner = f"""
        <div style="padding-top:10px;">
            <div class="page-title">Table of Contents</div>
            <table style="width:100%; border-collapse:collapse; margin-top:14px; border: 1.5px solid #e2eaf4; table-layout:fixed;">
                <thead style="display: table-header-group;">
                    <tr style="background:#4DB7D0;">
                        <th style="padding:12px 14px; text-align:left; font-size:12px; color:white; letter-spacing:0.5px; font-weight:700;">Section / Drug</th>
                        <th style="padding:12px 8px; text-align:right; font-size:12px; color:white; letter-spacing:0.5px; font-weight:700; width:55px;">Page</th>
                    </tr>
                </thead>
                <tbody>
                    {rows_html}
                </tbody>
            </table>
        </div>
        """
    return _wrap_page(inner, name, pg, page_id="tocpage")

def drug_crosslisted_template(
    drug_name,
    category,
    category_rows,
    patient_name,
    curr_pg,
    primary_category,
    primary_anchor,
    primary_rows=None,
    as_fragment=False,
):
    """Render a compact secondary occurrence for a cross-listed drug.

    The full pharmacogenomic monograph is rendered once under the first
    therapeutic category.  Additional category blocks are intended to sit
    immediately below that monograph for the same drug, preserving category-
    specific explanatory content while referring to the recommendation above.
    """

    def _safe(val):
        text = str(val).strip()
        return text if text.lower() not in {"", "nan", "none", "n/a"} else ""

    def _unique_texts(df, column):
        if df is None or df.empty or column not in df.columns:
            return []
        out = []
        for value in df[column].tolist():
            text = _safe(value)
            if text and text not in out:
                out.append(text)
        return out

    def _text_block(values):
        if not values:
            return '<div style="color:#6b7280; font-style:italic;">No description available.</div>'
        return ''.join(
            f'<div style="margin-bottom:{"7px" if i < len(values)-1 else "0"};">{text}</div>'
            for i, text in enumerate(values)
        )

    drug_id = sanitize_id(drug_name) + "__" + sanitize_id(category)
    fs_headers = "14px"
    fs_subheaders = "12px"
    fs_body = "11.5px"

    class_box = f"""
        <div class="crosslisted-category-heading" style="margin: 12px 0 10px 0; padding: 8px 14px; border-radius: 24px;
                    border: 1px solid #c9dff0; background: linear-gradient(90deg, #eef5fc 0%, #e0effa 100%);
                    text-align: center; font-size: {fs_headers}; font-weight: 700; color: #000000;
                    text-transform: uppercase; letter-spacing: 0.5px;">
            {category}
        </div>
    """

    drug_title_box = f"""
        <div class="drug-header"
             style="background: linear-gradient(90deg, #0984b6 0%, #63c0d3 100%); padding: 7px 12px; margin-bottom: 8px;">
            <h3 class="drug-title" style="font-size:{fs_subheaders}; font-weight:700; margin:0; color:#ffffff;">
                <span style="text-transform:uppercase;">{drug_name[0]}</span>{drug_name[1:]}
            </h3>
        </div>
    """

    about_values = _unique_texts(category_rows, "About this medication")
    about_block = f"""
        <div class="info-box crosslisted-about-box" style="margin-bottom:10px; border:1px solid #d0e4f7; border-left:4px solid #4DB7D0;
                    border-radius:8px; background:#f4f9ff; padding:9px 14px;">
            <div style="font-size:{fs_headers}; font-weight:700; color:#1a1a1a; margin-bottom:4px;">About this medication</div>
            <div style="font-size:{fs_body}; line-height:1.5; color:#333;">{_text_block(about_values)}</div>
        </div>
    """

    # Every therapeutic-category occurrence should carry its own explanatory
    # patient context when the GSI provides it.  Do NOT suppress this block just
    # because the text happens to be identical to the primary category.
    witm_col = "How this gene/phenotype affects the drug and what it means for you"
    current_witm = _unique_texts(category_rows, witm_col)
    impact_block = ""
    if current_witm:
        impact_block = f"""
        <div class="info-box crosslisted-impact-box" style="margin-bottom:10px; border:1px solid #e5e7eb; border-left:4px solid #94a3b8;
                    border-radius:8px; background:#ffffff; padding:9px 14px;">
            <div style="font-size:{fs_headers}; font-weight:700; color:#1a1a1a; margin-bottom:4px;">What it means for you</div>
            <div style="font-size:{fs_body}; line-height:1.5; color:#333;">{_text_block(current_witm)}</div>
        </div>
        """

    # Compact patient finding from the already patient-matched Step 4 rows.
    finding_pairs = []
    if category_rows is not None and not category_rows.empty:
        for row in category_rows.to_dict("records"):
            gene = _safe(row.get("Gene", ""))
            pheno = _safe(row.get("Phenotype", ""))
            if not gene or not pheno:
                continue
            pheno = re.sub(r'\(as:([^)]+)\)', lambda m: f'(AS:{m.group(1)})', pheno, flags=re.IGNORECASE)
            pair = f"<strong>{gene}</strong> — {pheno}"
            if pair not in finding_pairs:
                finding_pairs.append(pair)
    finding_html = "; &nbsp; ".join(finding_pairs[:6])
    if len(finding_pairs) > 6:
        finding_html += "; …"

    primary_label = str(primary_category).title()
    cross_ref = f"""
        <div style="margin-bottom:10px;">
            <div style="font-size:{fs_headers}; font-weight:700; color:#1a1a1a; margin-bottom:5px;">Clinical Recommendations</div>
            <div class="crosslisted-reference-box" style="padding:10px 14px; border:1px solid #bfdbfe; border-left:4px solid #2563eb;
                        border-radius:6px; background:#eff6ff; font-size:{fs_body}; line-height:1.55; color:#1e3a5f;">
                See the full pharmacogenomic recommendation above under
                <strong>{primary_label} — {drug_name}</strong>.
                {f'<div style="margin-top:7px; padding-top:7px; border-top:1px solid #dbeafe; color:#374151;"><strong>Key genetic finding:</strong> {finding_html}</div>' if finding_html else ''}
            </div>
        </div>
    """

    content = f"""
    <div class="crosslisted-category-block" id="{drug_id}" style="padding-top:6px;">
        {class_box}
        {drug_title_box}
        {about_block}
        {impact_block}
        {cross_ref}
    </div>
    """
    if as_fragment:
        return content
    return _wrap_page(content, patient_name, curr_pg, page_id=drug_id)


def drug_detail_template(drug_name, category, drug_rows, patient_name, curr_pg, is_section_start=False, coverage_categories=None, master_genes_df=None, drug_gene_map=None, drug_gene_catalog=None, rs12777823_in_step1: bool = False, per_drug_overrides=None, extra_category_blocks=None):
    # Note: "uncertain susceptibility" is a valid phenotype with recommendations from PharmCAT/GSI
    # (used for malignant hyperthermia genes), so it is NOT in the _BAD set
    _BAD = {"indeterminate", "no call", "no data available", "unknown", "no result", "nan", ""}
    # Phenotypes that indicate the gene simply wasn't covered by the consumer
    # platform (chip or pipeline never attempted it) — distinguished from
    # "indeterminate" which means the platform attempted it but couldn't call.
    _PLATFORM_LIMIT_PHENOS = {"no result", "no data available", "nan", ""}
    _extra_category_html = "".join(extra_category_blocks or [])

    # RYR1 genotype check for clinical override
    ryr1_status = None
    if master_genes_df is not None and not master_genes_df.empty:
        ryr1_match = master_genes_df[master_genes_df['Gene'].str.upper() == 'RYR1']
        if not ryr1_match.empty:
            pheno_val = str(ryr1_match.iloc[0].get('Phenotype', '')).strip().lower()
            is_valid_call = pheno_val and not any(term in pheno_val for term in ["indeterminate", "no call", "no result", "unknown", "uncallable"])
            if is_valid_call:
                if any(x in pheno_val for x in ["susceptibility", "positive", "malignant"]):
                    ryr1_status = 'positive'
                else:
                    ryr1_status = 'negative'

    def _safe(val):
        s = str(val).strip()
        return s if s.lower() not in ('', 'nan', 'none', 'n/a') else ''
        
    def style_phenotype(pheno):
        # Capitalise activity-score suffix: (as:1.0) → (AS:1.0)
        pheno = re.sub(r'\(as:([^)]+)\)', lambda m: f'(AS:{m.group(1)})', str(pheno))

        # Format metabolizer phenotypes: "normal metabolizer" → "Normal Metabolizer"
        pheno = re.sub(r'\bnormal\s+metabolizer\b', 'Normal Metabolizer', pheno, flags=re.IGNORECASE)
        pheno = re.sub(r'\bpoor\s+metabolizer\b', 'Poor Metabolizer', pheno, flags=re.IGNORECASE)
        pheno = re.sub(r'\bintermediate\s+metabolizer\b', 'Intermediate Metabolizer', pheno, flags=re.IGNORECASE)
        pheno = re.sub(r'\bultrarapid\s+metabolizer\b', 'Ultrarapid Metabolizer', pheno, flags=re.IGNORECASE)
        pheno = re.sub(r'\brapid\s+metabolizer\b', 'Rapid Metabolizer', pheno, flags=re.IGNORECASE)
        pheno = re.sub(r'\bextensive\s+metabolizer\b', 'Extensive Metabolizer', pheno, flags=re.IGNORECASE)

        p_lower = pheno.lower()
        if 'normal' in p_lower or 'metabolizer' in p_lower:
            return f'<span style="font-weight: 600; color: #374151;">{pheno}</span>'
        elif 'no data' in p_lower or 'indeterminate' in p_lower or 'unknown' in p_lower or 'unassigned' in p_lower:
            return f'<span style="background: #f3f4f6; color: #6b7280; padding: 3px 8px; border-radius: 6px; font-weight: 600; display: inline-block;">{pheno}</span>'
        else:
            return f'<span style="font-weight: 600; color: #374151;">{pheno}</span>'

    drug_id       = sanitize_id(drug_name) + "__" + sanitize_id(category)
    fs_headers    = "14px"
    fs_subheaders = "12px"
    fs_body       = "11.5px"

    # ── Derive Clinical Action badge ──────────────────────────────────────────
    # ── No-action phrases ─────────────────────────────────────────────────────
    # Must match the executive summary list so drug badges and summary buckets
    # stay consistent.
    _NO_ACTION_PHRASES = [
        "no specific action", "no action required", "no action is needed",
        "no change needed", "no change is needed",
        "no dose adjustment", "no dose adjustment is anticipated",
        "no need to avoid", "no indication to", "no reason to avoid",
        "no alternative needed", "no recommendation", "no adjustment",
        "no clinical action", "no clinically significant",
        "prescribe desired starting dose",
        "start with normal starting dose",     # thiopurine NM phrase
        "normal starting dose",
        "standard recommended dose", "standard recommended starting",
        "standard starting daily dose", "standard starting dose",
        "recommended starting dose", "recommended starting dosing",
        "initiate therapy with recommended", "initiate therapy with standard",
        "initiate standard dose", "initiate standard",
        "initiate therapy at standard dose",
        "use standard dose", "prescribe standard",
        "standard dose recommended",
        # Routine titration / TDM — NOT PGx-specific interventions
        "adjust doses based on disease",
        "adjust doses based on clinical",
        "adjust doses based on specific population",
        "based on disease-specific guidelines",
        "based on disease-specific and specific population",
        "use therapeutic drug monitoring to guide",
        "therapeutic drug monitoring to guide dose",
        "subsequent doses should be adjusted according to therapeutic",
        "standard prescribing", "standard dosing",
        "standard of care dosing",
        "routine dosing algorithms",
    ]

    # ── Action-required keywords — specific multi-word phrases ─────────────────
    # "adjust dose" and bare "alternative" removed — they match routine clinical
    # guidance ("adjust doses based on disease-specific guidelines") and create
    # false positives for Normal Metabolizer rows.
    _AR_KEYWORDS = [
        "contraindicated",
        "do not use", "do not prescribe",
        "consider an alternative",
        "use an alternative", "use alternative",
        "select an alternative",
        "alternative agent", "alternative drug",
        "dose reduction",
        "reduce dose", "reduce the dose",
        "decrease dose", "decreased dose",
        "decreased starting dose",
        "lower starting dose", "reduce starting dose",
        "not recommended",
        "increase the dose",      # ultrarapid metabolizer guidance
        "inadequate drug efficacy",
    ]

    # ── Use-with-caution keywords — specific phrases ───────────────────────────
    # Bare "monitor"/"consider"/"may increase" removed; they match routine TDM
    # instructions and dose-titration notes that are NOT PGx-specific cautions.
    _UWC_KEYWORDS = [
        "closely monitor", "close monitoring",
        "more frequent monitoring", "additional monitoring",
        "enhanced monitoring", "increased monitoring",
        "monitor for adverse", "monitor for toxicity",
        "use with caution", "prescribe with caution", "exercise caution",
        "proceed with caution",
        "increased risk of adverse", "increased risk of toxicity",
        "increased risk of side effect",
        "at increased risk of",
        "may increase the risk", "may increase toxicity",
        "may increase plasma", "may increase serum",
        "may decrease efficacy", "may decrease plasma",
        "lower dose", "start at a lower dose",
        "slower titration", "slower dose titration",
    ]

    # Strips "Other Considerations" / "Implications" section from rec text so
    # the no-action / AR checks only see the patient-specific guidance.
    _OC_SPLIT_RE = re.compile(
        r'\s*\b(?:Other Considerations|Implications|Affected Subgroup|'
        r'Submitted Genotype|Recommendation:)\b',
        re.IGNORECASE,
    )

    def _trim_rec(text: str) -> str:
        head = _OC_SPLIT_RE.split(text, maxsplit=1)[0]
        return re.sub(r'^\s*Recommendation\s*:?\s+', '', head, flags=re.IGNORECASE).strip()

    def _no_action_label(pheno: str, clf: str) -> str:
        """Return a phenotype-aware label when no PGx action is needed."""
        pheno = pheno.lower()
        clf   = clf.lower()
        if "normal metabolizer" in pheno:
            return "Normal metabolism — standard dosing"
        if "normal function" in pheno or pheno in ("normal",):
            return "Normal gene function — standard dosing"
        if "intermediate metabolizer" in pheno:
            return "Intermediate metabolism — standard dosing"
        if "ultrarapid metabolizer" in pheno or "ultrarapid" in pheno:
            return "Ultrarapid metabolism — follow clinical guidelines"
        if "poor function" in pheno:
            return "Reduced gene function — standard dosing"
        if "indeterminate" in pheno or "no result" in pheno or "no data" in pheno:
            return "Genotype undetermined — use clinical guidelines"
        if "no recommendation" in clf or clf == "unspecified":
            return "No PGx-specific guidance for your genotype"
        return "Standard prescribing applies"

    def _derive_action(rows):
        """
        Build action badge from drug_rows.

        Candidates: (clf_rank, act_priority, rec_text, clf_text, pheno_text)
          act_priority: 0=AR, 1=UWC, 2=no-action
        Sort ascending → most-actionable + strongest-evidence row wins.
        """
        clf_priority = {"strong": 0, "moderate": 1, "optional": 2}
        candidates = []
        for r in rows.to_dict("records"):
            clf   = str(r.get("Classification", "")).strip().lower()
            rec   = str(r.get("Recommendation", "")).strip()
            pheno = str(r.get("Phenotype", "")).strip()
            rank  = clf_priority.get(clf, 99)
            patient_rec = _trim_rec(rec).lower()
            is_no_action = any(p in patient_rec for p in _NO_ACTION_PHRASES)
            is_ar = not is_no_action and (
                any(k in patient_rec for k in _AR_KEYWORDS) or
                # "avoid" handled separately so "no need to avoid" can't slip through
                ("avoid" in patient_rec and "no need to avoid" not in patient_rec
                 and "need not avoid" not in patient_rec) or
                clf == "strong"
            )
            is_uwc = not is_no_action and not is_ar and (
                any(k in patient_rec for k in _UWC_KEYWORDS) or clf == "moderate"
            )
            act_priority = 0 if is_ar else (1 if is_uwc else 2)
            candidates.append((rank, act_priority, rec, clf, pheno))

        # Sort: actionability first, then classification strength.
        candidates.sort(key=lambda x: (x[1], x[0]))
        if not candidates:
            return "Standard Use", "Standard prescribing applies", "#15803d", "#f0fdf4"

        best_rec, best_clf, best_pheno = candidates[0][2], candidates[0][3], candidates[0][4]
        rec = _trim_rec(best_rec).lower()
        clf = best_clf.lower()

        # ── No action ──────────────────────────────────────────────────────────
        if any(p in rec for p in _NO_ACTION_PHRASES):
            return "Standard Use", _no_action_label(best_pheno, best_clf), "#15803d", "#f0fdf4"

        # ── Action required ───────────────────────────────────────────────────
        if (any(k in rec for k in _AR_KEYWORDS) or
                ("avoid" in rec and "no need to avoid" not in rec) or
                clf == "strong"):
            if "contraindicated" in rec:
                lbl = "Contraindicated"
            elif "do not use" in rec or ("avoid" in rec and "no need to avoid" not in rec):
                lbl = "Avoid — use alternative"
            elif "consider an alternative" in rec or "use an alternative" in rec:
                lbl = "Consider alternative drug"
            elif "not recommended" in rec:
                lbl = "Not recommended"
            elif "dose reduction" in rec or "reduce dose" in rec or "reduce the dose" in rec:
                lbl = "Dose reduction recommended"
            elif "alternative" in rec:
                lbl = "Alternative drug recommended"
            else:
                lbl = "Action required"
            return "Action Required", lbl, "#dc2626", "#fef2f2"

        # ── Use with caution ──────────────────────────────────────────────────
        if any(k in rec for k in _UWC_KEYWORDS) or clf == "moderate":
            if "closely monitor" in rec or "close monitoring" in rec:
                lbl = "Monitor closely"
            elif "increased monitoring" in rec or "additional monitoring" in rec:
                lbl = "Increased monitoring recommended"
            elif "caution" in rec:
                lbl = "Use with caution"
            elif "lower dose" in rec or "start at a lower dose" in rec:
                lbl = "Start with lower dose"
            elif "increased risk" in rec or "greater risk" in rec:
                lbl = "Increased risk — exercise caution"
            elif "slower titration" in rec:
                lbl = "Slower dose titration advised"
            else:
                lbl = "Use with caution"
            return "Use With Caution", lbl, "#d97706", "#fffbeb"

        return "Standard Use", _no_action_label(best_pheno, best_clf), "#15803d", "#f0fdf4"

    _action_bucket, _action_label_text, _action_color, _action_bg = _derive_action(drug_rows)

    header_html = ""
    if is_section_start:

        header_html = (
            '<div class="page-title" id="detailed_drug_reports">Drug-Specific Recommendations</div>'
            '<div style="margin-bottom:14px; padding:10px 16px; background:#f0fffe; '
            'border-left:4px solid #4DB7D0; border-radius:6px; font-size:12px; '
            'color:#374151; line-height:1.6;">'
            'This section provides detailed pharmacogenomic recommendations for the medications identified in your profile. '
            'All medications with pharmacogenomic relevance are shown. '
            '<em>Actions are indicated where applicable.</em>'
            '</div>'
        )

    class_box = f"""
        <div style="margin: 12px 0 10px 0; padding: 8px 14px; border-radius: 24px; border: 1px solid #c9dff0; background: linear-gradient(90deg, #eef5fc 0%, #e0effa 100%); text-align: center; font-size: {fs_headers}; font-weight: 700; color: #000000; text-transform: uppercase; letter-spacing: 0.5px;">
            {category}
        </div>
    """

    # --- Prodrug badge ---
    _is_prodrug = str(drug_rows.iloc[0].get("IsProdrug", "")).strip().lower() == "yes" if not drug_rows.empty else False
    prodrug_badge = ""
    if _is_prodrug:
        prodrug_badge = (
            '<span style="display:inline-block; margin-left:10px; padding:2px 8px; '
            'font-size:12px; font-weight:700; color:#0984b6; background:#e0f4fa; '
            'border:1px solid #63c0d3; border-radius:6px; vertical-align:middle; '
            'text-transform:uppercase; letter-spacing:0.5px;">Prodrug</span>'
        )

    drug_title_box = f"""
        <div class="drug-header" style="background: linear-gradient(90deg, #0984b6 0%, #63c0d3 100%); padding: 7px 12px; margin-bottom: 6px;">
            <h3 class="drug-title" style="font-size: {fs_subheaders}; font-weight: 700; margin: 0; color: #ffffff;"><span style="text-transform:uppercase;">{drug_name[0]}</span>{drug_name[1:]}{prodrug_badge}</h3>
        </div>
    """

    # --- 1. BUILD GENES TABLE (4 COLUMNS: Gene | Diplotype | Phenotype | Status) ---
    drug_key = str(drug_name).strip().lower()

    def _normalise_gene(g: str) -> str:
        """Lowercase rsid-like gene names so RS12777823 == rs12777823."""
        s = g.strip()
        if s.lower().startswith("rs") and s[2:].isdigit():
            return s.lower()
        return s

    genes_from_map = set()
    if drug_gene_map and drug_key in drug_gene_map:
        genes_from_map = {_normalise_gene(g) for g in drug_gene_map[drug_key]}

    # ── Separation: PharmCAT-analysed vs GSI-catalog-only ────────────────────
    genes_from_rows = set()
    for row in drug_rows.to_dict("records"):
        raw_g = _safe(row.get('Gene', ''))
        for sep in [';', '\n', ',']:
            if sep in raw_g:
                genes_from_rows.update([_normalise_gene(x) for x in raw_g.split(sep) if x.strip()])
                break
        else:
            if raw_g:
                genes_from_rows.add(_normalise_gene(raw_g))

    genes_from_catalog = set()
    if drug_gene_catalog and drug_key in drug_gene_catalog:
        genes_from_catalog = {_normalise_gene(g) for g in drug_gene_catalog[drug_key]}
    if drug_gene_map and drug_key in drug_gene_map:
        genes_from_catalog.update({_normalise_gene(g) for g in drug_gene_map[drug_key]})

    # ONLY show genes in the top table if Step 4 actively merged them for this drug
    genes_analyzed = genes_from_rows.copy()

    # ── Yellow-box genes (genes_catalog_only) ────────────────────────────────
    # We want to flag genes that are known to affect this drug (per GSI catalog
    # or step4 rows) but that PharmCAT did NOT specifically recommend for this
    # drug in its output (i.e. not in drug_gene_map for this drug).
    #
    # Previous (broken) logic:  genes_from_catalog − genes_from_rows
    #   → always empty because every gene in the catalog is also in drug_rows
    #     (step4 already merged them), so the subtraction gives nothing.
    #
    # Fixed logic:  (genes_from_rows ∪ genes_from_catalog) − pharmcat_genes
    #   pharmcat_genes = genes PharmCAT actually output recommendations for
    #                    (from step3 "2_ALL_RECOMMENDATIONS" sheet).
    #   The remainder = genes that came from GSI / catalog only.
    _pharmcat_genes_for_drug = (
        {_normalise_gene(g) for g in drug_gene_map.get(drug_key, set())}
        if drug_gene_map else set()
    )
    genes_catalog_only = (genes_from_rows | genes_from_catalog) - _pharmcat_genes_for_drug

    # Clinical Override: Suppress CACNA1S if RYR1 is positive (MHS)
    if ryr1_status == 'positive':
        genes_analyzed = {g for g in genes_analyzed if _normalise_gene(g).lower() != 'cacna1s'}
        genes_catalog_only = {g for g in genes_catalog_only if _normalise_gene(g).lower() != 'cacna1s'}

    # CRITICAL: Remove GSI-only genes from the main analyzed table
    # They should ONLY appear in the yellow box, NOT in the main "Genes Analyzed" table
    genes_analyzed = genes_analyzed - genes_catalog_only

    # --- CUSTOM WARFARIN OVERRIDE FOR rs12777823 ---
    if drug_key == "warfarin":
        # rs12777823_in_step1 is passed from build_report which reads the actual
        # step1 VCF file for this patient (results/step1_<sample_id>.vcf).
        # Route accordingly:

        if rs12777823_in_step1:
            # Found in step1 VCF → include in yellow info box alongside other GSI genes
            genes_catalog_only.add("rs12777823")
        else:
            # Not found in step1 VCF → keep hidden entirely
            genes_catalog_only.discard("rs12777823")

    # Genes shown in the "Genes Analyzed" table: ONLY genes that were analyzed by PharmCAT
    # GSI-only genes (in genes_catalog_only) do NOT appear in the main table, only in yellow box
    genes_to_evaluate = sorted(list(genes_analyzed))

    # Build per-gene status and phenotype maps
    gene_status_map    = {}   # gene -> "Analyzed" | "Requires Specialized Testing" | "Not Evaluated by Platform" | "Not Tested"
    gene_pheno_map     = {}   # gene -> list of (diplotype_str, phenotype_str)
    gene_collapsed_map = {}   # gene -> True if multiple diplotypes were collapsed into one row

    for g in genes_to_evaluate:
        if not g:
            continue
        if master_genes_df is not None and not master_genes_df.empty:
            match = master_genes_df[master_genes_df['Gene'].str.lower() == g.lower()]
            if match.empty:
                gene_status_map[g] = "Not Tested"
                gene_pheno_map[g]  = [("—", "—")]
            else:
                pheno_dict = {}
                for m_row in match.to_dict("records"):
                    raw_d = _safe(m_row.get('Diplotype', '')) or '—'
                    raw_p = _safe(m_row.get('Phenotype', '')) or 'No data available'
                    if ":" in raw_p and raw_p.lower().startswith(g.lower()):
                        raw_p = raw_p.split(":", 1)[1].strip()
                    if raw_p.lower() in ["no result", "not called", "unknown/unknown", "unknown", "unassigned", "n/a", "", "nan", "uncategorized"]:
                        raw_p = "No data available"
                    if raw_p not in pheno_dict:
                        pheno_dict[raw_p] = []
                    if raw_d and raw_d not in pheno_dict[raw_p] and raw_d != "—":
                        pheno_dict[raw_p].append(raw_d)

                if len(pheno_dict) > 1 and "No data available" in pheno_dict:
                    del pheno_dict["No data available"]

                # Collapse multi-row genes into ONE display row when multiple
                # phenotypes are present.  CYP4F2 unphased gives 6 distinct
                # diplotypes, NAT2 can give 70+, SLCO1B1 typically 5.  Showing
                # them as separate rows produces noise; one row with comma-
                # joined diplotypes is much more readable while preserving the
                # full information.
                if len(pheno_dict) > 1:
                    seen_d, all_dips = set(), []
                    all_phens = []
                    for p_str, d_vals in pheno_dict.items():
                        if p_str not in all_phens:
                            all_phens.append(p_str)
                        for d in d_vals:
                            if d and d not in seen_d:
                                seen_d.add(d)
                                all_dips.append(d)
                    collapsed_d = ", ".join(all_dips) if all_dips else "—"
                    collapsed_p = " / ".join(all_phens) if all_phens else "—"
                    pheno_list = [(collapsed_d, collapsed_p)]
                    gene_collapsed_map[g] = True
                else:
                    pheno_list = [
                        (", ".join(d_vals) if d_vals else "—", p_str)
                        for p_str, d_vals in pheno_dict.items()
                    ]
                    # Multiple diplotypes within a single phenotype is also
                    # a phasing-uncertainty signal (e.g. SLCO1B1 *37/*37, *37/*42
                    # both producing "Indeterminate").  Flag it.
                    gene_collapsed_map[g] = any(
                        len(d_vals) > 1 for d_vals in pheno_dict.values()
                    )

                # Apply VKORC1 / IFNL3 / rs12777823 pretty-label translation.
                # Strategy:
                #   1. Try diplotype-key lookup (VKORC1/IFNL3): if the diplotype
                #      changes, we got a match — use both pretty values.
                #   2. Try phenotype-key lookup (rs12777823): the diplotype is
                #      already readable ("G/G (reference)") but the phenotype is
                #      an ugly HGVS string — translate only the phenotype.
                def _translate_row(d_str, p_str):
                    pretty_d, pretty_p = pretty_genotype_pair(g, d_str)
                    if pretty_d != d_str:
                        # Diplotype-key lookup succeeded (VKORC1 / IFNL3 style)
                        return pretty_d, pretty_p
                    # Phenotype-key lookup: try the raw phenotype as the key
                    _pd2, _pp2 = pretty_genotype_pair(g, p_str)
                    if _pp2 != p_str:
                        # Keep existing diplotype; replace ugly HGVS phenotype
                        return d_str, _pp2
                    return d_str, p_str

                pheno_list = [_translate_row(d, p) for (d, p) in pheno_list]
                gene_pheno_map[g] = pheno_list

                all_phenos_lower = [pt[1].lower() for pt in pheno_list]
                if any(p in _BAD for p in all_phenos_lower):
                    # Distinguish: all bad phenos are "no result/no data" →
                    # platform simply didn't cover this gene (e.g. CYP3A4 on
                    # consumer arrays). Otherwise, the platform attempted
                    # calling but got an indeterminate result → specialized
                    # testing needed.
                    bad_phenos_hit = [p for p in all_phenos_lower if p in _BAD]
                    if all(p in _PLATFORM_LIMIT_PHENOS for p in bad_phenos_hit):
                        gene_status_map[g] = "Not Evaluated by Platform"
                    else:
                        gene_status_map[g] = "Requires Specialized Testing"
                else:
                    gene_status_map[g] = "Analyzed"
        else:
            gene_status_map[g] = "Not Tested"
            gene_pheno_map[g]  = [("—", "—")]

    # Sort genes so "Not Evaluated by Platform" appears at the bottom
    def _gene_sort_key(g_name):
        st = gene_status_map.get(g_name, "Not Tested")
        if st == "Not Evaluated by Platform":
            return (1, g_name)
        return (0, g_name)
    genes_to_evaluate.sort(key=_gene_sort_key)

    # Build HTML rows for the 4-column table
    # Special-case: hydralazine ordering — prefer Intermediate then Indeterminate
    try:
        if str(drug_name).strip().lower() == 'hydralazine':
            for gg, plist in list(gene_pheno_map.items()):
                # sort phenotypes so 'Intermediate' comes before 'Indeterminate'
                def _hyd_sort_key(item):
                    ph = (item[1] or '').lower()
                    if 'intermediate' in ph:
                        return 0
                    if 'indeterminate' in ph:
                        return 2
                    return 1
                gene_pheno_map[gg] = sorted(plist, key=_hyd_sort_key)
    except Exception:
        pass
    gene_rows_html = ""
    for g in genes_to_evaluate:
        if not g:
            continue
        status     = gene_status_map.get(g, "Not Tested")
        pheno_list = gene_pheno_map.get(g, [("—", "—")])

        if status == "Analyzed":
            status_badge = '<span style="color:#374151; font-size:11.5px; font-weight:700;">&#10003; Analyzed</span>'
        elif status == "Requires Specialized Testing":
            status_badge = '<span style="color:#374151; font-size:11.5px; font-weight:700;">&#9651; Specialized testing needed</span>'
        elif status == "Not Evaluated by Platform":
            status_badge = '<span style="background:#f0f4ff; color:#4b5da8; padding:2px 7px; border-radius:5px; font-size:11.5px; font-weight:700;">&#8857; Not evaluated by platform</span>'
        else:
            status_badge = '<span style="background:#f3f4f6; color:#6b7280; padding:2px 7px; border-radius:5px; font-size:11.5px; font-weight:700;">&#8212; Not Tested</span>'

        # For genes that couldn't be analyzed (Spec. Testing / Not Tested),
        # render exactly ONE row regardless of how many pheno_list entries
        # exist — extra rows would all be "— / —" noise with confusing
        # "(unphased)" labels that don't apply.
        rows_to_render = pheno_list if status == "Analyzed" else [("—", "—")]  # non-analyzed statuses always show single blank row

        is_collapsed = gene_collapsed_map.get(g, False) and status == "Analyzed"
        unphased_tag = (
            ' <span style="font-size:12px; color:#c2410c; font-style:italic; '
            'font-weight:500;">(unphased)</span>'
            if is_collapsed else ''
        )

        for (d_str, p_str) in rows_to_render:
            styled_p = style_phenotype(p_str) if status == "Analyzed" else f'<span style="color:#9ca3af; font-style:italic;">{p_str}</span>'
            gene_rows_html += (
                f'<tr>'
                f'<td style="font-weight:600; font-size:{fs_body}; color:#111827;">{g}{unphased_tag}</td>'
                f'<td style="font-size:{fs_body}; color:#374151; word-break:break-word;">{d_str}</td>'
                f'<td style="font-size:{fs_body}; color:#374151;">{styled_p}</td>'
                f'<td style="font-size:{fs_body};">{status_badge}</td>'
                f'</tr>'
            )
            # rows_to_render is now always 1 item (single row per gene after
            # the collapse step above), so we never produce additional rows here.

    if not gene_rows_html:
        gene_rows_html = f'<tr><td colspan="4" style="text-align:center; color:#9ca3af; font-size:{fs_body}; font-style:italic;">No data available</td></tr>'

    # ── Patient phenotype index (gene → set of normalised phenotypes) ─────────
    # Used below to filter Clinical Recommendation rows to only those matching
    # this patient's actual phenotype for each gene.  Phenotypes are normalised
    # to lowercase and stripped of "(AS:x.x)" suffixes so PharmCAT's plain
    # "Normal Metabolizer" matches GSI's "Normal Metabolizer (AS:2.0)".
    # NB: gene_pheno_map now contains pretty-translated values for VKORC1/IFNL3
    # ("-1639 G/G") which would NOT match the raw GSI rsID phenotype.  Rebuild
    # the index from the original master_genes_df so the filter compares the
    # raw PharmCAT phenotype against the raw GSI phenotype.  Also include the
    # pretty form as a fallback so manually-edited data with pretty labels
    # still matches.
    _patient_phen_per_gene = {}
    if master_genes_df is not None and not master_genes_df.empty:
        for _mr in master_genes_df.to_dict("records"):
            _g_raw = str(_mr.get("Gene", "")).strip()
            _p_raw = str(_mr.get("Phenotype", "")).strip()
            if not _g_raw:
                continue
            if ":" in _p_raw and _p_raw.lower().startswith(_g_raw.lower()):
                _p_raw = _p_raw.split(":", 1)[1].strip()
            _p_norm = re.sub(r"\s*\(as:[0-9.]+\)", "", _p_raw.lower().strip())
            _patient_phen_per_gene.setdefault(_g_raw.lower(), set()).add(_p_norm)
    # Append pretty-form values so a stripped-down test setup still matches
    for _g, _plist in gene_pheno_map.items():
        for (_d_str, _p_str) in _plist:
            for _frag in str(_p_str).split(" / "):
                _p_norm = re.sub(r"\s*\(as:[0-9.]+\)", "", _frag.lower().strip())
                if _p_norm:
                    _patient_phen_per_gene.setdefault(_g.lower(), set()).add(_p_norm)

    # Genes the patient cannot have personalised recommendations for
    # (Indeterminate / No Call / Spec Testing).  Used to skip those genes'
    # CPIC/DPWG/FDA rows entirely in the recommendations section.
    _spec_testing_genes_lower = {
        _g.lower() for _g, _st in gene_status_map.items()
        if _st in ("Requires Specialized Testing", "Not Tested", "Not Evaluated by Platform")
    }

    # Detect placeholder GSI text like "CPIC has no recommendations for X and Y."
    # or "CPIC provides no recommendations for warfarin and CYP2C9 Unspecified…".
    # These rows add noise without value — they aren't real recommendations.
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
        # Do NOT auto-suppress "not evaluated" rows — those carry the legitimate
        # "no guidance available" note from GSI (e.g. atazanavir+CYP2C19) and
        # should reach the renderer as muted info cards.
        if t in ("no recommendation", "no action required",
                 "not annotated", "n/a"):
            return True
        return False

    # ── Recommendation text formatter ────────────────────────────────────────
    # GSI Details often arrive as one run-on paragraph that mashes the CPIC
    # structure together, e.g.
    #   "There are recommendations for 3 different groups Patients with
    #    cardiovascular indications ... Recommendation Avoid standard dose ...
    #    Other Considerations For cardiovascular ... Implications CYP2C19: ...
    #    ──────── Patients with neurovascular ... Recommendation ..."
    # The formatter splits this into:
    #   • Subgroup blocks (when divider lines or "There are recommendations
    #     for N different groups" preamble are present)
    #   • Bold-prefixed sections for Recommendation / Implications /
    #     Other Considerations / Affected Subgroup / Submitted Genotype
    # Returns (html, is_structured).  When is_structured=True the caller
    # should NOT wrap the result in another "Recommendation:" prefix.
    # Case-SENSITIVE keyword match — GSI uses Title Case for section markers
    # ("Recommendation", "Implications", "Other Considerations").  Lower-case
    # uses inside content text (e.g. "No recommendation") must NOT trigger a
    # split.  The intro phrase regex stays case-insensitive (rare wording).
    _REC_KW_RE = re.compile(
        r'\b(Recommendation|Implications?|Other\s+Considerations|'
        r'Affected\s+Subgroup|Submitted\s+Genotype)\b\s*:?\s+'
    )
    _REC_DIVIDER_RE = re.compile(r'\s*[─━—―]{6,}\s*')
    _REC_INTRO_RE = re.compile(
        r'^\s*There\s+are\s+recommendations?\s+for\s+\d+\s+different\s+groups?\s*',
        re.IGNORECASE,
    )

    def _format_rec_text(text: str):
        if not text:
            return "", False
        text = str(text).strip()
        # Remove the "full FDA table" cruft that follows FDA Table quotes.
        text = re.sub(r'\s*"?\s*full\s+FDA\s+table\s*', ' ', text, flags=re.IGNORECASE)
        text = re.sub(r'\s{2,}', ' ', text).strip()

        if not _REC_KW_RE.search(text):
            return text, False

        is_multi = bool(_REC_INTRO_RE.match(text))
        body = _REC_INTRO_RE.sub('', text) if is_multi else text

        sections = _REC_DIVIDER_RE.split(body)
        rendered_sections = []
        for sec in sections:
            sec = sec.strip()
            if not sec:
                continue
            parts = _REC_KW_RE.split(sec)
            html_chunks = []
            preamble = parts[0].strip() if parts else ""
            if preamble:
                if is_multi:
                    # Multi-subgroup: preamble is a patient-group header
                    html_chunks.append(
                        f'<div style="font-weight:600; color:#1e3a5f; '
                        f'margin:0 0 6px 0; font-size:11.5px;">'
                        f'{preamble}</div>'
                    )
                else:
                    # Single recommendation: preamble IS the recommendation body
                    html_chunks.append(
                        f'<div style="margin-bottom:6px;">'
                        f'<strong style="color:#0f172a;">Recommendation:</strong> '
                        f'{preamble}</div>'
                    )
            for i in range(1, len(parts), 2):
                if i + 1 >= len(parts):
                    break
                kw = re.sub(r'\s+', ' ', parts[i]).strip().title()
                if kw.startswith("Implication"):
                    kw = "Implications"
                content = parts[i + 1].strip()
                if not content:
                    continue
                html_chunks.append(
                    f'<div style="margin-bottom:6px;">'
                    f'<strong style="color:#0f172a;">{kw}:</strong> {content}</div>'
                )
            if html_chunks:
                rendered_sections.append(''.join(html_chunks))

        if not rendered_sections:
            return text, False
        if len(rendered_sections) == 1:
            return rendered_sections[0], True
        joiner = ('<div style="height:1px; background:#cbd5e1; '
                  'margin:10px 0;"></div>')
        return joiner.join(rendered_sections), True

    # ── Optional: "Other Associated Genes" yellow info box ────────────────────
    # Lists genes known to affect this drug per GSI/CPIC/DPWG that are NOT in
    # PharmCAT's drug-specific recommendations. These genes are excluded from both
    # the Genes Analyzed table AND the Clinical Recommendations section to avoid
    # clutter and confusion — they're mentioned here instead.
    # Includes rs12777823 for warfarin when found in the patient's step1 VCF.

    _catalog_extra_html = ""

    # Yellow box should only show genes that are NOT already in the main Genes Analyzed table
    # (genes_catalog_only but not in genes_to_evaluate)
    _genes_for_yellow_box = genes_catalog_only - set(genes_to_evaluate)

    if _genes_for_yellow_box:
        _sorted_genes = sorted(_genes_for_yellow_box)

        # Format genes with proper English list formatting
        if len(_sorted_genes) == 1:
            _genes_str = f'<strong style="color:#92400e;">{_sorted_genes[0]}</strong>'

        elif len(_sorted_genes) == 2:
            _genes_str = (
                f'<strong style="color:#92400e;">{_sorted_genes[0]}</strong> and '
                f'<strong style="color:#92400e;">{_sorted_genes[1]}</strong>'
            )

        else:
            # 3+ genes: comma-separated with "and" before last
            _gene_parts = [
                f'<strong style="color:#92400e;">{g}</strong>'
                for g in _sorted_genes[:-1]
            ]

            _genes_str = (
                ', '.join(_gene_parts)
                + f' and <strong style="color:#92400e;">{_sorted_genes[-1]}</strong>'
            )

        # Single line statement combining all genes
        _is_plural = len(_sorted_genes) > 1
        _verb = 'are' if _is_plural else 'is'

        _genes_html = (
            f'<div style="margin-bottom:0px; margin-left:0px; '
            f'font-size:11.5px; line-height:1.6;">'
            f'{_genes_str} {_verb} known to affect this drug\'s metabolism '
            f'(GSI catalog), but was not analyzed by PharmCAT in this report.'
            f'</div>'
        )

        _catalog_extra_html = (
            f'<div style="margin-top:0; margin-bottom:10px; padding:6px 12px; '
            f'background:#fffbeb; border:1.5px solid #f59e0b; '
            f'border-left:4px solid #d97706; border-radius:6px; '
            f'font-size:11.5px; color:#78350f; line-height:1.5;">'
            f'<strong style="color:#92400e;">Other Associated Genes:</strong> '
            f'{_genes_html}'
            f'</div>'
        )


    # ── FDA DDI compounding alert ─────────────────────────────────────────────
    # IMPORTANT: this is OUTSIDE the "if _genes_for_yellow_box" block.
    # Therefore it can appear whether or not Other Associated Genes exists.

    _compounding_alert_html = generate_compounding_alert_html(
        drug_name=drug_name,
        gene_pheno_map=gene_pheno_map,
        gene_status_map=gene_status_map,
        fs_header=fs_headers,
        fs_body=fs_body,
    )


    genes_block = (
        f'<div style="font-size:{fs_headers}; font-weight:700; '
        f'color:#1a1a1a; margin:4px 0 3px 0;">Genes Analyzed</div>'

        f'<table style="margin-bottom:10px; border:1px solid #dde2ea; width:100%;">'

        f'<thead><tr>'
        f'<th style="width:15%; font-size:12px;">Gene</th>'
        f'<th style="width:30%; font-size:12px;">Diplotype</th>'
        f'<th style="width:35%; font-size:12px;">Phenotype</th>'
        f'<th style="width:20%; font-size:12px;">Status</th>'
        f'</tr></thead>'

        f'<tbody>{gene_rows_html}</tbody>'

        f'</table>'

        # First: Other Associated Genes, if applicable
        f'{_catalog_extra_html}'

        # Second: DDI block, if applicable
        f'{_compounding_alert_html}'
    )

    # --- 2. ABOUT THIS MEDICATION ---
    about = ""
    for r in drug_rows.to_dict("records"):
        if not about:
            about = _safe(r.get('About this medication', ''))
    about = about or 'No description available.'

    about_block = (
        f'<div class="info-box" style="margin-bottom:10px; border: 1px solid #d0e4f7; border-left: 4px solid #4DB7D0; border-radius: 8px; background: #f4f9ff; padding: 8px 14px;">'
        f'<div style="font-size:{fs_headers}; font-weight:700; color:#1a1a1a; margin-bottom:3px;">About this medication</div>'
        f'<div style="font-size:{fs_body}; line-height:1.5; color:#333;">{about}</div>'
        f'</div>'
    )

    # --- 3. WHAT IT MEANS FOR YOU (Clean Single Box with Flat Texts) ---
    analyzed_genes = [g for g in genes_to_evaluate if g and gene_status_map.get(g) == "Analyzed"]
    multi_gene = len(analyzed_genes) > 1

    impact_sections_html = ""
    has_unphased = False
    
    # Dictionary to group identical texts to prevent repeating paragraphs
    # Format: { "html content": ["CYP2C19", "CYP2D6"] }
    impact_groups = {}
    ordered_impacts = []

    for g in genes_to_evaluate:
        if not g: continue
        status = gene_status_map.get(g, "Not Tested")
        pheno_list = gene_pheno_map.get(g, [("—", "—")])

        if status == "Requires Specialized Testing":
            text = "The impact of this gene cannot be assessed because it requires specialized laboratory testing beyond standard DNA sequencing. Please refer to the Genes Requiring Specialized Testing section for details."
            block_html = f"<div style='margin-bottom:6px;'><span style='color:#333;'>{text}</span></div>"

            if block_html not in impact_groups:
                impact_groups[block_html] = []
                ordered_impacts.append(block_html)
            impact_groups[block_html].append(g)

        elif status == "Not Evaluated by Platform":
            text = ("This gene is not evaluated by consumer DNA genotyping platforms and was therefore not assessed in this analysis. "
                    "If genotyping of this gene is clinically relevant, a dedicated pharmacogenomic laboratory test may be requested through your healthcare provider.")
            block_html = f"<div style='margin-bottom:6px;'><span style='color:#333;'>{text}</span></div>"

            if block_html not in impact_groups:
                impact_groups[block_html] = []
                ordered_impacts.append(block_html)
            impact_groups[block_html].append(g)

        elif status == "Not Tested":
            text = "This gene was not evaluated in this analysis. Its potential impact has not been assessed."
            block_html = f"<div style='margin-bottom:6px;'><span style='color:#333;'>{text}</span></div>"

            if block_html not in impact_groups:
                impact_groups[block_html] = []
                ordered_impacts.append(block_html)
            impact_groups[block_html].append(g)

        else:
            is_unphased = gene_collapsed_map.get(g, False)
            if is_unphased: has_unphased = True

            # Build a map from identical explanatory text -> list of phenotypes/diplotypes
            # This consolidates repeated paragraphs (common for CYP4F2 and similar genes)
            how_map = {}
            for (d_str, p_str) in pheno_list:
                raw_phenos = [p.strip() for p in p_str.split(" / ")]
                individual_phenos = []
                for p in raw_phenos:
                    if p and p not in individual_phenos:
                        individual_phenos.append(p)

                for individual_p in individual_phenos:
                    how_text = ""
                    for dr in drug_rows.to_dict("records"):
                        dr_gene  = _safe(dr.get('Gene', ''))
                        dr_pheno = _safe(dr.get('Phenotype', ''))
                        if ":" in dr_pheno and dr_pheno.lower().startswith(dr_gene.lower()):
                            dr_pheno = dr_pheno.split(":", 1)[1].strip()
                        if dr_gene.lower() == g.lower() and dr_pheno.lower() == individual_p.lower():
                            how_text = _safe(dr.get('How this gene/phenotype affects the drug and what it means for you', ''))
                            if how_text: break

                    if not how_text:
                        for dr in drug_rows.to_dict("records"):
                            if _safe(dr.get('Gene', '')).lower() == g.lower():
                                how_text = _safe(dr.get('How this gene/phenotype affects the drug and what it means for you', ''))
                                if how_text: break

                    # Clinical Override: Strip "No RYR1 genotype is available. " if RYR1 is negative
                    if g.upper() == 'CACNA1S' and ryr1_status == 'negative':
                        how_text = re.sub(r"^no\s+ryr1\s+genotype\s+is\s+available\.\s*", "", how_text, flags=re.IGNORECASE)

                    how_text = how_text or 'No specific impact description available in the current database for this genotype.'

                    # Prevent 'refer to guideline' placeholder if gene truly has no recommendation
                    if "refer to guideline" in how_text.lower() or "refer to recommendation" in how_text.lower():
                        gene_has_rec = False
                        for dr_rec in drug_rows.to_dict("records"):
                            rg = _safe(dr_rec.get('Gene', ''))
                            rt = _safe(dr_rec.get('Recommendation', ''))
                            if rg.lower() == g.lower() and rt and not _is_placeholder_rec(rt):
                                gene_has_rec = True
                                break
                        if not gene_has_rec:
                            how_text = "This gene influences dosing algorithms, but standalone guidelines are not provided below. It must be calculated alongside your other relevant genes by your healthcare provider."

                    # Group identical paragraphs by a normalized text key to
                    # collapse near-duplicate explanatory paragraphs (e.g.,
                    # CYP4F2 unphased variants rendered slightly differently).
                    norm_key = re.sub(r"\b" + re.escape(g) + r"\b", "", how_text, flags=re.IGNORECASE)
                    norm_key = re.sub(r"\s+", " ", norm_key).strip().lower()
                    how_map.setdefault(norm_key, {'text': how_text, 'phenos': []})
                    how_map[norm_key]['phenos'].append(individual_p)

            # Build consolidated sub_blocks: list phenotypes/diplotypes then the single paragraph
            sub_blocks = ""
            for k, v in how_map.items():
                how_text = v['text']
                phenos = v['phenos']
                # Preserve original order and deduplicate phenotypes
                seen = []
                for p in phenos:
                    if p not in seen:
                        seen.append(p)
                phenos_str = ", ".join(seen)
                label_display = f"<strong>If {phenos_str}:</strong> " if is_unphased else ""
                sub_blocks += f"<div style='margin-bottom:6px;'><span style='color:#333;'>{label_display}{how_text}</span></div>"
            
            if sub_blocks not in impact_groups:
                impact_groups[sub_blocks] = []
                ordered_impacts.append(sub_blocks)
            impact_groups[sub_blocks].append(g)

    # Build the final HTML by grouping shared texts under combined headers
    for block_html in ordered_impacts:
        genes = impact_groups[block_html]
        genes_str = " / ".join(genes)
        if multi_gene:
            impact_sections_html += f"<div style='margin-bottom:5px;'><strong style='color:#0D3B7A;'>{genes_str}</strong><br>{block_html}</div>"
        else:
            impact_sections_html += f"<div style='margin-bottom:5px;'>{block_html}</div>"

    if not impact_sections_html:
        impact_sections_html = f'<div style="font-size:{fs_body}; color:#6b7280; font-style:italic;">Refer to guideline recommendations below.</div>'

    unphased_note_global = ""
    if has_unphased:
        unphased_note_global = f'<div style="font-size:11.5px; color:#c2410c; font-style:italic; margin-bottom:8px;">&#9888; Phasing was not performed for one or more genes. Multiple possible interpretations are provided below.</div>'

    impact_block = (
        f'<div class="info-box" style="margin-bottom:10px; border: 1px solid #d0e4f7; border-left: 4px solid #4DB7D0; border-radius: 8px; background: #f4f9ff; padding: 8px 14px; page-break-inside: auto; break-inside: auto;">'
        f'<div style="font-size:{fs_headers}; font-weight:700; color:#1a1a1a; margin-bottom:6px;">What it means for you</div>'
        f'{unphased_note_global}'
        f'{impact_sections_html}'
        f'</div>'
    )

    # --- 4. RECOMMENDATIONS ---
    ev_parts = []
    for label, col in [('CPIC', 'CPIC Level'), ('PharmGKB', 'PharmGKB LoE'), ('FDA', 'PGx on FDA Label')]:
        for r in drug_rows.to_dict("records"):
            v = _safe(r.get(col, ''))
            if v:
                ev_parts.append(f'<span style="background:#f3f4f6; border:1px solid #e5e7eb; color:#4b5563; padding:2px 6px; border-radius:4px; font-size:12px; font-weight:600; margin-left:4px;">{label}: {v}</span>')
                break
    ev_html = ''.join(ev_parts)


    # ── Source Status Summary bar ─────────────────────────────────────────────
    # Shows a compact badge for each source that evaluated this drug-gene combo.
    # Gives the patient/doctor clear visibility into WHY a source has no guidance.
    _STATUS_BADGE_STYLES = {
        "dosing info":     {"bg": "#dcfce7", "color": "#166534", "border": "#86efac"},
        "alternate drug":  {"bg": "#dcfce7", "color": "#166534", "border": "#86efac"},
        "other guidance":  {"bg": "#dcfce7", "color": "#166534", "border": "#86efac"},
        "no action":       {"bg": "#dbeafe", "color": "#1e40af", "border": "#93c5fd"},
        "no recommendation": {"bg": "#f3f4f6", "color": "#4b5563", "border": "#d1d5db"},
    }

    def _status_badge(source_label, status_val):
        if not status_val:
            return ""
        s_lower = str(status_val).lower()
        style = next(
            (st for key, st in _STATUS_BADGE_STYLES.items() if key in s_lower),
            None,
        )
        if style is None:
            return ""
        return (
            f'<span style="background:{style["bg"]}; color:{style["color"]}; '
            f'border:1px solid {style["border"]}; padding:2px 6px; border-radius:4px; '
            f'font-size:12px; font-weight:600; white-space:nowrap;">'
            f'{source_label}: {status_val}</span>'
        )

    _src_statuses = {}
    for _src in ['CPIC', 'DPWG', 'FDA_Label', 'FDA_Table']:
        for _r in drug_rows.to_dict("records"):
            _sv = _safe(_r.get(f'{_src} Status', ''))
            if _sv:
                _src_statuses[_src] = _sv
                break

    # Source display labels → column prefix mapping (GSI uses underscores)
    source_styles = {
        'CPIC':       {'text': '#1e40af', 'bg': '#eff6ff', 'border': '#bfdbfe'},
        'DPWG':       {'text': '#0f766e', 'bg': '#f0fdfa', 'border': '#99f6e4'},
        'FDA_Label':  {'text': '#c2410c', 'bg': '#fff7ed', 'border': '#fed7aa'},
        'FDA_Table':  {'text': '#c2410c', 'bg': '#fff7ed', 'border': '#fed7aa'},
    }
    source_display = {
        'CPIC':      'CPIC',
        'DPWG':      'DPWG',
        'FDA_Label': 'FDA Label',
        'FDA_Table': 'FDA Table',
        'PharmCAT':  'CPIC',  # PharmCAT fallback uses CPIC recommendations
    }

    source_status_bar = ""
    all_recs  = []
    seen_recs = set()

    # Loop through all sources independently—keeping 1, 2, 3, or all 4 guidelines if present
    for source in ['CPIC', 'DPWG', 'FDA_Label', 'FDA_Table']:
        for r in drug_rows.to_dict("records"):
            # GSI schema: single Details column; fallback to old column names for compatibility
            details = _safe(r.get(f'{source} Details', ''))
            rec     = _safe(r.get(f'{source} Recommendation', ''))   # PharmCAT / legacy
            impl    = _safe(r.get(f'{source} Implications', ''))      # PharmCAT / legacy
            other   = _safe(r.get(f'{source} Other Considerations', ''))  # legacy
            annotation = _safe(r.get(f'{source} AnnotationLink', ''))

            # Merge: if GSI Details present, use it as the recommendation body
            if details and not rec:
                rec = details
            if not any([rec, impl, other, details]):
                continue

            g = _safe(r.get('Gene', ''))

            # Clinical Override: Skip CACNA1S recommendations if RYR1 is positive (MHS)
            if g.upper() == 'CACNA1S' and ryr1_status == 'positive':
                continue

            # Clinical Override: Strip "No RYR1 genotype is available. " if RYR1 is negative
            if g.upper() == 'CACNA1S' and ryr1_status == 'negative':
                rec = re.sub(r"^no\s+ryr1\s+genotype\s+is\s+available\.\s*", "", rec, flags=re.IGNORECASE)
                details = re.sub(r"^no\s+ryr1\s+genotype\s+is\s+available\.\s*", "", details, flags=re.IGNORECASE)
                impl = re.sub(r"^no\s+ryr1\s+genotype\s+is\s+available\.\s*", "", impl, flags=re.IGNORECASE)
            p = _safe(r.get('Phenotype', ''))
            if ":" in p and p.lower().startswith(g.lower()):
                p = p.split(":", 1)[1].strip()

            # Format phenotype: capitalize metabolizer phenotypes and convert (as:X) to (AS:X)
            p = re.sub(r'\(as:([^)]+)\)', lambda m: f'(AS:{m.group(1)})', str(p))
            p = re.sub(r'\bnormal\s+metabolizer\b', 'Normal Metabolizer', p, flags=re.IGNORECASE)
            p = re.sub(r'\bpoor\s+metabolizer\b', 'Poor Metabolizer', p, flags=re.IGNORECASE)
            p = re.sub(r'\bintermediate\s+metabolizer\b', 'Intermediate Metabolizer', p, flags=re.IGNORECASE)
            p = re.sub(r'\bultrarapid\s+metabolizer\b', 'Ultrarapid Metabolizer', p, flags=re.IGNORECASE)
            p = re.sub(r'\brapid\s+metabolizer\b', 'Rapid Metabolizer', p, flags=re.IGNORECASE)
            p = re.sub(r'\bextensive\s+metabolizer\b', 'Extensive Metabolizer', p, flags=re.IGNORECASE)

            # ── FILTER 0: Skip genes that are in the yellow box ────────────────
            # These genes are "known to affect the drug" but PharmCAT has no
            # specific guidance for them. They are displayed in the yellow box
            # below with explanatory text. Skip them here to avoid redundant
            # "no guidance" cards in the Clinical Recommendations section.
            if g and _normalise_gene(g) in genes_catalog_only:
                continue

            # ── FILTER 1: Skip genes that require specialised testing ─────────
            # If this patient's phenotype for the gene is Indeterminate / No Call,
            # we cannot personalise the recommendation.  The drug should rely on
            # other called genes, OR appear in the Specialized Testing section.
            if g and g.lower() in _spec_testing_genes_lower:
                continue

            # ── FILTER 2: Only show recommendations for the patient's phenotype ─
            # GSI rows include guidance for every phenotype (Normal, Intermediate,
            # Poor, etc.).  Without this filter, a single drug page would show
            # CPIC text for ALL phenotypes — confusing the reader.  Skip rows
            # whose phenotype doesn't match the patient's call for that gene.
            _patient_phens = _patient_phen_per_gene.get(g.lower(), set())
            _row_phen_norm = re.sub(r"\s*\(as:[0-9.]+\)", "", p.lower().strip())
            if _patient_phens and _row_phen_norm and _row_phen_norm not in _patient_phens:
                continue

            # ── FILTER 3: Suppress placeholder "X has no recommendations for Y" ─
            # These add noise without conveying actionable information.  When
            # CPIC has no guidance for a particular gene-drug pair, the cleanest
            # presentation is silence (or a future muted "Coverage Notes" line).
            # EXCEPTION: Do NOT filter out FDA/DPWG placeholder entries, as they
            # may contain meaningful guidance when CPIC lacks actionable content.
            # Apply per-drug overrides: e.g. {'brivaracetam': {'suppress_sources': ['DPWG']}}
            drug_key_local = str(drug_name).strip().lower()
            suppress_sources = []
            if per_drug_overrides and isinstance(per_drug_overrides, dict):
                sup = per_drug_overrides.get(drug_key_local, {}).get('suppress_sources')
                if sup:
                    suppress_sources = [s.strip().lower() for s in sup]

            # If this source is configured to be suppressed for this drug, skip it
            if source_display.get(source, source).lower() in suppress_sources:
                continue

            if source not in ('FDA_Label', 'FDA_Table', 'DPWG'):
                if _is_placeholder_rec(rec) and _is_placeholder_rec(details):
                    continue
                if _is_placeholder_rec(rec) and not impl and not other:
                    continue

            context_str = f"{g} ({p})" if g and p else g or p or "General Recommendation"
            # Deduplicate by context+source: allow multiple sources for same gene-phenotype,
            # but skip if we already have THIS SOURCE for this context
            context_source_key = (context_str, source_display.get(source, source))
            if context_source_key in seen_recs:
                continue
            seen_recs.add(context_source_key)

            # Tag this entry as informational (muted card) vs actionable.
            # Informational = the source ran the lookup but has truly no guidance:
            #   • Status = "Not Evaluated" (e.g. atazanavir + CYP2C19 — "No CPIC guidance.")
            # Status values like "No Action" and "No Recommendation" frequently
            # carry substantive guidance (e.g. voriconazole's "Initiate therapy with
            # recommended standard of care dosing…", rosuvastatin's "Based on ABCG2
            # status, prescribe desired starting dose…").  Render those as normal
            # actionable cards.  Pure-placeholder text (e.g. "CPIC has no
            # recommendations for atorvastatin and ABCG2") is already filtered
            # upstream by _is_placeholder_rec.
            _status_val = _safe(r.get(f'{source} Status', '')).lower()
            _is_info = 'not evaluated' in _status_val

            all_recs.append({
                'source':   source_display.get(source, source),
                'style':    source_styles[source],
                'context':  context_str,
                'impl':     impl,
                'rec':      rec,
                'other':    other,
                'is_info':  _is_info,
                'status':   _status_val,
            })

    # ── Merge FDA_Label / FDA_Table into a single 'FDA' source ──────────────
    # User requirement: FDA Label and FDA Table are displayed as one "FDA" row.
    # If both sources have content for the same gene-phenotype context:
    # - Prefer FDA_Table if FDA_Label status is "No Action" (not actionable)
    # - Otherwise use FDA_Label (appears first in loop)
    # Entries for distinct contexts are all kept.
    _fda_merged_style = {'text': '#c2410c', 'bg': '#fff7ed', 'border': '#fed7aa'}
    _merged_all_recs  = []
    _fda_by_ctx       = {}  # Store FDA entries by context, preferring actionable ones

    # First pass: collect FDA entries, preferring FDA_Table over non-actionable FDA_Label
    for _item in all_recs:
        if _item['source'] in ('FDA Label', 'FDA Table'):
            _ctx = _item['context']
            if _ctx not in _fda_by_ctx:
                _fda_by_ctx[_ctx] = _item
            else:
                # Both FDA_Label and FDA_Table exist for this context
                # Prefer FDA_Table if FDA_Label status is "No Action"
                _existing = _fda_by_ctx[_ctx]
                _existing_status = _existing.get('status', '').lower()
                if _existing['source'] == 'FDA Label' and 'no action' in _existing_status:
                    _fda_by_ctx[_ctx] = _item

    # Second pass: build merged records
    _fda_ctx_added = set()
    for _item in all_recs:
        if _item['source'] in ('FDA Label', 'FDA Table'):
            _ctx = _item['context']
            if _ctx not in _fda_ctx_added:
                _fda_ctx_added.add(_ctx)
                _merged = {**_fda_by_ctx[_ctx], 'source': 'FDA', 'style': _fda_merged_style}
                _merged_all_recs.append(_merged)
            # else: already added the best FDA entry for this context
        else:
            _merged_all_recs.append(_item)
    all_recs = _merged_all_recs

    # ── Source Priority Filter: CPIC > DPWG/FDA ───────────────────────────────
    # If CPIC guidelines exist with ACTIONABLE guidance (not just "not evaluated"),
    # show ONLY CPIC recommendations. If CPIC has no actionable content (all entries
    # marked is_info=True, meaning "Gene analyzed, but no guidance available"), then
    # fall back to DPWG and/or FDA (whichever is available), regardless of their
    # is_info status. This ensures FDA/DPWG recommendations are shown when CPIC has
    # no guidance, even if FDA/DPWG entries are marked as informational.

    _has_cpic_actionable = any(
        r['source'] == 'CPIC' and not r.get('is_info', False)
        for r in all_recs
    )
    if _has_cpic_actionable:
        # Filter: keep only CPIC and non-guideline sources (PharmCAT)
        # Note: this also removes CPIC is_info cards so they don't clutter the page
        all_recs = [r for r in all_recs if r['source'] in ('CPIC', 'PharmCAT')]
    else:
        # No actionable CPIC: filter to only DPWG/FDA recommendations (any status)
        # This ensures FDA/DPWG guidance is shown when CPIC has no recommendations.
        _has_dpwg_fda = any(
            r['source'] in ('DPWG', 'FDA')
            for r in all_recs
        )
        if _has_dpwg_fda:
            all_recs = [r for r in all_recs if r['source'] in ('DPWG', 'FDA')]
        # else: keep all (DPWG/FDA empty, only PharmCAT remains)

    # Remove all "is_info" entries — these are placeholder cards that say
    # "Gene analyzed, but no current clinical guidance available" and add
    # no value to the report.  They just create noise.
    all_recs = [r for r in all_recs if not r.get('is_info', False)]

    # PharmCAT fallback: if every GSI source returned no recommendation (e.g.
    # CPIC "No Recommendation" for DPYD+fluorouracil even though PharmCAT
    # embeds a CPIC-Level-A "Avoid" call), surface the PharmCAT Recommendation
    # text directly so the Clinical Recommendations box is never empty when
    # there IS actionable guidance.
    if not all_recs:
        _pharmcat_style = {'text': '#6b21a8', 'bg': '#faf5ff', 'border': '#d8b4fe'}
        for r in drug_rows.to_dict("records"):
            g = _safe(r.get('Gene', ''))
            if g and g.lower() in _spec_testing_genes_lower:
                continue
            _pc_rec  = _safe(r.get('Recommendation', ''))
            _pc_impl = _safe(r.get('Implication', ''))
            if not _pc_rec or _is_placeholder_rec(_pc_rec):
                continue
            p = _safe(r.get('Phenotype', ''))
            if ":" in p and p.lower().startswith(g.lower()):
                p = p.split(":", 1)[1].strip()
            context_str = f"{g} ({p})" if g and p else g or p or "General Recommendation"
            rec_key = f"PharmCAT_{context_str}_{_pc_rec[:30]}"
            if rec_key in seen_recs:
                continue
            seen_recs.add(rec_key)
            all_recs.append({
                'source':  'PharmCAT',
                'style':   _pharmcat_style,
                'context': context_str,
                'impl':    _pc_impl,
                'rec':     _pc_rec,
                'other':   '',
                'is_info': False,
                'status':  '',
            })

# ================================================================
# DYNAMIC DRUG-PAGE PAGINATION
#
# Keep all recommendations inside one logical drug section.
# Chromium decides the physical page breaks from the ACTUAL
# rendered height instead of recommendation count.
# ================================================================

    pages_html = []

    rec_chunks = (
        [all_recs]
        if all_recs
        else []
    )

    def _links(raw: str) -> str:
        return ' '.join(
            f'<a href="{u}" style="color:#1a56db; font-size:12px; word-break:break-all; text-decoration:none;">{u}</a>'
            for u in raw.split() if u.startswith('http')
        )

    urls, citations = "", ""
    for r in drug_rows.to_dict("records"):
        if not urls:      urls      = _safe(r.get('URLs', ''))
        if not citations: citations = _safe(r.get('citation links', ''))

    links_html = ""
    if urls or citations:
        parts = []
        if urls:      parts.append(f"<strong style='color:#374151;'>CPIC:</strong> {_links(urls)}")
        if citations: parts.append(f"<strong style='color:#374151;'>Citations:</strong> {_links(citations)}")
        links_html = (
            f'<div style="margin-top:2px; padding:4px 8px; font-size:12px; border:1px solid #e5e7eb; '
            f'border-radius:6px; background:#ffffff; color:#6b7280; line-height:1.5;">'
            f'{" &nbsp;&nbsp;|&nbsp;&nbsp; ".join(parts)}</div>'
        )

    # ── Build pages ──────────────────────────────────────────────────────────────
    # Identify spec-testing genes that affect THIS drug — used to pick a more
    # informative fallback message when no actionable recommendations were
    # produced.  Otherwise the user sees a misleading "Standard Prescribing
    # Apply" even when the drug's primary PGx-driver gene (e.g. SLCO1B1 for
    # statins) couldn't be reliably called.
    _spec_for_this_drug = [
        g for g in genes_to_evaluate
        if g and g.lower() in _spec_testing_genes_lower
    ]
    # Separate "needs specialized lab test" from "not on consumer platform"
    _spec_test_genes   = [g for g in _spec_for_this_drug
                          if gene_status_map.get(g) == "Requires Specialized Testing"]
    _platform_not_eval = [g for g in _spec_for_this_drug
                          if gene_status_map.get(g) == "Not Evaluated by Platform"]

    if not rec_chunks:
        subline = f"<span style='float:right;'>{ev_html}</span>" if ev_html else ""
        if _spec_test_genes:
            _gene_list = ", ".join(_spec_test_genes)
            _is_plural = len(_spec_test_genes) > 1
            _plat_note = ""
            if _platform_not_eval:
                _plat_genes = ", ".join(_platform_not_eval)
                _plat_note = (
                    f" Additionally, <strong>{_plat_genes}</strong> "
                    f"{'is' if len(_platform_not_eval) == 1 else 'are'} not evaluated "
                    f"by consumer DNA platforms and could not be assessed."
                )
            fallback_body = (
                f'<strong style="color:#92400e;">Personalised Recommendation Unavailable:</strong> '
                f'A definitive pharmacogenomic recommendation for {drug_name.title()} could not be '
                f'generated because '
                f'<strong>{_gene_list}</strong>{" requires" if not _is_plural else " require"} '
                f'specialised laboratory testing beyond standard DNA sequencing. '
                f'Other relevant genes that <em>were</em> analysed for this medication are '
                f'shown above with their interpreted phenotypes. '
                f'Please refer to the <em>Genes Requiring Specialized Testing</em> section for '
                f'details on the recommended follow-up assay. Standard CPIC dosing applies '
                f'pending those results.{_plat_note}'
            )
            _fb_bg = "#fffbeb"
            _fb_border = "#f59e0b"
        elif _platform_not_eval:
            # Only "Not Evaluated by Platform" genes — no specialized testing needed,
            # but the recommendation is incomplete because the platform doesn't cover all genes.
            _gene_list = ", ".join(_platform_not_eval)
            _is_plural = len(_platform_not_eval) > 1
            fallback_body = (
                f'<strong style="color:#4b5da8;">Partial Evaluation — Gene Not on Platform:</strong> '
                f'<strong>{_gene_list}</strong> '
                f'{"is" if not _is_plural else "are"} not evaluated by consumer DNA genotyping '
                f'platforms and could not be assessed in this analysis. '
                f'Other relevant genes that <em>were</em> analysed for this medication are '
                f'shown above with their interpreted phenotypes. '
                f'If a complete pharmacogenomic assessment is required, a dedicated clinical '
                f'laboratory test covering this gene may be requested through your healthcare provider.'
            )
            _fb_bg = "#f0f4ff"
            _fb_border = "#818cf8"
        else:
            fallback_body = (
                '<strong>Standard Prescribing Guidelines Apply:</strong> Based on the genetic '
                'variants analyzed, there are no specific pharmacogenomic dose adjustments '
                'required or available. Prescribe as directed by standard clinical guidelines, '
                'taking into account patient age, weight, organ function, and drug-drug '
                'interactions.'
            )
            _fb_bg = "#f9fafb"
            _fb_border = "#e2eaf4"

        fallback_rec = f"""
        <div style="margin-bottom: 14px;">
            <div style="margin-bottom: 10px; padding-left: 2px; overflow:hidden;">
                <span style="font-size: {fs_headers}; font-weight: 700; color: #1a1a1a;">Clinical Recommendations</span>
                {subline}
            </div>
            <div style="padding:12px 16px; border:1px solid {_fb_border}; border-left:4px solid {_fb_border}; border-radius:6px; background:{_fb_bg}; color:#4b5563; font-size:12px; line-height:1.55;">
                {fallback_body}
            </div>
        </div>
        """
        content = f"""<div style="padding-top: 10px;">{header_html}{class_box}{drug_title_box}{genes_block}{source_status_bar}{about_block}{impact_block}{fallback_rec}{links_html}{_extra_category_html}</div>"""
        pages_html.append(_wrap_page(content, patient_name, curr_pg, page_id=drug_id))
        curr_pg += 1
    else:
        for i, chunk in enumerate(rec_chunks):
            grouped_chunk = {}
            for item in chunk:
                if item['source'] not in grouped_chunk:
                    grouped_chunk[item['source']] = {'style': item['style'], 'items': []}
                grouped_chunk[item['source']]['items'].append(item)

            guidelines_block = ""
            if grouped_chunk:
                g_html = ""
                for src, data in grouped_chunk.items():
                    style      = data['style']
                    items_html = ""
                    for idx, item in enumerate(data['items']):
                        # Deduplicate identical recommendation blocks per-source
                        b_bottom = (
                            "border-bottom: 1.5px solid #cbd5e1; padding-bottom: 14px; margin-bottom: 14px;"
                            if idx < len(data['items']) - 1 else "padding-bottom: 4px;"
                        )
                        # Skip duplicate item text to avoid repeated paragraphs
                        # Phenotype/context is part of the identity.  Two patient-compatible
                        # phenotypes can legitimately share identical CPIC text; do not hide one
                        # merely because the recommendation wording is the same.
                        _item_key = '||'.join([str(item.get(k, '')).strip() for k in ('context', 'rec', 'impl', 'other')])
                        if not hasattr(grouped_chunk[src], '_seen_texts'):
                            grouped_chunk[src]['_seen_texts'] = set()
                        if _item_key in grouped_chunk[src]['_seen_texts']:
                            continue
                        grouped_chunk[src]['_seen_texts'].add(_item_key)
                        ctx_html = (
                            f"<div style='font-size:11.5px; font-weight:700; color:#111827; margin-bottom:8px;'>{item['context']}</div>"
                            if item['context'] else ""
                        )
                        if item.get('is_info'):
                            rec_html = (
                                f"<div style='background-color:#f9fafb; border:1px dashed #d1d5db; "
                                f"border-radius:6px; padding:7px 10px; font-size:11.5px; color:#6b7280; "
                                f"font-style:italic; line-height:1.5;'>"
                                f"<strong style='color:#4b5563; font-style:normal;'>{src}:</strong> "
                                f"Gene analyzed, but no current clinical guidance available for this specific drug-gene pair."
                                f"</div>"
                            )
                            impl_html, other_html = "", ""
                        else:
                            # Format the rec text with bold keyword headings and per-
                            # subgroup splits when the GSI text contains them.
                            _formatted_rec, _is_structured = _format_rec_text(item['rec'])
                            if _is_structured:
                                # Formatter has already added bold "Recommendation:" /
                                # "Implications:" / "Other Considerations:" headings.
                                # Wrap in the same card style as before.
                                rec_html = (
                                    f"<div style='background-color:#f8fafc; "
                                    f"border:1px solid #e2e8f0; border-radius:6px; "
                                    f"padding:10px 12px; font-size:12px; color:#0f172a; "
                                    f"margin-bottom:8px; line-height:1.55;'>"
                                    f"{_formatted_rec}</div>"
                                ) if _formatted_rec else ""
                                # Suppress legacy impl/other slots — the formatter
                                # already rendered them inline from the GSI text.
                                impl_html, other_html = "", ""
                            else:
                                # Plain unstructured text — keep the original layout
                                # so legacy PharmCAT-only data (separate impl/other
                                # columns) still render correctly.
                                impl_html  = f"<div style='font-size:12px; color:#4b5563; margin-bottom:8px; line-height:1.5;'><strong style='color:#374151;'>Implication:</strong> {item['impl']}</div>" if item['impl'] else ""
                                rec_html   = f"<div style='background-color:#f8fafc; border:1px solid #e2e8f0; border-radius:6px; padding:8px 10px; font-size:12px; color:#0f172a; margin-bottom:8px; line-height:1.5;'><strong style='color:#0f172a;'>Recommendation:</strong> {_formatted_rec}</div>" if _formatted_rec else ""
                                other_html = f"<div style='font-size:11.5px; color:#6b7280; line-height:1.4;'><strong style='color:#4b5563;'>Other Considerations:</strong> {item['other']}</div>" if item['other'] else ""
                        items_html += f"""
                        <div class="recommendation-entry" style="{b_bottom} page-break-inside:avoid; break-inside:avoid;">
                            {ctx_html}{impl_html}{rec_html}{other_html}
                        </div>
                        """
                    # Tightened spacing: outer margin 14→4px, padding 12 16→6 12,
                    # source-tag bottom margin 12→4px.  Reduces wasted vertical
                    # space between recommendation cards by ~30–40% so multi-gene
                    # drugs (Atorvastatin, Fluvastatin, Lovastatin) fit on one
                    # page instead of bleeding to a near-empty continuation page.
                    # NOTE: page-break-inside:avoid REMOVED from the outer source
                    # container — for multi-gene drugs (warfarin, clopidogrel) the
                    # entire source block can be too tall, causing huge gaps when
                    # the browser pushes it to the next page.
                    g_html += f"""
                    <div style="margin-bottom:4px; padding:6px 12px 1px 12px; border-left:4px solid {style['text']}; border-radius:6px; background:#ffffff; box-shadow: 0 1px 2px rgba(0,0,0,0.04); border-top:1px solid #f3f4f6; border-right:1px solid #f3f4f6; border-bottom:1px solid #f3f4f6;">
                        <div class="recommendation-source-label" style="margin-bottom:4px; page-break-after:avoid; break-after:avoid;">
                            <span style="background-color:{style['bg']}; color:{style['text']}; border:1px solid {style['border']}; padding:2px 6px; border-radius:12px; font-size:12px; font-weight:700; letter-spacing:0.3px; text-transform:uppercase;">{src}</span>
                        </div>
                        {items_html}
                    </div>"""

                subline = f"<span style='float:right;'>{ev_html}</span>" if ev_html else ""
                
                # --- WARFARIN CUSTOM OVERRIDE ---
                warfarin_extra = ""
                if drug_name.strip().lower() == "warfarin":
                    warfarin_extra = f"""
                    <div style="margin-bottom:14px; padding:10px 14px; border-left:4px solid #1e40af; border-radius:6px; background:#eff6ff; border-top:1px solid #bfdbfe; border-right:1px solid #bfdbfe; border-bottom:1px solid #bfdbfe;">
                        <div style="font-size:12px; font-weight:700; color:#1e3a5f; margin-bottom:4px;">Warfarin Dosing Algorithm</div>
                        <div style="font-size:{fs_body}; color:#333; line-height:1.5;">
                            The CPIC guidelines reference a specific dosing algorithm to calculate Warfarin dosage based on CYP2C9 and VKORC1 variants, as well as patient age, weight, and clinical factors. Doctors can access the validated, interactive dosing calculator directly at <strong><a href="http://www.warfarindosing.org" style="color:#1a56db; text-decoration:none;">WarfarinDosing.org</a></strong>.
                        </div>
                    </div>
                    """

                guidelines_block = f"""
                <div style="margin-bottom: 8px;">
                    <div style="margin-bottom: 6px; padding-left: 2px; overflow:hidden; page-break-after: avoid; break-after: avoid;">
                        <span style="font-size: {fs_headers}; font-weight: 700; color: #1a1a1a;">Clinical Recommendations</span>
                        {subline}
                    </div>
                    {g_html}
                    {warfarin_extra}
                </div>
                """

                if i == 0:
                    content = f"""<div style="padding-top: 8px;">{header_html}{class_box}{drug_title_box}{genes_block}{source_status_bar}{about_block}{impact_block}{guidelines_block}{links_html if len(rec_chunks) == 1 else ""}{_extra_category_html if len(rec_chunks) == 1 else ""}</div>"""
                    pages_html.append(_wrap_page(content, patient_name, curr_pg, page_id=drug_id))
                else:
                    cont_title_box = f"""<div class="drug-header" id="{drug_id}_cont_{i}" style="background: linear-gradient(90deg, #0984b6 0%, #63c0d3 100%); padding: 7px 12px; border-radius: 12px; margin-bottom: 6px;"><h3 class="drug-title" style="font-size: {fs_subheaders}; font-weight: 700; margin: 0; color: #ffffff; display:inline-block;"><span style="text-transform:uppercase;">{drug_name[0]}</span>{drug_name[1:]} (Continued)</h3></div>"""
                    content = f"""<div style="padding-top: 8px;">{cont_title_box}{guidelines_block}{links_html if i == len(rec_chunks) - 1 else ""}{_extra_category_html if i == len(rec_chunks) - 1 else ""}</div>"""
                    pages_html.append(_wrap_page(content, patient_name, curr_pg, page_id=f"{drug_id}_{i}"))
                curr_pg += 1
        return pages_html, curr_pg


def other_evaluated_medicines_template(
    df,
    name,
    pg,
    master_genes_df,
    drug_gene_map
):
    """
    Other Evaluated Medications

    Groups medications by:

        Gene + patient's genetic result

    instead of showing one table row per medication.

    Important:
    A medication associated with more than one gene can appear in
    more than one gene card. This is intentional because each
    gene-drug relationship was evaluated separately.
    """

    # ========================================================================
    # Helpers
    # ========================================================================

    def _safe(val):
        s = str(val).strip()

        return (
            s
            if s.lower() not in (
                '',
                'nan',
                'none',
                'n/a'
            )
            else ''
        )


    def _is_unknown_result(value):
        """
        Return True when the phenotype/diplotype does not represent
        a usable patient genetic result.
        """

        v = str(
            value or ''
        ).strip().lower()

        return v in {
            '',
            'nan',
            'none',
            'n/a',
            'unknown',
            'unknown/unknown',
            'unknown / unknown',
            'no result',
            'no data available',
            'not called',
            'unassigned',
            'uncategorized',
        }


    # ========================================================================
    # Empty OEM dataframe
    # ========================================================================

    if (
        df is None
        or df.empty
        or "Drug Name" not in df.columns
    ):

        content = """
        <div style="padding-top:10px;">

            <div
                class="page-title"
                id="other_evaluated"
                style="margin-bottom:14px;"
            >
                Other Evaluated Medications
            </div>

            <div style="
                padding:14px;
                border:1px solid #e5e7eb;
                border-radius:6px;
                background:#f9fafb;
                text-align:center;
                color:#6b7280;
                font-size:11.5px;
            ">
                No other evaluated medications.
            </div>

        </div>
        """

        return (
            [
                _wrap_page(
                    content,
                    name,
                    pg,
                    page_id="other_evaluated"
                )
            ],
            pg + 1
        )


    # ========================================================================
    # STEP 1
    #
    # Build one patient-result record per gene from master_genes_df.
    #
    # Example:
    #
    # CYP2C9 ->
    #     diplotype = *1/*11
    #     phenotype = Intermediate Metabolizer
    #
    # CYP2D6 ->
    #     diplotype = Unknown/Unknown
    #     phenotype = No data available
    # ========================================================================

    gene_result_map = {}


    if (
        master_genes_df is not None
        and not master_genes_df.empty
        and "Gene" in master_genes_df.columns
    ):

        master_df = (
            master_genes_df.copy()
        )

        master_df["Gene"] = (
            master_df["Gene"]
            .astype(str)
            .str.strip()
        )


        for gene in (
            master_df["Gene"]
            .dropna()
            .unique()
        ):

            gene = str(
                gene
            ).strip()

            if not gene:
                continue


            matches = master_df[
                master_df["Gene"]
                .str.lower()
                == gene.lower()
            ]


            diplotypes = []
            phenotypes = []


            for row in (
                matches
                .to_dict("records")
            ):

                raw_d = _safe(
                    row.get(
                        "Diplotype",
                        ""
                    )
                )

                raw_p = _safe(
                    row.get(
                        "Phenotype",
                        ""
                    )
                )


                # Some master files contain:
                #
                # CYP2C9: Intermediate Metabolizer
                #
                # Strip the gene prefix.
                if (
                    ":"
                    in raw_p
                    and raw_p.lower().startswith(
                        gene.lower()
                    )
                ):

                    raw_p = (
                        raw_p
                        .split(
                            ":",
                            1
                        )[1]
                        .strip()
                    )


                if (
                    raw_d
                    and raw_d
                    not in diplotypes
                ):
                    diplotypes.append(
                        raw_d
                    )


                if (
                    raw_p
                    and raw_p
                    not in phenotypes
                ):
                    phenotypes.append(
                        raw_p
                    )


            # ------------------------------------------------------------
            # Remove "no result" phenotype if another usable phenotype
            # exists for the same gene.
            # ------------------------------------------------------------

            usable_phenotypes = [
                p
                for p in phenotypes
                if not _is_unknown_result(p)
            ]


            if usable_phenotypes:

                phenotypes = (
                    usable_phenotypes
                )


            # ------------------------------------------------------------
            # Collapse multiple possible diplotypes / phenotypes.
            #
            # This prevents one uncertain gene from generating several
            # duplicate OEM cards.
            # ------------------------------------------------------------

            if diplotypes:

                diplotype_display = (
                    " / ".join(
                        diplotypes
                    )
                )

            else:

                diplotype_display = (
                    "Unknown"
                )


            if phenotypes:

                phenotype_display = (
                    " / ".join(
                        phenotypes
                    )
                )

            else:

                phenotype_display = (
                    "No data available"
                )


            # ------------------------------------------------------------
            # A gene is considered uncallable if we have no meaningful
            # phenotype.
            # ------------------------------------------------------------

            is_uncallable = (
                not usable_phenotypes
            )


            gene_result_map[
                gene.upper()
            ] = {

                "gene":
                    gene,

                "diplotype":
                    diplotype_display,

                "phenotype":
                    phenotype_display,

                "is_uncallable":
                    is_uncallable,
            }


    # ========================================================================
    # STEP 2
    #
    # Determine which genes belong to each OEM drug.
    #
    # We preserve your CURRENT behavior:
    #
    #     genes from drug_gene_map
    #       +
    #     genes present directly in the OEM rows
    #
    # ========================================================================

    grouped = {}


    unique_drugs = sorted(
        df["Drug Name"]
        .dropna()
        .astype(str)
        .unique(),
        key=lambda x: x.lower()
    )


    for drug in unique_drugs:

        drug_name = (
            str(drug)
            .strip()
        )

        drug_key = (
            drug_name.lower()
        )


        drug_rows = df[
            df["Drug Name"]
            .astype(str)
            .str.strip()
            .str.lower()
            == drug_key
        ]


        genes_for_drug = set()


        # ------------------------------------------------------------
        # Genes from your Step 3 / Step 5 drug-gene map
        # ------------------------------------------------------------

        if drug_gene_map:

            for gene in (
                drug_gene_map.get(
                    drug_key,
                    []
                )
            ):

                gene = str(
                    gene
                ).strip()

                if gene:

                    genes_for_drug.add(
                        gene
                    )


        # ------------------------------------------------------------
        # Genes stored directly in OEM dataframe rows
        # ------------------------------------------------------------

        if "Gene" in drug_rows.columns:

            for row in (
                drug_rows
                .to_dict("records")
            ):

                raw_gene = _safe(
                    row.get(
                        "Gene",
                        ""
                    )
                )


                if not raw_gene:
                    continue


                # Supports:
                #
                # CYP2C9
                # CYP2C9; G6PD
                # CYP2C9, G6PD
                # CYP2C9\nG6PD

                gene_parts = re.split(
                    r'[;,\n]+',
                    raw_gene
                )


                for gene in gene_parts:

                    gene = (
                        gene.strip()
                    )

                    if gene:

                        genes_for_drug.add(
                            gene
                        )


        # ------------------------------------------------------------
        # Drug with no identifiable gene
        # ------------------------------------------------------------

        if not genes_for_drug:

            genes_for_drug.add(
                "Gene Not Available"
            )


        # ====================================================================
        # STEP 3
        #
        # Add this medication underneath EACH relevant gene/result.
        #
        # Multi-gene drugs intentionally appear in multiple cards.
        # ====================================================================

        for gene in sorted(
            genes_for_drug,
            key=lambda x: x.upper()
        ):

            gene_upper = (
                gene.upper()
            )


            result = (
                gene_result_map.get(
                    gene_upper
                )
            )


            if result is None:

                result = {

                    "gene":
                        gene,

                    "diplotype":
                        "Unknown",

                    "phenotype":
                        "No data available",

                    "is_uncallable":
                        True,
                }


            # ------------------------------------------------------------
            # Gene + patient result is the grouping key.
            # ------------------------------------------------------------

            group_key = (

                result["gene"].upper(),

                result["diplotype"],

                result["phenotype"],

                result["is_uncallable"],
            )


            if group_key not in grouped:

                grouped[
                    group_key
                ] = {

                    "gene":
                        result["gene"],

                    "diplotype":
                        result["diplotype"],

                    "phenotype":
                        result["phenotype"],

                    "is_uncallable":
                        result["is_uncallable"],

                    "drugs":
                        [],
                }


            display_drug = (
                drug_name.title()
            )


            if (
                display_drug
                not in grouped[
                    group_key
                ]["drugs"]
            ):

                grouped[
                    group_key
                ]["drugs"].append(
                    display_drug
                )


    # ========================================================================
    # STEP 4
    #
    # Convert grouped data into cards.
    #
    # Larger medication groups appear first.
    # ========================================================================

    cards = list(
        grouped.values()
    )


    cards.sort(
        key=lambda card: (
            -len(
                card["drugs"]
            ),
            card["gene"].upper(),
        )
    )


    # ========================================================================
    # Card renderer
    # ========================================================================

    def _render_card(card):

        gene = (
            card["gene"]
        )

        diplotype = (
            card["diplotype"]
        )

        phenotype = (
            card["phenotype"]
        )

        medications = (
            card["drugs"]
        )


        is_uncallable = (
            card["is_uncallable"]
        )


        # ------------------------------------------------------------
        # Header
        # ------------------------------------------------------------

        if is_uncallable:

            # Example:
            #
            # CYP2D6 — Unknown / Unknown

            header_text = (
                f'{gene} — '
                f'{diplotype}'
            )


            # Don't repeat:
            #
            # Unknown — No data available
            #
            # The explanatory note underneath is enough.

            subtitle_html = (
                '<div style="'
                'font-size:10.5px; '
                'font-style:italic; '
                'color:#6b7280; '
                'line-height:1.4; '
                'margin-top:2px; '
                'margin-bottom:5px;'
                '">'
                'This gene could not be confidently determined '
                'from your uploaded raw DNA file.'
                '</div>'
            )


            background = (
                "#fff8e8"
            )

            border_color = (
                "#d99a2b"
            )

            header_color = (
                "#92400e"
            )


        else:

            # Example:
            #
            # G6PD — B(reference)/B(reference) — Normal

            header_text = (
                f'{gene} — '
                f'{diplotype}'
            )


            if (
                phenotype
                and not _is_unknown_result(
                    phenotype
                )
            ):

                header_text += (
                    f' — {phenotype}'
                )


            subtitle_html = ""


            background = (
                "#eef6fc"
            )

            border_color = (
                "#2f86bd"
            )

            header_color = (
                "#023D79"
            )


        # Keep the compact visual design, but avoid handing Paged.js one giant
        # inline fragment.  Split the comma-separated list into short blocks that
        # can move together to the next page if needed, while preserving the
        # exact same appearance in the PDF.
        chunk_lists = _comma_safe_chunks(medications, max_chars=180)
        flat_items = [item for chunk in chunk_lists for item in chunk]
        meds_chunks = []
        flat_index = 0

        for chunk in chunk_lists:
            chunk_html = []
            for med in chunk:
                flat_index += 1
                suffix = ", " if flat_index < len(flat_items) else ""
                chunk_html.append(
                    f'<span class="pdf-qa-oem-medication oem-medication-item">'
                    f'{html_lib.escape(str(med))}'
                    f'</span>{suffix}'
                )
            meds_chunks.append(
                '<div class="oem-medication-chunk">' + ''.join(chunk_html) + '</div>'
            )

        meds_text = ''.join(meds_chunks)

        med_count = (
            len(
                medications
            )
        )


        count_text = (
            f'{med_count} medication'
            f'{"s" if med_count != 1 else ""} affected'
        )

        return f"""
        <div
            class="oem-card"
            style="
                background:{background};
                border-left:4px solid {border_color};
                border-radius:5px;
                padding:8px 12px 7px 12px;
                margin-bottom:8px;
            "
        >

            <div class="oem-card-head">

                <div style="
                    font-size:13px;
                    font-weight:700;
                    color:{header_color};
                    line-height:1.3;
                ">
                    {header_text}
                </div>

                {subtitle_html}

            </div>

            <div class="oem-medication-list" style="
                font-size:11.5px;
                color:#1f2937;
                line-height:1.4;
                margin-top:4px;
                margin-bottom:4px;
            ">
                {meds_text}
            </div>

            <div class="oem-card-count" style="
                font-size:10.5px;
                color:#6b7280;
                font-weight:500;
            ">
                {count_text}
            </div>

        </div>
        """

    # ========================================================================
    # STEP 5
    #
    # Let Chromium paginate using ACTUAL rendered card heights.
    # No hard-coded page budgets or card counts.
    # ========================================================================

    cards_html = "\n".join(
        _render_card(card)
        for card in cards
    )

    content = f"""
    <div
        class="flow-section"
        style="padding-top:10px;"
    >

        <div
            class="page-title"
            id="other_evaluated"
            style="margin-bottom:10px;"
        >
            Other Evaluated Medications
        </div>

        <p style="
            font-size:11.5px;
            color:#555;
            line-height:1.5;
            margin-bottom:14px;
        ">
            These medications were evaluated using your genetic
            profile but do not have a patient-specific clinical
            recommendation in this report. To make this section
            easier to review, medications are grouped below by
            the gene and genetic result involved in their evaluation.
        </p>

        {cards_html}

    </div>
    """

    return (
        [
            _wrap_page(
                content,
                name,
                pg,
                page_id="other_evaluated"
            )
        ],
        pg + 1
    )

# ============================================================================
# COMPLETE MEDICATION PANEL
# ============================================================================

def panel_coverage_template(all_drugs_df, detail_drug_pages, name, pg):
    """
    Renders the "Complete Medication Panel" section — a categorised list of
    EVERY drug associated with your genes that had a meaningful
    evaluation status in the clinical database.

    Parameters
    ----------
    all_drugs_df : pd.DataFrame
        Columns: Drug Name, Drug Category, Genes, About this medication,
                 Best Status, Has Guidance
    detail_drug_pages : dict
        {drug_name_lower: page_number} for drugs that have full detail pages.
    name : str  — patient name (for page footer)
    pg   : int  — starting page number
    """
    import math

    total_drugs = len(all_drugs_df)
    total_cats  = all_drugs_df["Drug Category"].nunique()

    # ── Badge helpers ──────────────────────────────────────────────────────────
    def _status_badge(row):
        drug_lower = str(row.get("Drug Name", "")).strip().lower()
        has_guid   = bool(row.get("Has Guidance", False))
        status     = str(row.get("Best Status", "")).strip()

        if has_guid:
            detail_pg = detail_drug_pages.get(drug_lower)
            if detail_pg:
                return (
                    f'<span style="background:#E7F7FB; color:#023D79; padding:2px 7px; '
                    f'border-radius:4px; font-size:12px; font-weight:700; white-space:nowrap;">'
                    f'&#128269; Detailed Report &rarr; p.{detail_pg}</span>'
                )
            return (
                '<span style="background:#E7F7FB; color:#023D79; padding:2px 7px; '
                'border-radius:4px; font-size:12px; font-weight:700;">&#128269; Detailed Report Available</span>'
            )
        elif status == "No Specific Action":
            return (
                '<span style="background:#e8f5e9; color:#1e7e34; padding:2px 7px; '
                'border-radius:4px; font-size:12px; font-weight:700; white-space:nowrap;">'
                '&#10003; No Action Needed</span>'
            )
        elif status == "No Action":
            return (
                '<span style="background:#f3f4f6; color:#6b7280; padding:2px 7px; '
                'border-radius:4px; font-size:12px; font-weight:600; white-space:nowrap;">'
                '&#9711; Not Gene-Dependent</span>'
            )
        elif status == "Gene Not Available":
            return (
                '<span style="background:#fef9c3; color:#92400e; padding:2px 7px; '
                'border-radius:4px; font-size:12px; font-weight:600; white-space:nowrap;">'
                '&#9888; Gene Not Analyzed</span>'
            )
        return (
            '<span style="background:#f3f4f6; color:#6b7280; padding:2px 7px; '
            'border-radius:4px; font-size:12px; white-space:nowrap;">—</span>'
        )

    # ── Build table rows, grouping by Drug Category ───────────────────────────
    cat_colors = ['#e3f2fd', '#fce4ec', '#e8f5e9', '#fff8e1', '#f3e5f5',
                  '#e0f7fa', '#fbe9e7', '#e8eaf6', '#f9fbe7', '#e0f2f1']
    cat_idx = 0
    rows_html = ""

    for cat, cat_grp in all_drugs_df.groupby("Drug Category"):
        bg = cat_colors[cat_idx % len(cat_colors)]
        cat_idx += 1
        n_drugs    = len(cat_grp)
        n_guidance = int(cat_grp["Has Guidance"].sum())
        rows_html += f"""
        <tr>
            <td colspan="3" style="padding:7px 12px; background:{bg};
                border-left:3px solid #4DB7D0; font-size:11.5px; font-weight:700;
                color:#023D79; letter-spacing:0.3px; text-transform:uppercase;">
                {cat}
                <span style="font-size:12px; font-weight:400; color:#555; margin-left:8px;">
                    {n_drugs} medication{'s' if n_drugs != 1 else ''}
                    {f'&nbsp;&bull;&nbsp;<strong style="color:#023D79">{n_guidance} with detailed report</strong>' if n_guidance else ''}
                </span>
            </td>
        </tr>"""

        for row in cat_grp.to_dict("records"):
            drug_name  = str(row.get("Drug Name", "")).strip().title()
            genes_str  = str(row.get("Genes", "")).strip()
            badge_html = _status_badge(row)
            rows_html += f"""
        <tr style="border-bottom:1px solid #f0f4f8;">
            <td style="padding:5px 12px 5px 22px; font-size:11.5px; font-weight:600;
                color:#111827; vertical-align:middle;">{drug_name}</td>
            <td style="padding:5px 12px; font-size:12px; color:#6b7280;
                vertical-align:middle;">{genes_str}</td>
            <td style="padding:5px 12px; vertical-align:middle;">{badge_html}</td>
        </tr>"""

    # ── Legend ─────────────────────────────────────────────────────────────────
    legend = """
    <div style="display:flex; gap:18px; margin-bottom:12px; padding:7px 12px;
        background:#f8fafc; border-radius:6px; flex-wrap:wrap;">
        <span style="font-size:12px; color:#333;">
            <strong style="color:#023D79;">&#128269; Detailed Report Available</strong>
            &mdash; Specific guidance found for your genotype
        </span>
        <span style="font-size:12px; color:#333;">
            <strong style="color:#1e7e34;">&#10003; No Action Needed</strong>
            &mdash; Evaluated; your genotype requires no specific action
        </span>
        <span style="font-size:12px; color:#333;">
            <strong style="color:#6b7280;">&#9711; Not Gene-Dependent</strong>
            &mdash; This medication is not significantly affected by genetics
        </span>
        <span style="font-size:12px; color:#333;">
            <strong style="color:#92400e;">&#9888; Gene Not Analyzed</strong>
            &mdash; This gene was not called in your DNA file; testing may be needed
        </span>
    </div>"""

    content = f"""
    <div style="padding-top:10px;">
        <div class="page-title" id="complete_panel">Complete Medication Panel</div>
        <div class="page-subtitle" style="margin-bottom:14px;">
            Your DNA results were evaluated against <strong>{total_drugs}</strong>
            medications across <strong>{total_cats}</strong> drug categories.
        </div>
        {legend}
        <table style="width:100%; border-collapse:collapse; border:1.5px solid #e2eaf4;
            font-family: 'DM Sans', sans-serif;">
            <thead>
                <tr style="background:#4DB7D0;">
                    <th style="padding:8px 12px; text-align:left; font-size:12px;
                        color:white; font-weight:700; width:38%;">Medication</th>
                    <th style="padding:8px 12px; text-align:left; font-size:12px;
                        color:white; font-weight:700; width:25%;">Genes Tested</th>
                    <th style="padding:8px 12px; text-align:left; font-size:12px;
                        color:white; font-weight:700; width:37%;">Evaluation Result</th>
                </tr>
            </thead>
            <tbody>
                {rows_html}
            </tbody>
        </table>
    </div>
    """
    return [_wrap_page(content, name, pg, page_id="complete_panel")], pg + 1


# ============================================================================
# GENOTYPE SUMMARY
# ============================================================================

def genotype_summary_template(df, name, pg, qc_df=None):
    for col in ['Gene', 'Diplotype', 'Phenotype']:
        if col not in df.columns:
            df[col] = 'N/A'
    if 'Activity Score' not in df.columns:
        df['Activity Score'] = ''
    if 'PharmCAT Phenotype' not in df.columns:
        df['PharmCAT Phenotype'] = ''
            
    gene_groups = {}
    gene_activity = {}  # gene → activity score (for display alongside phenotype)
    if not df.empty:
        df = df.sort_values('Gene', na_position='last')
        for row in df.to_dict("records"):
            g = str(row.get('Gene', '')).strip()
            d = str(row.get('Diplotype', '')).strip()
            # Prefer PharmCAT Phenotype (original from report) over the GSI-merge key
            pp = str(row.get('PharmCAT Phenotype', '')).strip()
            p_raw = str(row.get('Phenotype', '')).strip()
            p = pp if (pp and pp.lower() not in ('', 'nan', 'n/a', 'none', 'indeterminate')) else p_raw
            # If PharmCAT Phenotype is blank/indeterminate, fall back to Phenotype column
            if not p or p.lower() in ('nan', 'none'):
                p = p_raw
            a = str(row.get('Activity Score', '')).strip()
            
            if not g or g.lower() in ['nan', 'n/a', 'none']: continue
                
            if ":" in p and p.lower().startswith(g.lower()): 
                p = p.split(":", 1)[1].strip()
                
            if p.lower() in ["no result", "not called", "unknown/unknown", "unknown", "unassigned", "n/a", "", "nan", "uncategorized"]:
                p = "No data available"
            
            # Store activity score for this gene (first non-empty wins)
            if a and a.lower() not in ('', 'nan', 'none', 'n/a', 'no result'):
                if g not in gene_activity:
                    gene_activity[g] = a
            
            if g not in gene_groups:
                gene_groups[g] = {}
            if p not in gene_groups[g]:
                gene_groups[g][p] = []
                
            if d and d.lower() not in ['nan', 'n/a', 'none'] and d not in gene_groups[g][p]:
                gene_groups[g][p].append(d)

    valid_rows = []
    for g in sorted(gene_groups.keys()):
        phenos = gene_groups[g]

        if len(phenos) > 1 and "No data available" in phenos:
            del phenos["No data available"]

        # Collapse multiple phenotype/diplotype combos into a SINGLE row per
        # gene.  Without this, CYP4F2 unphased produces 6 rows, NAT2 ~70 rows,
        # SLCO1B1 5 rows — each a separate line.  One row per gene with all
        # diplotypes and phenotypes joined is dramatically more readable.
        all_dips, all_phens, seen_d = [], [], set()
        for p, d_list in phenos.items():
            if p not in all_phens:
                all_phens.append(p)
            for d in d_list:
                d = str(d).strip()
                if d and d not in seen_d:
                    seen_d.add(d)
                    all_dips.append(d)
        is_collapsed = len(phenos) > 1 or len(all_dips) > 1

        # Apply VKORC1 / IFNL3 pretty-label translation
        translated_dips, translated_phens = [], []
        for d in (all_dips or ["-"]):
            pd_, pp_ = pretty_genotype_pair(g, d)
            if pd_ != d:
                translated_dips.append(pd_)
                if pp_ not in translated_phens:
                    translated_phens.append(pp_)
            else:
                translated_dips.append(d)
        if translated_phens:
            # Genotype-only gene — replace phenotype list with translated
            phen_str = " / ".join(translated_phens)
        else:
            phen_str = " / ".join(all_phens) if all_phens else "-"
        # Append Activity Score if available (matching PharmCAT report display)
        _as = gene_activity.get(g, "")
        if _as and phen_str != "-" and "No data" not in phen_str:
            phen_str += f' <span style="font-size:12px; color:#6b7280;">(Activity Score: {_as})</span>'
        d_str = ", ".join(translated_dips) if translated_dips else "-"

        gene_label = g + (
            ' <span style="font-size:12px; color:#c2410c; font-style:italic; '
            'font-weight:500;">(unphased)</span>'
            if is_collapsed else ''
        )
        valid_rows.append((gene_label, d_str, phen_str))

    if not valid_rows:
        valid_rows = [("No data available", "-", "-")]

    def style_phenotype(pheno):
        # Format metabolizer phenotypes: "normal metabolizer" → "NormalMetabolizer"
        pheno = str(pheno)
        pheno = re.sub(r'\bnormal\s*metabolizer\b', 'Normal Metabolizer', pheno, flags=re.IGNORECASE)
        pheno = re.sub(r'\bpoor\s*metabolizer\b', 'Poor Metabolizer', pheno, flags=re.IGNORECASE)
        pheno = re.sub(r'\bintermediate\s*metabolizer\b', 'Intermediate Metabolizer', pheno, flags=re.IGNORECASE)
        pheno = re.sub(r'\bultrarapid\s*metabolizer\b', 'Ultrarapid Metabolizer', pheno, flags=re.IGNORECASE)
        pheno = re.sub(r'\brapid\s*metabolizer\b', 'Rapid Metabolizer', pheno, flags=re.IGNORECASE)
        pheno = re.sub(r'\bextensive\s*metabolizer\b', 'Extensive Metabolizer', pheno, flags=re.IGNORECASE)

        p_lower = pheno.lower()
        if 'normal' in p_lower or 'metabolizer' in p_lower:
            return f'<span style="font-weight: 600; color: #374151;">{pheno}</span>'
        elif 'no data' in p_lower or 'indeterminate' in p_lower or 'unknown' in p_lower or 'unassigned' in p_lower:
            return f'<span style="background: #f3f4f6; color: #6b7280; padding: 3px 8px; border-radius: 6px; font-weight: 600; display: inline-block;">{pheno}</span>'
        else:
            return f'<span style="font-weight: 600; color: #374151;">{pheno}</span>'

    note = '<div class="info-box" style="margin-top:16px; border: 1px solid #d0e4f7; border-left: 4px solid #1a73e8; background: #f4f9ff; padding: 14px 18px; border-radius: 8px;"><div style="font-weight:700; margin-bottom:8px; font-size: 11.5px; color:#1a1a1a;">Note</div><div style="font-size:11.5px; line-height: 1.65; color: #333;">Genotypes are determined based on the variants in your submitted DNA file, matched against PharmCAT allele definitions. Always consult a qualified healthcare provider before making any medication decisions.</div></div>'

    # ── QC callout ────────────────────────────────────────────────────────────
    # Aggressive dedup: the QC sheet has one row per (gene, diplotype) so a
    # gene with 70 unphased combos (NAT2) produces 70 identical "Uncalled
    # variants: *41; *43; *45..." rows.  Group by gene, take the first
    # non-empty message.  Also extract human-readable text from any JSON
    # blob the GSI sometimes carries in Messages (matching PharmCAT report
    # exception_type = note / footnote entries).  Drop pure "*1 reference
    # allele characterized by absence of variants" footnotes — they appear
    # for every gene that matched a reference allele and add no signal.
    import json as _json
    per_gene = {}  # gene -> first_meaningful_message
    if qc_df is not None and not qc_df.empty:
        _BAD_QC = {"", "nan", "none", "no data"}

        def _clean_qc_message(raw: str) -> str:
            raw = (raw or "").strip()
            if not raw:
                return ""
            # If it looks like JSON, try to extract the 'message' field
            if raw.startswith("{"):
                try:
                    obj = _json.loads(raw)
                    msg = (obj.get("message") or "").strip()
                    if msg:
                        return msg
                except Exception:
                    pass
            return raw

        _BORING_PATTERNS = (
            "characterized by the absence of variants",
            "is on the negative",
            "gene is on the ne",
        )

        for qr in qc_df.to_dict("records"):
            gene = str(qr.get("Gene", "")).strip()
            if not gene or gene.lower() in _BAD_QC or gene in per_gene:
                continue
            msg = _clean_qc_message(str(qr.get("Messages", "")))
            uncalled = str(qr.get("Uncalled", "")).strip()
            if msg.lower() in _BAD_QC:
                msg = ""
            if any(p in msg.lower() for p in _BORING_PATTERNS):
                msg = ""  # drop reference-allele footnotes — pure noise
            detail_parts = []
            if msg:
                detail_parts.append(msg[:300])
            if uncalled and uncalled.lower() not in _BAD_QC:
                detail_parts.append(f"Uncalled variants: {uncalled[:120]}")
            if detail_parts:
                per_gene[gene] = " · ".join(detail_parts)


    # ========================================================================
    # CONTINUOUS GENOTYPE SUMMARY - PAGINATION SAFE
    #
    # IMPORTANT: do not use one large native HTML <table> here. Paged.js can
    # lose the first fragment of a <tr> at a physical page boundary.  That is
    # exactly what caused the SLCO1B1 row to lose its gene label and opening
    # diplotypes in a real report.
    #
    # Instead, every genotype entry is an independent block/grid row.  Rows
    # that fit on a page are kept intact and move to the next physical page when
    # needed.  The surrounding topic itself remains completely fluid and can
    # span any number of pages for any provider/raw-DNA file.
    # ========================================================================

    genotype_rows_html = ""
    for g, d, p in valid_rows:
        styled_p = style_phenotype(p)
        genotype_rows_html += f"""
        <div class="genotype-grid-row">
            <div class="genotype-grid-cell pdf-qa-genotype-gene"
                 style="font-weight:600; color:#111827;">
                {g}
            </div>
            <div class="genotype-grid-cell pdf-qa-genotype-diplotype"
                 style="color:#374151;">
                {d}
            </div>
            <div class="genotype-grid-cell pdf-qa-genotype-phenotype">
                {styled_p}
            </div>
        </div>
        """

    genotype_grid = f"""
    <div class="genotype-grid" style="width:100%; margin:0;">
        <div class="genotype-grid-header">
            <div class="genotype-grid-cell" style="font-weight:700;">Gene</div>
            <div class="genotype-grid-cell" style="font-weight:700;">Diplotype(s)</div>
            <div class="genotype-grid-cell" style="font-weight:700;">Phenotype</div>
        </div>
        {genotype_rows_html}
    </div>
    """

    analytical_html = ""
    if per_gene:
        qc_rows_html = ""
        for gene in sorted(per_gene):
            msg = per_gene[gene]
            qc_rows_html += f"""
            <div class="analytical-grid-row">
                <div class="analytical-grid-cell"
                     style="font-weight:700; font-size:12px; color:#111827;">
                    {gene}
                </div>
                <div class="analytical-grid-cell"
                     style="font-size:12px; color:#374151; line-height:1.4;">
                    {msg}
                </div>
            </div>
            """

        analytical_html = f"""
        <div
            id="genotype_analytical_notes"
            class="analytical-notes-box"
            style="
                margin-top:12px;
                background:#fffbeb;
                border:1px solid #f59e0b;
                border-left:4px solid #f59e0b;
                border-radius:6px;
                padding:10px 12px 8px 12px;
            "
        >
            <div
                class="analytical-notes-heading"
                style="font-size:12px; font-weight:700; color:#92400e; margin-bottom:5px;"
            >
                Analytical Notes
            </div>

            <div style="font-size:11.5px; color:#374151; line-height:1.5; margin-bottom:8px;">
                The genes below produced specific analytical notes during diplotype
                calling (uncalled variants, ambiguous matches, etc.). These do not
                change the recommendations above but may be relevant for clinical
                follow-up.
            </div>

            <div class="analytical-grid" style="width:100%;">
                <div class="analytical-grid-header">
                    <div class="analytical-grid-cell" style="font-weight:700; font-size:12px;">Gene</div>
                    <div class="analytical-grid-cell" style="font-weight:700; font-size:12px;">Note</div>
                </div>
                {qc_rows_html}
            </div>
        </div>
        """

    content = f"""
    <div class="flow-section genotype-summary-section" style="padding-top:10px;">
        <div class="page-title" style="margin-bottom:12px;">
            Genotype Summary
        </div>

        {genotype_grid}
        {analytical_html}
        {note}
    </div>
    """

    return (
        [
            _wrap_page(
                content,
                name,
                pg,
                page_id="genotype_summary"
            )
        ],
        pg + 1
    )


# ============================================================================
# PHARMCAT NO-GUIDANCE DRUG LIST
# ============================================================================

def pharmcat_no_guidance_template(coverage_df, patient_name, pg):
    import html as _html

    _NO_REAL_CATEGORY = {"medications evaluated — no current guidance", "no guideline available", "medications evaluated — has guidance", "nan", ""}

    valid_rows = []
    if coverage_df is not None and not coverage_df.empty:
        for row in coverage_df.sort_values("Drug Name").to_dict("records"):
            drug     = str(row.get("Drug Name", "")).strip()
            category = str(row.get("Drug Category", "")).strip()
            if not drug or drug.lower() in ("nan", ""):
                continue
            display_cat = category if category.lower() not in _NO_REAL_CATEGORY else "—"
            valid_rows.append((drug, display_cat))

    pages = []
    curr_pg = pg
    total_rows = len(valid_rows)
    
    if total_rows == 0:
        content = f"""
        <div
            class="flow-section"
            style="padding-top:10px;"
        >
            <div class="page-title" id="no_guideline_available" style="margin-bottom: 10px;">
                No Guideline Available
            </div>
            <div style="font-size:11.5px; color:#374151; line-height:1.6; margin-bottom:12px;">
                The medications listed below are associated with genes analysed in this report. However, current clinical pharmacogenomics
                guidelines do not yet include dosing recommendations for these specific combinations.
            </div>
            <div class="no-guideline-grid" style="width:100%;">
                <div class="no-guideline-grid-header">
                    <div class="no-guideline-grid-cell" style="font-size:12px; font-weight:700;">Medication</div>
                    <div class="no-guideline-grid-cell" style="font-size:12px; font-weight:700;">Drug Class / Category</div>
                </div>
                <div style="padding:10px; text-align:center; font-size:12px; color:#6b7280;">
                    No additional drugs identified.
                </div>
            </div>
        </div>
        """
        pages.append(_wrap_page(content, patient_name, curr_pg, page_id="no_guideline_available"))
        return pages, pg + 1


    rows_html = ""
    for i, (drug, display_cat) in enumerate(valid_rows):
        bg = "#f9fafb" if i % 2 == 0 else "#ffffff"
        rows_html += (
            f'<div class="keep-together no-guideline-grid-row" style="background:{bg};">'
            f'<div class="no-guideline-grid-cell pdf-qa-no-guideline-medication" '
            f'style="font-size:12px; color:#1e3a5f;">'
            f'<strong>{_html.escape(drug.title())}</strong></div>'
            f'<div class="no-guideline-grid-cell" style="font-size:12px; color:#374151;">'
            f'{_html.escape(display_cat)}</div>'
            f'</div>'
        )
        
    title_html = f"""
    <div class="page-title" id="no_guideline_available" style="margin-bottom: 10px;">
        No Guideline Available
    </div>
    <div style="font-size:11.5px; color:#374151; line-height:1.6; margin-bottom:12px;">
        The medications listed below are associated with genes analysed in this report. However, current clinical pharmacogenomics
        guidelines do not yet include dosing recommendations for these specific combinations.
    </div>
    """
    page_id = "no_guideline_available"
        
    content = f"""
    <div
        class="flow-section"
        style="padding-top:10px;"
    >
        {title_html}
        <div class="no-guideline-grid" style="width:100%;">
            <div class="no-guideline-grid-header">
                <div class="no-guideline-grid-cell" style="font-size:12px; font-weight:700;">Medication</div>
                <div class="no-guideline-grid-cell" style="font-size:12px; font-weight:700;">Drug Class / Category</div>
            </div>
            {rows_html}
        </div>
    </div>
    """
    pages.append(_wrap_page(content, patient_name, curr_pg, page_id=page_id))
    return pages, pg + 1

# ============================================================================
# GENES REQUIRING SPECIALIZED TESTING
# ============================================================================

# Gene-specific explanations — keyed by uppercase gene symbol.
# Each entry has a consumer-facing plain-language explanation and a
# clinician-facing technical reason. Used as primary text in the report;
# any non-generic data-driven reason from the pipeline is appended to the
# clinician note as "Analysis finding: …".
_GENE_REASONS = {
    "HLA-A": {
        "consumer": (
            "Your HLA-A gene type requires a specialized high-resolution test. "
            "HLA genes direct how your immune system reacts to certain medications. "
            "Standard DNA analysis cannot provide the level of detail needed — "
            "a dedicated HLA typing assay is required for an accurate result."
        ),
        "clinician": (
            "HLA allele-level resolution (e.g., HLA-A*31:01) requires dedicated "
            "high-resolution HLA typing. Standard SNP microarray or short-read WGS "
            "cannot reliably phase HLA haplotypes to the 4-digit allele level "
            "required for CPIC hypersensitivity risk assessment (carbamazepine, "
            "oxcarbazepine)."
        ),
    },
    "HLA-B": {
        "consumer": (
            "Your HLA-B gene type requires a specialized high-resolution test. "
            "Certain HLA-B variants are associated with serious medication reactions "
            "(such as severe skin reactions to carbamazepine or abacavir). "
            "A dedicated HLA typing assay is needed for a conclusive result."
        ),
        "clinician": (
            "HLA allele-level resolution (e.g., HLA-B*57:01, HLA-B*15:02, "
            "HLA-B*58:01) requires dedicated high-resolution HLA typing. Standard "
            "SNP microarray or short-read WGS cannot reliably phase HLA haplotypes "
            "to 4-digit allele resolution required for CPIC SJS/TEN and "
            "hypersensitivity reaction risk assessment."
        ),
    },
    "CYP2D6": {
        "consumer": (
            "CYP2D6 is a highly variable gene whose structure can include extra "
            "copies, deletions, or rearrangements that standard DNA tests cannot "
            "reliably detect. Because these structural changes directly determine "
            "how your body processes many common medications, a specialized test "
            "is required for a conclusive result."
        ),
        "clinician": (
            "CYP2D6 structural variation — copy number variation, gene deletion/"
            "duplication, and hybrid alleles (e.g., CYP2D6*36+*10, *68+*4) — "
            "cannot be reliably resolved from standard short-read WGS without "
            "dedicated CNV analysis. Long-read sequencing or a dedicated "
            "CYP2D6 genotyping panel is recommended for definitive diplotype "
            "and activity score assignment."
        ),
    },
    "SLCO1B1": {
        "consumer": (
            "Multiple variants were detected in your SLCO1B1 gene, but standard "
            "testing could not determine how they are arranged across your "
            "chromosomes — a concept called 'phasing'. This arrangement is "
            "critical: it determines whether your risk of statin-related muscle "
            "side effects is elevated. A specialized test can resolve this."
        ),
        "clinician": (
            "SLCO1B1 diplotype assignment requires phasing to distinguish compound "
            "heterozygous from heterozygous configurations (e.g., *1a/*5 vs. "
            "*15/*1a). Without long-read sequencing or family-based phasing, a "
            "definitive functional diplotype cannot be assigned and CPIC-guided "
            "statin myopathy risk stratification cannot be completed."
        ),
    },
    "CYP2C19": {
        "consumer": (
            "Your CYP2C19 result could not be fully determined from this DNA "
            "analysis. The specific combination of variants detected, or the "
            "format of the DNA data provided, requires additional testing to "
            "confirm how the variants are arranged — information needed to "
            "accurately classify your metabolizer type."
        ),
        "clinician": (
            "CYP2C19 diplotype could not be resolved due to ambiguous phasing "
            "or a novel/rare variant preventing definitive star-allele assignment. "
            "Consider long-read sequencing or an alternative clinical genotyping "
            "assay for actionable phenotype classification relevant to clopidogrel, "
            "PPI, and SSRI dosing guidance."
        ),
    },
    "NAT2": {
        "consumer": (
            "The NAT2 gene variant combination detected in your DNA does not "
            "match a standard well-characterized pattern. This may be due to "
            "an unusual or rare combination of variants. A specialist review "
            "or alternative test is needed to accurately determine your "
            "NAT2 metabolizer status."
        ),
        "clinician": (
            "NAT2 haplotype assignment is ambiguous: the detected variant "
            "combination does not correspond to a recognized diplotype in CPIC/"
            "PharmVar nomenclature. Manual clinical genetics review or alternative "
            "genotyping is recommended before issuing phenotype-based dosing "
            "guidance for isoniazid, hydralazine, or other NAT2-metabolized drugs."
        ),
    },
    "CYP4F2": {
        "consumer": (
            "The CYP4F2 gene region could not be reliably analyzed from your "
            "sample, likely due to data quality limitations in that region. "
            "CYP4F2 contributes to warfarin dose calculations — your clinician "
            "should note this gap if anticoagulant therapy is under review."
        ),
        "clinician": (
            "CYP4F2 could not be genotyped from this sample due to insufficient "
            "sequencing coverage or variant quality metrics below analytical "
            "thresholds at the CYP4F2 locus. Warfarin dosing algorithms requiring "
            "CYP4F2*3 (rs2108622) input should proceed without this variant's "
            "contribution, or targeted re-genotyping should be performed."
        ),
    },
    "IFNL3": {
        "consumer": (
            "The IFNL3 (also known as IL28B) gene is used to predict response "
            "to certain hepatitis C treatments. A definitive result could not "
            "be obtained from your DNA analysis. Your clinician can order a "
            "targeted test if hepatitis C treatment is being considered."
        ),
        "clinician": (
            "IFNL3/IL28B (rs12979860, rs8099917) genotyping was indeterminate "
            "from this sample. Clinical correlation with hepatitis C treatment "
            "history is required; targeted re-testing is recommended if "
            "interferon-based or direct-acting antiviral therapy selection "
            "depends on IFNL3 status."
        ),
    },
    "IL28B": {
        "consumer": (
            "The IL28B (IFNL3) gene is used to predict response to certain "
            "hepatitis C treatments. A definitive result could not be obtained "
            "from your DNA analysis. Your clinician can order a targeted test "
            "if hepatitis C treatment is being considered."
        ),
        "clinician": (
            "IFNL3/IL28B (rs12979860, rs8099917) genotyping was indeterminate "
            "from this sample. Clinical correlation with hepatitis C treatment "
            "context is required; targeted re-testing is recommended if "
            "treatment decisions depend on this gene."
        ),
    },
    "G6PD": {
        "consumer": (
            "G6PD is an X-linked gene that affects how red blood cells respond "
            "to certain medications (such as primaquine or dapsone) and foods "
            "like fava beans. Accurately interpreting G6PD results requires "
            "information that could not be fully resolved from this DNA analysis. "
            "A G6PD enzyme activity blood test is often more reliable and "
            "is recommended."
        ),
        "clinician": (
            "G6PD interpretation requires zygosity assessment — hemizygous (XY) "
            "vs. heterozygous (XX) classification significantly affects phenotype "
            "assignment. WHO spectrophotometric enzyme activity measurement is "
            "more clinically reliable than genotyping alone, especially for "
            "heterozygous females, and should be the primary method for clinical "
            "decision-making."
        ),
    },
    "MT-RNR1": {
        "consumer": (
            "This gene is located in your mitochondrial DNA — a small, separate "
            "set of genetic instructions inherited from your mother. Analyzing it "
            "for medication-relevant variants (particularly aminoglycoside antibiotic "
            "sensitivity) requires specialized mtDNA sequencing not included in "
            "standard DNA analysis."
        ),
        "clinician": (
            "MT-RNR1 is mitochondrial; heteroplasmy detection for m.1555A>G and "
            "related variants requires specialized mtDNA sequencing with sufficient "
            "allele-fraction resolution. Standard nuclear DNA pipelines do not "
            "reliably quantify mitochondrial variant heteroplasmy for aminoglycoside "
            "ototoxicity risk assessment."
        ),
    },
}

_GENE_REASONS_DEFAULT = {
    "consumer": (
        "A reliable result for this gene could not be determined from your DNA "
        "analysis. This may be due to technical limitations, a complex gene "
        "structure, or data quality in this region. Your clinician can advise "
        "on whether additional specialized testing is available or needed."
    ),
    "clinician": (
        "Genotyping for this gene was indeterminate from this sample. Review raw "
        "variant data and sequencing quality metrics; an alternative targeted assay "
        "or clinical genetics consultation may be required for actionable "
        "phenotype assignment."
    ),
}

# Reason strings from the pipeline that are too generic to add value
_GENERIC_REASONS = {
    "gene could not be called in this pipeline",
    "gene could not be called",
    "cannot be called",
    "not available",
    "n/a",
    "none",
    "",
}


def specialized_genes_template(df, name, pg):
    def _safe(val):
        s = str(val).strip()
        return s if s.lower() not in ('', 'nan', 'none', 'n/a') else ''

    gene_data = {}
    for row in df.to_dict("records"):
        gene = _safe(row.get("Gene", ""))
        if not gene: continue
        raw_reason = _safe(row.get("Cannot Analyze Reason", ""))

        if gene not in gene_data:
            gene_data[gene] = {"raw_reason": raw_reason, "categories": {}}

        cat  = _safe(row.get("Drug Category", ""))
        drug = _safe(row.get("Drug Name", ""))

        if not drug or not cat or cat == "No entries in current database": continue

        if cat not in gene_data[gene]["categories"]:
            gene_data[gene]["categories"][cat] = []
        if drug not in gene_data[gene]["categories"][cat]:
            gene_data[gene]["categories"][cat].append(drug)

    if not gene_data:
        return [], pg

    gene_list = sorted(gene_data.keys())
    total_genes = len(gene_list)

    pages = []
    curr_pg = pg

    note_html = """
    <div class="keep-together" style="background:#eff6ff; border:1px solid #bfdbfe; border-radius:6px; padding:10px 14px; font-size:11.5px; color:#1e40af; margin-top:8px;">
        <strong>Clinical Note:</strong> The genes listed above could not be reliably genotyped from this sample. Drug-specific recommendations for these medications therefore cannot be personalized based on your DNA. Consult your clinician or a pharmacogenomics specialist for guidance — additional specialized testing (e.g., high-resolution HLA assay, long-read sequencing for CYP2D6, enzyme activity testing for G6PD) may be required.
    </div>"""


    cards_html = ""
    for gene in gene_list:
        raw_reason = gene_data[gene]["raw_reason"]
        cats       = gene_data[gene]["categories"]

        reasons        = _GENE_REASONS.get(gene.upper(), _GENE_REASONS_DEFAULT)
        consumer_text  = reasons["consumer"]
        clinician_text = reasons["clinician"]

        if raw_reason and raw_reason.lower().strip() not in _GENERIC_REASONS:
            clinician_text += (
                f"<br><span style='color:#374151; font-weight:600;'>"
                f"Analysis finding:</span> {raw_reason}"
            )

        if cats:
            cat_groups_html = ""
            for cat in sorted(cats.keys()):
                drugs_in_cat = sorted(cats[cat])

                # Keep the visual format as normal comma-separated prose,
                # but do not give Paged.js one huge inline box to fragment.
                # Short chunks are layout-only safety units; there are no bullets
                # and no provider/gene/drug-specific page rules.
                chunk_lists = _comma_safe_chunks(
                    drugs_in_cat,
                    max_chars=220,
                )

                chunk_html = []
                flat_items = [item for chunk in chunk_lists for item in chunk]
                flat_index = 0

                for chunk in chunk_lists:
                    pieces = []
                    for d in chunk:
                        flat_index += 1
                        suffix = ", " if flat_index < len(flat_items) else ""
                        pieces.append(
                            f'<span class="pdf-qa-specialized-medication specialized-med-item">'
                            f'{html_lib.escape(str(d).title())}'
                            f'</span>{suffix}'
                        )
                    chunk_html.append(
                        '<div class="specialized-med-chunk">'
                        + ''.join(pieces)
                        + '</div>'
                    )

                first_chunk = chunk_html[0] if chunk_html else ''
                remaining_chunks = ''.join(chunk_html[1:])

                cat_groups_html += f"""
                <div class="specialized-med-group">
                    <div class="specialized-med-first-fragment">
                        <div class="specialized-med-category">
                            {html_lib.escape(str(cat))}
                        </div>
                        <div class="specialized-med-list">
                            {first_chunk}
                        </div>
                    </div>
                    <div class="specialized-med-list">
                        {remaining_chunks}
                    </div>
                </div>
                """

            drug_section_html = f"""
<div
    class="specialized-affected-section"
    style="
        padding:6px 12px 8px 12px;
        border-top:1px solid #dde2ea;
    "
>

    <div
        class="section-block-header keep-with-next"
        style="
            font-size:12px;
            font-weight:700;
            color:#374151;
            text-transform:uppercase;
            letter-spacing:0.5px;
            margin-bottom:6px;
        "
    >
        Affected Medications
    </div>

    {cat_groups_html}

</div>
"""
        else:
            drug_section_html = """
            <div style="padding:8px 12px 10px 12px; border-top:1px solid #dde2ea; font-size:11.5px; color:#6b7280; font-style:italic;">
                Drug interaction data for this gene is not yet available in our database. Consult external clinical resources for up-to-date recommendations.
            </div>"""

        cards_html += f"""
<div
    class="specialized-gene-card"
    style="
        margin-bottom:14px;
        border:1px solid #dde2ea;
    "
>

    <div class="specialized-gene-top">

        <div style="
            background:#f3f5f8;
            padding:7px 12px;
            font-size:11.5px;
            font-weight:700;
            color:#1a1a1a;
            border-bottom:1px solid #dde2ea;
            display:flex;
            justify-content:space-between;
            align-items:center;
        ">
            <span style="
                font-size:11.5px;
                font-weight:800;
                color:#0D3B7A;
            ">
                {gene}
            </span>

            <span style="
                background:#fef3c7;
                color:#92400e;
                border:1px solid #fcd34d;
                border-radius:4px;
                padding:2px 8px;
                font-size:12px;
                font-weight:700;
            ">
                Cannot Be Called
            </span>
        </div>

        <div style="
            padding:8px 12px 6px 12px;
            font-size:11.5px;
            color:#374151;
            line-height:1.55;
        ">
            {consumer_text}
        </div>

        <div style="
            padding:6px 12px 8px 12px;
            background:#f8fafc;
            font-size:12px;
            color:#6b7280;
            line-height:1.5;
            border-top:1px solid #f0f2f5;
            border-bottom:1px solid #dde2ea;
        ">
            <span style="
                font-size:12px;
                font-weight:700;
                color:#9ca3af;
                text-transform:uppercase;
                letter-spacing:0.4px;
            ">
                For Clinicians &mdash;
            </span>

            {clinician_text}
        </div>

    </div>

    {drug_section_html}

</div>
"""

    title_html = f"""
    <div class="page-title" id="specialized_genes" style="margin-bottom: 18px;">Genes Requiring Specialized Testing</div>
    """
    page_id = "specialized_genes"

    content = f"""
    <div
        class="flow-section"
        style="padding-top:10px;"
    >
        {title_html}
        {cards_html}
        {note_html}
    </div>
    """
    pages.append(_wrap_page(content, name, curr_pg, page_id=page_id))
    return pages, pg + 1