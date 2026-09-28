# PGx Xcode Life Report Generator

This project creates a personalized pharmacogenomics report from a raw DNA input file. It processes genotype data, runs PharmCAT, merges drug-gene recommendations, and produces a final PDF report.

## How the report pipeline works

### Step 1: Convert raw DNA to VCF
This step reads the uploaded raw genotype file (for example 23andMe or AncestryDNA output) and converts it into a VCF file that PharmCAT can understand.

It extracts relevant genetic variants and keeps only the values needed for pharmacogenomic analysis.

### Step 2: Run PharmCAT
This step runs the PharmCAT Java-based pharmacogenomics engine.

PharmCAT analyzes the patient genotype and generates structured pharmacogenomic recommendations and gene-specific result data.

### Step 2b: Warfarin calculator input
This is a supporting step that creates a Warfarin-relevant CSV file from the VCF and PharmCAT JSON output.

It captures markers such as:

- CYP2C9
- CYP4F2
- VKORC1
- other warfarin-related genotype markers

This is used for Warfarin-specific calculation and reporting support.

### Step 3: JSON to summary Excel
This step converts PharmCAT output into a summary Excel workbook that can be processed downstream.

It organizes patient genotype details and recommendation data into structured forms for the report engine.

### Step 4: Add all sources recommendations
This step merges multiple sources of pharmacogenomic guidance, including PharmCAT, GSI, and additional recommendation sources.

It combines data into the final drug-gene recommendation set used by the report.

### Step 5: Build the final PDF report
This is the last stage.

The project assembles the report HTML, renders it using Paged.js, and then merges it with the front/back cover pages to produce the final PDF.

## Required setup

### Python dependencies

```bash
pip install -r requirements.txt
```

Current Python packages in the project include:

- pandas
- openpyxl
- PyMuPDF
- PyPDF2

### Java requirement

The PharmCAT step requires Java.

Check whether Java is installed:

```bash
java -version
```

If Java is missing, install a JDK (Java 21 is a good option) and make sure it is on your PATH.

### Node / Paged.js requirement

The PDF render step uses the Paged.js CLI.

Install the project’s Node dependency:

```bash
npm install
```

This installs:

- pagedjs-cli@0.4.3

## Run the full pipeline

From the project root, run:

```bash
python main.py 01_inputs/sample_AncestryDNA.txt --name ancestry
```

You can replace the input file and patient name with your own values.

Example:

```bash
python main.py 01_inputs/sample_23andMe_v3.txt --name "John Smith"
```

## Output

Final PDFs are generated in:

```text
results/reports_drugwise_pdf/
```

## Optional debug mode

To show more internal logs during generation, use:

Windows:

```bash
set PGX_VERBOSE=1
```

macOS/Linux:

```bash
export PGX_VERBOSE=1
```

Then rerun the pipeline.
