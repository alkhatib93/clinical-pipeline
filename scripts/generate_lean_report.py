import sys
import pandas as pd
from cyvcf2 import VCF

vcf_in = sys.argv[1]
coverage_summary = sys.argv[2]
r1r2_summary = sys.argv[3]
fr_summary = sys.argv[4]
tsv_out = sys.argv[5]

# --- Load exon coverage summary ---
coverage_data = {}
with open(coverage_summary) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 5:
            continue
        region = parts[0]  # e.g. 1:17018000-17019000
        cov = {}
        for p in parts[1:]:
            key, val = p.split(":")
            cov[key.strip()] = val.strip().replace("%","")
        coverage_data[region] = cov

def get_exon_cov(chrom, pos):
    for region, cov in coverage_data.items():
        rchrom, coords = region.split(":")
        start, end = map(int, coords.split("-"))
        if str(chrom) == rchrom and start <= pos <= end:
            return cov
    return {">=20x":"NA", ">=30x":"NA", ">=50x":"NA", ">=100x":"NA"}

# --- Load R1/R2 summary ---
r1r2_data = {}
with open(r1r2_summary) as f:
    for line in f:
        chrom, start, end, r1, r2, ratio = line.strip().split("\t")
        r1r2_data[(chrom, int(start), int(end))] = (r1, r2, ratio)

def get_r1r2(chrom, pos):
    for (c, s, e), vals in r1r2_data.items():
        if str(chrom) == c and s <= pos <= e:
            return vals
    return ("NA","NA","NA")

# --- Load Forward/Reverse summary ---
fr_data = {}
with open(fr_summary) as f:
    for line in f:
        chrom, start, end, fwd, rev, ffrac, rfrac = line.strip().split("\t")
        fr_data[(chrom, int(start), int(end))] = (fwd, rev, ffrac, rfrac)

def get_fr(chrom, pos):
    for (c, s, e), vals in fr_data.items():
        if str(chrom) == c and s <= pos <= e:
            return vals
    return ("NA","NA","NA","NA")

# --- Parse VCF ---
vcf = VCF(vcf_in)

# Extract CSQ format
csq_line = [h for h in vcf.raw_header.split("\n") if h.startswith("##INFO=<ID=CSQ")][0]
csq_format = csq_line.split("Format: ")[1].rstrip('">').split("|")

records = []
for variant in vcf:
    chrom, pos, ref, alt = variant.CHROM, variant.POS, variant.REF, ",".join(variant.ALT)

    dp = variant.format("DP")[0][0] if "DP" in variant.FORMAT else variant.INFO.get("DP")
    vaf = variant.format("VAF")[0][0] if "VAF" in variant.FORMAT else None

    gt = variant.genotypes[0][:2]
    zygosity = "het" if gt[0] != gt[1] else ("hom" if gt[0] == 1 else "ref")

    if "AD" in variant.FORMAT:
        ad_ref, ad_alt = variant.format("AD")[0]
    else:
        ad_ref, ad_alt = None, None

    gene = transcript = hgvsc = hgvsp = consequence = clinvar = gnomad_af = revel = spliceai = None
    csq_entries = variant.INFO.get("CSQ")
    if csq_entries:
        first_csq = csq_entries.split(",")[0].split("|")
        ann = dict(zip(csq_format, first_csq))
        gene = ann.get("SYMBOL")
        transcript = ann.get("Feature")
        hgvsc = ann.get("HGVSc")
        hgvsp = ann.get("HGVSp")
        consequence = ann.get("Consequence")
        gnomad_af = ann.get("gnomADg_AF") or ann.get("gnomADe_AF")
        revel = ann.get("REVEL")
        spliceai = ann.get("SpliceAI")
        clinvar = ann.get("CLIN_SIG")
    else:
        ann = {}

    cov = get_exon_cov(chrom, pos)
    r1r2 = get_r1r2(chrom, pos)
    fr = get_fr(chrom, pos)
    try:
        r1 = int(r1r2[0])
        r2 = int(r1r2[1])
        r1r2_abs = round(r1 / r2, 2) if r2 > 0 else None
    except:
        r1r2_abs = None

    try:
        fwd = int(fr[0])
        rev = int(fr[1])
        fr_abs = round(fwd / rev, 2) if rev > 0 else None
    except:
        fr_abs = None

    flag = "PASS"
    try:
        cov20 = float(cov[">=20x"]) if cov and cov[">=20x"] not in [None, ""] else None
    except (ValueError, KeyError):
        cov20 = None
    if cov20 is not None and cov20 < 95.0:
        flag = "LowCov"

    records.append({
        "Location": f"{chrom}:{pos}-{pos+len(ref)-1}",
        "Ref": ref,
        "Alt": alt,
        "Gene": gene,
        "Transcript": transcript,
        "Consequence": consequence,
        "Impact": ann.get("IMPACT"),
        "HGVSc": hgvsc,
        "HGVSp": hgvsp,
        "Depth": dp,
        "VAF": vaf,
        "Zygosity": zygosity,
        "ExonCov20": cov[">=20x"],
        "ExonCov30": cov[">=30x"],
        "ExonCov50": cov[">=50x"],
        "ExonCov100": cov[">=100x"],
        "R1_count": r1r2[0],
        "R2_count": r1r2[1],
        "R1R2_ratio": r1r2[2],
        "R1R2_abs_ratio": r1r2_abs,
        "FWD_count": fr[0],
        "REV_count": fr[1],
        "FWD_frac": fr[2],
        "REV_frac": fr[3],
        "FR_abs_ratio": fr_abs,
        #"QC_Flag": flag,
        "SpliceAI": spliceai,
        "REVEL": revel,
        "gnomAD_AF": gnomad_af,
        "ClinVar": clinvar,
    })

df = pd.DataFrame(records)
with pd.ExcelWriter(tsv_out) as writer:
    # Full table
    df.to_excel(writer, index=False, sheet_name="All Variants")

    # High confidence
    high_conf = df[
        ((df["Zygosity"] == "het") & (df["VAF"].between(0.35, 0.65))) |
        ((df["Zygosity"] == "hom") & (df["VAF"] >= 0.9))
    ]
    high_conf.to_excel(writer, index=False, sheet_name="High Confidence")

    # Medium confidence
    medium_conf = df[
        (df["VAF"].between(0.1, 0.34)) | (df["VAF"].between(0.65, 0.89))
    ]
    medium_conf.to_excel(writer, index=False, sheet_name="Medium Confidence")

    # Low confidence
    low_conf = df[df["VAF"] < 0.2]
    low_conf.to_excel(writer, index=False, sheet_name="Low Confidence")

    # Coverage summary (optional)
    cov_cols = ["ExonCov20", "ExonCov30", "ExonCov50", "ExonCov100"]
    df[cov_cols] = df[cov_cols].apply(pd.to_numeric, errors="coerce")
    cov_summary = df[cov_cols].describe().transpose()
    cov_summary.to_excel(writer, sheet_name="Coverage Summary")

    # Gene summary (optional)
    gene_summary = df.groupby("Gene").agg({
        "Location": "count",
        "VAF": ["mean", "min", "max"],
        "Consequence": lambda x: ','.join(set(x.dropna()))
    }).sort_values(("Location", "count"), ascending=False)
    gene_summary.to_excel(writer, sheet_name="Gene Summary")

    # ClinVar
    clinvar_path = df[df["ClinVar"].str.contains("pathogenic", case=False, na=False)]
    clinvar_path.to_excel(writer, index=False, sheet_name="ClinVar Pathogenic")

#df.to_csv(tsv_out, sep="\t", index=False)
