#!/usr/bin/env python3
import sys
import pandas as pd
from cyvcf2 import VCF

if len(sys.argv) != 6:
    sys.stderr.write(
        "Usage: generate_lean_report.py <vcf.vcf.gz> <exon_cov.txt> <r1r2.tsv> <frstrand.tsv> <out.xlsx>\n"
    )
    sys.exit(1)

vcf_in, coverage_summary, r1r2_summary, fr_summary, xlsx_out = sys.argv[1:]

# --- Load exon coverage summary ---
coverage_data = {}
with open(coverage_summary) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue
        region = parts[0]  # e.g. 1:17018000-17019000
        cov = {}
        for p in parts[1:]:
            if ":" in p:
                key, val = p.split(":", 1)
                cov[key.strip()] = val.strip().replace("%", "")
        coverage_data[region] = cov

def get_exon_cov(chrom, pos):
    # pos is 1-based; regions are half-open? We used inclusive start/end earlier.
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
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 6:  # chrom start end r1 r2 ratio
            continue
        chrom, start, end, r1, r2, ratio = parts[:6]
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
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 6:  # chrom start end fwd rev frac balance
            continue
        chrom, start, end, fwd, rev, ffrac = parts[:6]
        # Some versions include a 7th "balance" col; keep robust:
        rfrac = parts[6] if len(parts) > 6 else None
        fr_data[(chrom, int(start), int(end))] = (fwd, rev, ffrac, rfrac)

def get_fr(chrom, pos):
    for (c, s, e), vals in fr_data.items():
        if str(chrom) == c and s <= pos <= e:
            return vals
    return ("NA","NA","NA","NA")

# --- Parse VCF ---
vcf = VCF(vcf_in)

# Extract CSQ format (VEP)
csq_format = []
for h in vcf.raw_header.split("\n"):
    if h.startswith("##INFO=<ID=CSQ"):
        try:
            csq_format = h.split("Format: ")[1].rstrip('">').split("|")
        except Exception:
            pass
        break

def first_csq_dict(variant):
    csq_entries = variant.INFO.get("CSQ")
    if not csq_entries or not csq_format:
        return {}
    first = csq_entries.split(",")[0].split("|")
    # pad to length
    if len(first) < len(csq_format):
        first = first + [""] * (len(csq_format) - len(first))
    return dict(zip(csq_format, first))

records = []
for v in vcf:
    chrom, pos, ref = v.CHROM, v.POS, v.REF
    alts = v.ALT
    alt = ",".join(alts) if alts else ""

    # depth & VAF (from FORMAT, else INFO)
    dp = None
    if "DP" in v.FORMAT:
        try:
            dp = int(v.format("DP")[0][0])
        except Exception:
            dp = None
    if dp is None:
        dp = v.INFO.get("DP")

    vaf = None
    if "VAF" in v.FORMAT:
        try:
            vaf = float(v.format("VAF")[0][0])
        except Exception:
            vaf = None

    # zygosity
    gt = v.genotypes[0][:2] if v.genotypes else (None, None)
    if gt[0] is None or gt[1] is None:
        zyg = "NA"
    else:
        zyg = "het" if gt[0] != gt[1] else ("hom" if gt[0] == 1 else "ref")

    # AD (ref, alt) if available
    ad_ref = ad_alt = None
    if "AD" in v.FORMAT:
        try:
            ad_ref, ad_alt = v.format("AD")[0][:2]
        except Exception:
            pass

    ann = first_csq_dict(v)
    gene       = ann.get("SYMBOL") or None
    transcript = ann.get("Feature") or None
    hgvsc      = ann.get("HGVSc") or None
    hgvsp      = ann.get("HGVSp") or None
    consequence= ann.get("Consequence") or None
    impact     = ann.get("IMPACT") or None
    gnomad_af  = ann.get("gnomADg_AF") or ann.get("gnomADe_AF") or None
    revel      = ann.get("REVEL") or None
    spliceai   = ann.get("SpliceAI") or None
    clinvar    = ann.get("CLIN_SIG") or None
    # Try to grab ClinVar star rating when present
    clinstar   = ann.get("CLIN_SIG_CONF") or None

    cov = get_exon_cov(chrom, pos)
    r1r2 = get_r1r2(chrom, pos)
    fr   = get_fr(chrom, pos)

    # absolute ratios from counts
    def safe_ratio(a, b):
        try:
            a = float(a); b = float(b)
            return round(a/b, 2) if b > 0 else None
        except Exception:
            return None
    r1r2_abs = safe_ratio(r1r2[0], r1r2[1])
    fr_abs   = safe_ratio(fr[0],   fr[1])

    # combined HGVS
    hgvs_combined = None
    if gene or hgvsc or hgvsp:
        cg = hgvsc or ""
        pg = hgvsp or ""
        hgvs_combined = f"{gene or ''}: {cg or ''} / {pg or ''}".strip()

    # Run incidence within the run (single-sample now)
    run_count = 1

    row = {
        "Chrom": chrom,
        "Start": pos,
        "End": pos + len(ref) - 1 if ref else pos,
        "Location": f"{chrom}:{pos}-{pos + len(ref) - 1 if ref else pos}",
        "Ref": ref,
        "Alt": alt,
        "Gene": gene,
        "Transcript": transcript,
        "Consequence": consequence,
        "Impact": impact,
        "HGVSc": hgvsc,
        "HGVSp": hgvsp,
        "HGVS_combined": hgvs_combined,
        "Depth": dp,
        "VAF": vaf,
        "Zygosity": zyg,
        "ExonCov20": cov.get(">=20x", "NA"),
        "ExonCov30": cov.get(">=30x", "NA"),
        "ExonCov50": cov.get(">=50x", "NA"),
        "ExonCov100": cov.get(">=100x", "NA"),
        "R1_count": r1r2[0],
        "R2_count": r1r2[1],
        "R1R2_ratio": r1r2[2],
        "R1R2_abs_ratio": r1r2_abs,
        "FWD_count": fr[0],
        "REV_count": fr[1],
        "FWD_frac": fr[2],
        "REV_frac": fr[3],
        "FR_abs_ratio": fr_abs,
        "RunCount": run_count,
        "gnomAD_AF": gnomad_af,
        "SpliceAI": spliceai,
        "REVEL": revel,
        "ClinVar": clinvar,
        "ClinVar_Stars": clinstar,
    }
    records.append(row)

df = pd.DataFrame(records)

# Coerce numerics safely
for col in ["VAF", "Depth", "gnomAD_AF", "ExonCov20", "ExonCov30", "ExonCov50", "ExonCov100"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# High-confidence logic
is_het_high = (df["Zygosity"] == "het") & df["VAF"].between(0.35, 0.65, inclusive="both")
is_hom_high = (df["Zygosity"] == "hom") & (df["VAF"] >= 0.90)
high_conf = is_het_high | is_hom_high

# Population rarity
is_rare = (df["gnomAD_AF"].fillna(0) < 0.01)
is_common = (df["gnomAD_AF"].fillna(0) >= 0.01)

# Buckets
df_unique_high_rare   = df[ high_conf & is_rare ].copy()
df_pop_high_common    = df[ high_conf & is_common ].copy()
df_medium_conf        = df[ df["VAF"].between(0.20, 0.34, inclusive="both") | df["VAF"].between(0.65, 0.89, inclusive="both") ].copy()
df_low_conf           = df[ (df["VAF"] < 0.20) ].copy()

# Coverage summary tab
cov_cols = ["ExonCov20", "ExonCov30", "ExonCov50", "ExonCov100"]
cov_summary = df[cov_cols].describe().transpose()

# Gene summary
def uniq_join(x):
    vals = sorted(set([s for s in x.dropna().astype(str) if s]))
    return ", ".join(vals)[:2000]
gene_summary = (
    df.groupby("Gene", dropna=False)
      .agg(Variants=("Location","count"),
           Mean_VAF=("VAF","mean"),
           Min_VAF=("VAF","min"),
           Max_VAF=("VAF","max"),
           Consequences=("Consequence", uniq_join))
      .sort_values("Variants", ascending=False)
)

# ClinVar pathogenic bucket (broad match)
clinvar_path = df[df["ClinVar"].str.contains("pathogenic", case=False, na=False)].copy()

# Order columns nicely for exports
cols_order = [
    "Chrom","Start","End","Location","Ref","Alt",
    "Gene","Transcript","Consequence","Impact",
    "HGVSc","HGVSp","HGVS_combined",
    "Depth","VAF","Zygosity",
    "ExonCov20","ExonCov30","ExonCov50","ExonCov100",
    "R1_count","R2_count","R1R2_ratio","R1R2_abs_ratio",
    "FWD_count","REV_count","FWD_frac","REV_frac","FR_abs_ratio",
    "RunCount","gnomAD_AF","SpliceAI","REVEL","ClinVar","ClinVar_Stars"
]
df = df.reindex(columns=[c for c in cols_order if c in df.columns])

with pd.ExcelWriter(xlsx_out, engine="openpyxl") as writer:
    df.to_excel(writer, index=False, sheet_name="All Variants")
    df_unique_high_rare.reindex(columns=df.columns).to_excel(writer, index=False, sheet_name="Unique HighConf Rare")
    df_pop_high_common.reindex(columns=df.columns).to_excel(writer, index=False, sheet_name="Population HighConf")
    df_medium_conf.reindex(columns=df.columns).to_excel(writer, index=False, sheet_name="Medium Confidence")
    df_low_conf.reindex(columns=df.columns).to_excel(writer, index=False, sheet_name="Low Confidence")
    cov_summary.to_excel(writer, sheet_name="Coverage Summary")
    gene_summary.to_excel(writer, sheet_name="Gene Summary")
    clinvar_path.reindex(columns=df.columns).to_excel(writer, index=False, sheet_name="ClinVar Pathogenic")

print(f"Wrote Excel report → {xlsx_out}")
