#!/usr/bin/env python3
import sys, os, re, argparse
import pandas as pd
from cyvcf2 import VCF

# -------------------------
# CLI
# -------------------------
p = argparse.ArgumentParser(description="Build LEAN Excel deliverable")
p.add_argument("vcf", help="VCF (bgzip+tabix) with VEP+ClinVar if available")
p.add_argument("coverage_summary", help="exon_coverage_summary.sorted.txt")
p.add_argument("r1r2_summary", help="r1r2_per_exon.tsv")
p.add_argument("fr_summary", help="frstrand_per_exon.tsv")
p.add_argument("xlsx_out", help="output Excel path, e.g. results/sample.xlsx")

# Optional QC inputs for Sample Summary (all optional)
p.add_argument("--sample-id", default=None, help="Sample ID override (else inferred from VCF header or filename)")
p.add_argument("--assay", default="", help="WES/WGS (free text)")
p.add_argument("--build", default="GRCh38", help="Reference build")
p.add_argument("--flagstat", default=None, help="samtools flagstat output")
p.add_argument("--stats", default=None, help="samtools stats output")
p.add_argument("--picard-align", default=None, help="Picard CollectAlignmentSummaryMetrics.txt")
p.add_argument("--picard-insert", default=None, help="Picard CollectInsertSizeMetrics.txt")
p.add_argument("--mosdepth-summary", default=None, help="mosdepth *.summary.txt (plain, not gz)")
p.add_argument("--verifybamid", default=None, help="verifyBamID2 *.selfSM or .summary output (FREEMIX)")
p.add_argument("--sexcheck", default=None, help="file with inferred sex (one line like: Sex: Male)")
p.add_argument("--sf-genes", default=None, help="ACMG SF gene list (one gene SYMBOL per line)")
p.add_argument("--acmg-thresholds", default=None,
               help="mosdepth thresholds.bed(.gz) run with ACMG BED (cols: chrom start end [region] 20X 30X)")
p.add_argument("--gaps20", default=None, help="annotated gaps <20x BED (chrom start end RegionLabel)")
p.add_argument("--gaps30", default=None, help="annotated gaps <30x BED (chrom start end RegionLabel)")



args = p.parse_args()

# -------------------------
# Helpers: safe parsers
# -------------------------
def safe_float(x):
    try: return float(x)
    except: return None

def sample_from_vcf(vcf_path):
    try:
        v = VCF(vcf_path)
        if v.samples: return v.samples[0]
    except: pass
    # fallback to filename stem
    return os.path.basename(vcf_path).split(".")[0]

def parse_flagstat(path):
    """Return dict with Total_Reads, Mapped_Reads, Mapped_Percent, Duplicate_Reads, Duplicate_Percent"""
    out = {}
    if not path or not os.path.exists(path): return out
    with open(path) as f:
        for line in f:
            s = line.strip()
            # total line is often like: "123456 + 0 in total ..."
            if " in total " in s and "QC-passed" in s:
                n = s.split()[0]
                if n.isdigit(): out["Total_Reads"] = int(n)
            # mapped line: "123456 + 0 mapped (99.97% : N/A)"
            elif re.search(r"\bmapped\s*\(", s):
                n = s.split()[0]
                out["Mapped_Reads"] = int(n) if n.isdigit() else None
                m = re.search(r"\(([\d\.]+)%", s)
                if m: out["Mapped_Percent"] = safe_float(m.group(1))
            # duplicates line (may be absent if not marked)
            elif re.search(r"\bduplicates\b", s) and "(" in s:
                n = s.split()[0]
                out["Duplicate_Reads"] = int(n) if n.isdigit() else None
                m = re.search(r"\(([\d\.]+)%", s)
                if m: out["Duplicate_Percent"] = safe_float(m.group(1))
    return out

def parse_samtools_stats(path):
    """Return Ti/Tv, Het/Hom (approx) from samtools stats if present"""
    out = {}
    if not path or not os.path.exists(path): return out
    ti = tv = het = homalt = mean_insert_size = 0
    with open(path) as f:
        for line in f:
            # SN  number of first fragments:   123
            if line.startswith("SN"):
                if "number of heterozygous sites" in line:
                    het = int(line.strip().split()[-1])
                if "number of homozygous-ALT sites" in line:
                    homalt = int(line.strip().split()[-1])
                if "TSTV" in line:
                    # some builds print: "SN\tts/tv\t2.02"
                    parts = line.strip().split()
                    try: out["TiTv"] = float(parts[-1])
                    except: pass
                if "insert size average" in line:
                    mean_insert_size = float(line.strip().split()[-1])  
            # For bcftools stats you’d parse different keys; kept minimal here.
    if het or homalt:
        try: out["Het_Hom_Ratio"] = round(het / max(homalt,1), 3)
        except: pass
    if mean_insert_size:
        out["Mean_Insert_Size"] = mean_insert_size
    return out

def parse_picard_alignment(path):
    """Picard CollectAlignmentSummaryMetrics: %Q30, %mapped, mean read length, etc."""
    out = {}
    if not path or not os.path.exists(path): return out
    # Picard has a header section; data lines often start with CATEGORY
    # We'll pick the FIRST row for 'PAIR' or 'UNPAIRED' if present
    df = None
    try:
        with open(path) as f:
            lines = [l.strip() for l in f if l.strip() and not l.startswith("#")]
        # Find header line (tab-delimited)
        hdr_idx = next(i for i,l in enumerate(lines) if l.startswith("CATEGORY"))
        header = lines[hdr_idx].split("\t")
        row = lines[hdr_idx+1].split("\t")
        d = dict(zip(header,row))
        # Common keys (presence depends on Picard version)
        out["Pct_Q30"] = safe_float(d.get("PF_READS_ALIGNED_Q20_BASES", ""))  # fallback if Q30 missing
        if "PCT_PF_READS_ALIGNED" in d: out["Pct_Mapped"] = safe_float(d["PCT_PF_READS_ALIGNED"])*100 if d["PCT_PF_READS_ALIGNED"] not in (None,"") else None
        if "MEAN_READ_LENGTH" in d: out["Mean_Read_Length"] = safe_float(d["MEAN_READ_LENGTH"])
        if "MEAN_INSERT_SIZE" in d: out["Mean_Insert_Size"] = safe_float(d["MEAN_INSERT_SIZE"])
        if "PCT_ADAPTER" in d: out["Pct_Adapter"] = safe_float(d["PCT_ADAPTER"])*100 if d["PCT_ADAPTER"] not in (None,"") else None
    except Exception:
        pass
    return out

def parse_picard_insert(path):
    """Picard InsertSizeMetrics: mean insert size if not captured above."""
    out = {}
    if not path or not os.path.exists(path): return out
    try:
        with open(path) as f:
            lines = [l.strip() for l in f if l.strip()]
        hdr_idx = next(i for i,l in enumerate(lines) if l.startswith("MEDIAN_INSERT_SIZE"))
        header = lines[hdr_idx].split("\t")
        row = lines[hdr_idx+1].split("\t")
        d = dict(zip(header,row))
        if "MEAN_INSERT_SIZE" in d:
            out["Mean_Insert_Size"] = safe_float(d["MEAN_INSERT_SIZE"])
    except Exception:
        pass
    return out

def parse_mosdepth_summary(path):
    """Return Mean_Coverage and %>=10x/20x/30x if mosdepth summary available."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    try:
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            return out
        row = df.iloc[0]  # assume one sample per file

        out["Mean_Coverage"] = float(row.get("MeanDepth", 0))

        for col in ["Pct>=10x","Pct>=20x","Pct>=30x","Pct>=50x","Pct>=100x"]:
            if col in df.columns:
                try:
                    out[col] = float(row[col])
                except Exception:
                    pass
    except Exception:
        # silently ignore parse errors
        return out
    return out


def read_acmg_thresholds(th_path, sf_genes=None):
    """
    Read mosdepth thresholds generated WITH ACMG BED (thresholds 20,30).
    If 'region' column is missing, synthesize RegionLabel as 'chrom:start-end'.
    Returns per-exon rows with RegionLabel, Chrom, ExonStart, ExonEnd, ExonLen, Pct>=20x, Pct>=30x.
    """
    import pandas as pd, os
    if not th_path or not os.path.exists(th_path):
        return pd.DataFrame(columns=["RegionLabel","Chrom","ExonStart","ExonEnd","ExonLen","Pct>=20x","Pct>=30x","Gene","Transcript","Exon"])

    df = pd.read_csv(th_path, sep="\t", comment="#", low_memory=False)
    # Normalize column names
    cols = {c.lower(): c for c in df.columns}
    def col(name, default=None): return cols.get(name.lower(), default)

    # If header is missing / different, coerce
    if col("chrom") is None or col("start") is None or col("end") is None:
        # assume fixed layout: chrom start end region 20X 30X ...
        df.columns = ["chrom","start","end","region","20X","30X"] + list(df.columns[6:])
        cols = {c.lower(): c for c in df.columns}

    chrom = col("chrom") or "chrom"
    start = col("start") or "start"
    end   = col("end")   or "end"
    region= col("region")  # may be None if no label in BED

    # Synthesize RegionLabel if missing/blank
    if (region is None) or (region not in df.columns):
        df["RegionLabel"] = df[chrom].astype(str) + ":" + df[start].astype(int).astype(str) + "-" + df[end].astype(int).astype(str)
    else:
        # If region exists but is empty/unknown, still fall back to coords
        df["RegionLabel"] = df[region].astype(str)
        mask_empty = df["RegionLabel"].isin(["", ".", "unknown", "UNKNOWN", "None"])
        df.loc[mask_empty, "RegionLabel"] = df.loc[mask_empty, chrom].astype(str) + ":" + df.loc[mask_empty, start].astype(int).astype(str) + "-" + df.loc[mask_empty, end].astype(int).astype(str)

    # Optional gene parsing (will be '.' if no names)
    parts = df["RegionLabel"].astype(str).str.split("|", n=2, expand=True)
    df["Gene"] = parts[0] if parts.shape[1] > 0 else "."
    df["Transcript"] = parts[1] if parts.shape[1] > 1 else "."
    df["Exon"] = parts[2] if parts.shape[1] > 2 else "."

    if sf_genes:
        df = df[df["Gene"].isin(sf_genes)]

    df["ExonLen"] = df[end] - df[start]
    if "20X" not in df.columns or "30X" not in df.columns:
        # if thresholds file has more cols, pick those with digits '20' / '30'
        for t in ("20","30"):
            cand = [c for c in df.columns if "".join(ch for ch in str(c) if ch.isdigit()) == t]
            if cand:
                df[f"{t}X"] = df[cand[0]]

    df["Pct>=20x"] = (df["20X"] / df["ExonLen"] * 100).round(2)
    df["Pct>=30x"] = (df["30X"] / df["ExonLen"] * 100).round(2)

    out = df.rename(columns={chrom:"Chrom", start:"ExonStart", end:"ExonEnd"})[
        ["RegionLabel","Gene","Transcript","Exon","Chrom","ExonStart","ExonEnd","ExonLen","Pct>=20x","Pct>=30x"]
    ].drop_duplicates(subset=["RegionLabel","Chrom","ExonStart","ExonEnd"])
    return out


def load_gap_bed_for_agg(path, threshold_tag, sf_genes=None):
    """
    Read gaps BED (chrom start end [RegionLabel]).
    If RegionLabel is missing, synthesize it from chrom:start-end of the *exon* is not known —
    but since we intersected with the exon BED when annotating, col4 should already be RegionLabel=exon coords.
    Aggregates per RegionLabel: counts, total bp, and intervals list.
    """
    import pandas as pd, os
    if not path or not os.path.exists(path):
        return pd.DataFrame(columns=["RegionLabel", f"Gaps{threshold_tag}_n", f"Gaps{threshold_tag}_bp", f"Gaps{threshold_tag}_intervals"])

    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype={0:str})
    # If only 3 cols, synthesize RegionLabel from the *gap* coords (less ideal). We prefer 4th col from annotated step.
    if df.shape[1] >= 4:
        df.columns = ["Chrom","Start","End","RegionLabel"] + [f"extra{i}" for i in range(df.shape[1]-4)]
    else:
        df.columns = ["Chrom","Start","End"]
        df["RegionLabel"] = df["Chrom"].astype(str) + ":" + df["Start"].astype(int).astype(str) + "-" + df["End"].astype(int).astype(str)

    # Optional filter by SF genes if RegionLabel encodes gene (it won't here, so skip)
    df["GapLen"] = pd.to_numeric(df["End"], errors="coerce") - pd.to_numeric(df["Start"], errors="coerce")
    df = df.dropna(subset=["Start","End","GapLen"])

    def _intervals(g):
        return ";".join(f"{r.Chrom}:{int(r.Start)}-{int(r.End)}" for _, r in g.sort_values(["Chrom","Start","End"]).iterrows())

    agg = (df.groupby("RegionLabel", as_index=False)
             .agg(**{
                 f"Gaps{threshold_tag}_n": ("GapLen","size"),
                 f"Gaps{threshold_tag}_bp": ("GapLen","sum"),
             }))
    ints = (df.groupby("RegionLabel").apply(_intervals)
             .reset_index(name=f"Gaps{threshold_tag}_intervals"))
    out = agg.merge(ints, on="RegionLabel", how="left")
    out[f"Gaps{threshold_tag}_n"] = out[f"Gaps{threshold_tag}_n"].astype("Int64")
    out[f"Gaps{threshold_tag}_bp"] = out[f"Gaps{threshold_tag}_bp"].astype("Int64")
    return out


def parse_verifybamid(path):
    """Return FREEMIX contamination if VerifyBamID2 present."""
    out = {}
    if not path or not os.path.exists(path): return out
    try:
        with open(path) as f:
            for l in f:
                if "FREEMIX" in l:
                    m = re.search(r"FREEMIX[:=\s]+([\d\.eE+-]+)", l)
                    if m: out["FREEMIX"] = safe_float(m.group(1))
                # selfSM format: first line often: SAMPLE\tFREEMIX\t... ; handle tabular too
                if "\t" in l and l.upper().startswith("SAMPLE"):
                    # next line contains number
                    row = next(f).strip().split("\t")
                    for i,h in enumerate(l.strip().split("\t")):
                        if h.upper()=="FREEMIX":
                            out["FREEMIX"] = safe_float(row[i])
                            break
                    break
    except Exception:
        pass
    return out

def parse_sexcheck(path):
    if not path or not os.path.exists(path): return {}
    txt = open(path).read().strip()
    m = re.search(r"sex\s*[:=]\s*(\w+)", txt, re.IGNORECASE)
    return {"Sex_Check": m.group(1)} if m else {"Sex_Check": txt}

# -------------------------
# Coverage summary loader (your format)
# -------------------------
coverage_data = {}
with open(args.coverage_summary) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 5: continue
        region = parts[0]             # e.g. "1:17018000-17019000"
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

# R1/R2 and FWD/REV
r1r2_data = {}
with open(args.r1r2_summary) as f:
    for line in f:
        chrom, start, end, r1, r2, ratio = line.strip().split("\t")
        r1r2_data[(chrom, int(start), int(end))] = (r1, r2, ratio)

def get_r1r2(chrom, pos):
    for (c,s,e), vals in r1r2_data.items():
        if str(chrom)==c and s<=pos<=e: return vals
    return ("NA","NA","NA")

fr_data = {}
with open(args.fr_summary) as f:
    for line in f:
        chrom, start, end, fwd, rev, ffrac, rfrac = line.strip().split("\t")
        fr_data[(chrom, int(start), int(end))] = (fwd, rev, ffrac, rfrac)

def get_fr(chrom, pos):
    for (c,s,e), vals in fr_data.items():
        if str(chrom)==c and s<=pos<=e: return vals
    return ("NA","NA","NA","NA")

# Optional SF gene list
sf_genes = set()
if args.sf_genes and os.path.exists(args.sf_genes):
    with open(args.sf_genes) as f:
        for l in f:
            g = l.strip().split()[0]
            if g: sf_genes.add(g)

# -------------------------
# Parse VCF to dataframe
# -------------------------
vcf = VCF(args.vcf)
sample_id = args.sample_id or sample_from_vcf(args.vcf)

# CSQ fields (if present)
csq_format = []
for h in vcf.raw_header.split("\n"):
    if h.startswith("##INFO=<ID=CSQ"):
        try:
            csq_format = h.split("Format: ")[1].rstrip('">').split("|")
        except: pass
        break

records = []
for var in vcf:
    chrom, pos, ref = var.CHROM, var.POS, var.REF
    alt = ",".join(var.ALT) if var.ALT else ""
    filt = var.FILTER or "PASS"

    dp = var.format("DP")[0][0] if "DP" in var.FORMAT else var.INFO.get("DP")
    vaf = var.format("VAF")[0][0] if "VAF" in var.FORMAT else None
    gq  = var.format("GQ")[0][0] if "GQ" in var.FORMAT else None

    gt = var.genotypes[0][:2]
    zyg = "het" if gt[0] != gt[1] else ("hom" if gt[0] == 1 else "ref")

    ad_ref = ad_alt = None
    if "AD" in var.FORMAT:
        try: ad_ref, ad_alt = var.format("AD")[0]
        except: pass

    gene = transcript = hgvsc = hgvsp = consequence = exon = intron = impact = None
    clinvar = gnomad_af = revel = spliceai = cadd = None
    ann = {}
    if csq_format:
        csq = var.INFO.get("CSQ")
        if csq:
            first = csq.split(",")[0].split("|")
            ann = dict(zip(csq_format, first))
            gene = ann.get("SYMBOL")
            transcript = ann.get("Feature")
            hgvsc = ann.get("HGVSc")
            hgvsp = ann.get("HGVSp")
            consequence = ann.get("Consequence")
            exon = ann.get("EXON")
            intron = ann.get("INTRON")
            impact = ann.get("IMPACT")
            gnomad_af = ann.get("gnomADg_AF") or ann.get("gnomADe_AF") or ann.get("AF")
            revel = ann.get("REVEL")
            spliceai = ann.get("SpliceAI")
            clinvar = ann.get("CLIN_SIG")
            cadd = ann.get("CADD_PHRED") or ann.get("CADDraw")

    cov = get_exon_cov(chrom, pos)
    r1r2 = get_r1r2(chrom, pos)
    fr = get_fr(chrom, pos)
    try:
        r1 = int(r1r2[0]); r2 = int(r1r2[1])
        r1r2_abs = round(r1 / r2, 2) if r2 > 0 else None
    except: r1r2_abs = None
    try:
        fwd = int(fr[0]); rev = int(fr[1])
        fr_abs = round(fwd / rev, 2) if rev > 0 else None
    except: fr_abs = None

    rec = {
        "Sample": sample_id,
        "Chrom": chrom, "Pos": pos, "Ref": ref, "Alt": alt,
        "Variant": f"{chrom}:{pos}:{ref}:{alt}",
        "FILTER": filt,
        "Gene": gene, "Transcript": transcript,
        "Consequence": consequence, "Impact": impact,
        "Exon": exon, "Intron": intron,
        "HGVSc": hgvsc, "HGVSp": hgvsp,
        "GT": f"{gt[0]}/{gt[1]}",
        "AD_Ref": ad_ref, "AD_Alt": ad_alt,
        "DP": dp, "GQ": gq, "VAF": vaf, "Zygosity": zyg,
        "ExonCov20": cov.get(">=20x","NA"),
        "ExonCov30": cov.get(">=30x","NA"),
        "ExonCov50": cov.get(">=50x","NA"),
        "ExonCov100": cov.get(">=100x","NA"),
        "R1": r1r2[0], "R2": r1r2[1], "R1R2_frac": r1r2[2], "R1R2_abs": r1r2_abs,
        "FWD": fr[0], "REV": fr[1], "FWD_frac": fr[2], "REV_frac": fr[3], "FR_abs": fr_abs,
        "gnomAD_AF": gnomad_af, "REVEL": revel, "SpliceAI": spliceai, "CADD": cadd,
        #"ClinVar": ann.get("CLIN_SIG") if ann else None,
        "ClinVar": clinvar,
        "HGVS_full": None
    }
    # Combined HGVS (Gene c. p.)
    if gene or hgvsc or hgvsp:
        rec["HGVS_full"] = f"{gene or ''} {hgvsc or ''} {hgvsp or ''}".strip()
    records.append(rec)

df = pd.DataFrame(records)

# -------------------------
# Build tabs
# -------------------------
with pd.ExcelWriter(args.xlsx_out) as xw:

    # 1) Sample Summary (single row)
    ss = {
        "Sample_ID": sample_id,
        "Assay": args.assay,
        "Build": args.build
    }
    ss.update(parse_flagstat(args.flagstat))
    ss.update(parse_samtools_stats(args.stats))
    # Insert size & mapping/Q30 from Picard if present
    tmp = parse_picard_alignment(args.picard_align);  ss.update(tmp)
    tmp = parse_picard_insert(args.picard_insert);    ss.update({k:v for k,v in tmp.items() if v is not None})
    # Mosdepth means
    ss.update(parse_mosdepth_summary(args.mosdepth_summary))
    # Contamination & sex
    ss.update(parse_verifybamid(args.verifybamid))
    ss.update(parse_sexcheck(args.sexcheck))

    # Derive Ti/Tv if not provided (rough, SNP only)
    try:
        snps = df[(df["Ref"].str.len()==1) & (df["Alt"].str.len()==1)]
        transitions = {("A","G"),("G","A"),("C","T"),("T","C")}
        ti = sum((r.Ref,r.Alt) in transitions for _,r in snps.iterrows())
        tv = len(snps)-ti
        if tv>0 and "TiTv" not in ss: ss["TiTv"] = round(ti/tv,3)
    except: pass

    # Het/Hom summary from VCF if not in stats
    try:
        het = (df["Zygosity"]=="het").sum()
        hom = (df["Zygosity"]=="hom").sum()
        if hom>0 and "Het_Hom_Ratio" not in ss: ss["Het_Hom_Ratio"] = round(het/hom,3)
    except: pass

    pd.DataFrame([ss]).to_excel(xw, index=False, sheet_name="Sample Summary")

    # 2) ACMG SF (P/LP) — filter by ClinVar and (optionally) SF gene list
    acmg = df.copy()
    acmg = acmg[acmg["ClinVar"].str.contains("pathogenic", case=False, na=False)]
    if sf_genes:
        acmg = acmg[acmg["Gene"].isin(sf_genes)]
    acmg_cols = [
        "Gene","Variant","HGVSc","HGVSp","Zygosity","GT","AD_Ref","AD_Alt","DP","GQ",
        "Consequence","Exon","Intron","ClinVar","gnomAD_AF","REVEL","SpliceAI","CADD","HGVS_full"
    ]
    acmg[acmg_cols].to_excel(xw, index=False, sheet_name="ACMG SF (P-LP)")

    # 3) Coverage gaps in ACMG SF genes (or all genes if no list)
    combined_df = pd.DataFrame(columns=[
        "Gene","Transcript","Exon","RegionLabel","Chrom","ExonStart","ExonEnd","ExonLen",
        "Pct>=20x","Pct>=30x",
        "Gaps<20x_n","Gaps<20x_bp","Gaps<20x_intervals",
        "Gaps<30x_n","Gaps<30x_bp","Gaps<30x_intervals"
    ])

    try:
        th_exon = read_acmg_thresholds(args.acmg_thresholds, sf_genes=sf_genes)  # RegionLabel = chrom:start-end if no names
        g20_agg = load_gap_bed_for_agg(args.gaps20, "<20x", sf_genes=sf_genes)
        g30_agg = load_gap_bed_for_agg(args.gaps30, "<30x", sf_genes=sf_genes)

        combined_df = (th_exon
                       .merge(g20_agg, on="RegionLabel", how="left")
                       .merge(g30_agg, on="RegionLabel", how="left"))

        # ensure all expected columns exist
        for c in ["Gaps<20x_n","Gaps<20x_bp","Gaps<20x_intervals","Gaps<30x_n","Gaps<30x_bp","Gaps<30x_intervals"]:
            if c not in combined_df.columns: combined_df[c] = pd.NA

        combined_df = (combined_df[
            ["Gene","Transcript","Exon","RegionLabel","Chrom","ExonStart","ExonEnd","ExonLen",
             "Pct>=20x","Pct>=30x",
             "Gaps<20x_n","Gaps<20x_bp","Gaps<20x_intervals",
             "Gaps<30x_n","Gaps<30x_bp","Gaps<30x_intervals"]]
            .sort_values(["Chrom","ExonStart","ExonEnd"])
        )
    except Exception:
        pass

    combined_df.to_excel(xw, index=False, sheet_name="Coverage gaps")


    # 4) PASS variant table
    pass_cols = [
        "Variant","Gene","Transcript","Consequence","GT","AD_Ref","AD_Alt","DP","GQ",
        "FILTER","VAF","gnomAD_AF","ClinVar","REVEL","SpliceAI","CADD","HGVS_full"
    ]
    df[df["FILTER"]=="PASS"][pass_cols].to_excel(xw, index=False, sheet_name="PASS variants")

    # Your existing tabs
    #df.to_excel(xw, index=False, sheet_name="All Variants")

    #high_conf = df[
    #    ((df["Zygosity"]=="het") & df["VAF"].between(0.35,0.65, inclusive="both")) |
    #    ((df["Zygosity"]=="hom") & (df["VAF"]>=0.90))
    #]
    #high_conf.to_excel(xw, index=False, sheet_name="High Confidence")

    #medium_conf = df[ df["VAF"].between(0.10,0.34) | df["VAF"].between(0.65,0.89) ]
    #medium_conf.to_excel(xw, index=False, sheet_name="Medium Confidence")

    #low_conf = df[ df["VAF"] < 0.20 ]
    #low_conf.to_excel(xw, index=False, sheet_name="Low Confidence")

    # Coverage summary (from per-exon fields already in df)
    #cov_cols = ["ExonCov20","ExonCov30","ExonCov50","ExonCov100"]
    #for c in cov_cols: df[c] = pd.to_numeric(df[c], errors="coerce")
    #df[cov_cols].describe().transpose().to_excel(xw, sheet_name="Coverage Summary")

    #gene_summary = (df
    #    .groupby("Gene", dropna=True)
    #    .agg(Variants=("Variant","count"),
    #         VAF_mean=("VAF","mean"),
    #         VAF_min=("VAF","min"),
    #         VAF_max=("VAF","max"),
    #         Consequences=("Consequence", lambda x: ",".join(sorted(set([y for y in x.dropna()])))))
    #    .sort_values("Variants", ascending=False))
    #gene_summary.to_excel(xw, sheet_name="Gene Summary")

print(f"Wrote Excel report → {args.xlsx_out}")
