#!/usr/bin/env python3

import argparse
import os
import subprocess
import hail as hl

#all of us paths
MT_PATH  = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/splitMT/hail.mt"
ANC_PATH = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"

#pops to run (META = all samples)
POPS = {
    "META": None,
#    "EUR":  "eur",
#    "AFR":  "afr",
#    "AMR":  "amr",
}

def sh(cmd):
    print("\n$ " + " ".join(cmd))
    subprocess.check_call(cmd)

def normalize_chr(c: str) -> str:
    c = str(c)
    if c in {"X", "23", "chr23"}:
        return "chrX"
    if not c.startswith("chr"):
        return "chr" + c
    return c

def main():
    ap = argparse.ArgumentParser(description="AoU: export region PLINK + compute LD (gz), subset to sumstats IDs.")
    ap.add_argument("--gene", required=True, help="Label for output (e.g., IL4R)")
    ap.add_argument("--phecode", required=True, help="Phenotype code used in sumstats filename (e.g., EM_239)")
    ap.add_argument("--chr", required=True, help="Chromosome (e.g., chr1, 1, chrX, X, 23)")
    ap.add_argument("--left", type=int, required=True, help="Left bp (inclusive)")
    ap.add_argument("--right", type=int, required=True, help="Right bp (inclusive)")
    ap.add_argument("--threads", type=int, default=8, help="PLINK threads (default 8)")
    ap.add_argument("--max_samples", type=int, default=50000, help="Downsample to at most this many samples per pop (default 50000)")
    args = ap.parse_args()

    bucket = os.getenv("WORKSPACE_BUCKET")
    if not bucket:
        raise SystemExit("ERROR: WORKSPACE_BUCKET is not set.")
    gproj = os.getenv("GOOGLE_PROJECT")
    if not gproj:
        raise SystemExit("ERROR: GOOGLE_PROJECT is not set (needed for requester-pays gsutil).")

    chrom = normalize_chr(args.chr)
    left, right = args.left, args.right

    def ensure_hail_initialized(**kwargs):
        try:
            hl.current_backend()
        except Exception:
            hl.init(**kwargs)

    ensure_hail_initialized()

    #load MT + ancestry and join
    mt = hl.read_matrix_table(MT_PATH)

    anc = hl.import_table(ANC_PATH, impute=True).key_by("research_id")  # research_id is int32
    mt = mt.annotate_cols(anc=anc[hl.int32(mt.s)])                      # mt.s is numeric string

    #build MT row ID to match sumstats "ID" like chr1:10857561:C:T (chrX supported)
    mt = mt.annotate_rows(
        varid=hl.format(
            "%s:%d:%s:%s",
            mt.locus.contig,
            mt.locus.position,
            mt.alleles[0],
            mt.alleles[1],
        )
    )

    for pop, anc_code in POPS.items():
        print(f"\n==== {pop} ====")

        #filter samples
        mt_pop = mt if anc_code is None else mt.filter_cols(mt.anc.ancestry_pred == anc_code)

        # downsample
        n0 = mt_pop.count_cols()
        print("Samples (before downsample):", n0)
        if args.max_samples and n0 > args.max_samples:
            frac = args.max_samples / n0
            mt_pop = mt_pop.sample_cols(frac, seed=1)
            print("Samples (after downsample):", mt_pop.count_cols())
        else:
            print("Samples (after downsample):", n0)

        #region filter
        mt_reg = mt_pop.filter_rows(
            (mt_pop.locus.contig == chrom) &
            (mt_pop.locus.position >= left) &
            (mt_pop.locus.position <= right)
        )

        # ---- CHANGED: import only needed sumstats fields, filter to region BEFORE keying ----
        ss_path = f"{bucket}/data/{pop}_coloc_mesa_{args.phecode}.tsv"
        print("Sumstats:", ss_path)

        ss = hl.import_table(
            ss_path,
            delimiter="\t",
            impute=False,
            types={"ID": hl.tstr, "CHR": hl.tstr, "POS": hl.tint32},
            missing=""
        )

        # filter sumstats to this region first (reduces shuffle a ton)
        ss = ss.filter((ss.CHR == chrom) & (ss.POS >= left) & (ss.POS <= right))

        # now key only the (much smaller) filtered table
        ss = ss.select("ID").key_by("ID")
        # -------------------------------------------------------------------------------

        #subset MT variants to those present in sumstats IDs
        mt_reg = mt_reg.filter_rows(hl.is_defined(ss[mt_reg.varid]))
        print("Variants (region ∩ sumstats):", mt_reg.count_rows())

        #export PLINK to bucket
        prefix = f"{pop}_{args.gene}_{args.phecode}_{chrom}_{left}_{right}"
        out_gcs = f"{bucket}/data/{prefix}"
        print("Exporting PLINK:", out_gcs)
        hl.export_plink(mt_reg, out_gcs)

        #copy locally
        sh(["gsutil", "-u", gproj, "-m", "cp", f"{out_gcs}.bed", f"{out_gcs}.bim", f"{out_gcs}.fam", "."])

        #compute dense LD
        sh([
            "plink",
            "--threads", str(args.threads),
            "--bfile", prefix,
            "--r", "square",
            "--out", f"{prefix}_LD"
        ])

        #always gzip LD and upload
        sh(["gzip", "-f", f"{prefix}_LD.ld"])
        sh(["gsutil", "-u", gproj, "cp", f"{prefix}_LD.ld.gz", f"{bucket}/data/"])

        # Local cleanup
        for fn in [
            f"{prefix}.bed", f"{prefix}.bim", f"{prefix}.fam",
            f"{prefix}.log", f"{prefix}.nosex",
            f"{prefix}_LD.log", f"{prefix}_LD.ld.gz"
        ]:
            if os.path.exists(fn):
                os.remove(fn)

        print(f"Done: {bucket}/data/{prefix}_LD.ld.gz")

if __name__ == "__main__":
    main()
