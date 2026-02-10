#!/usr/bin/env bash
set -euo pipefail

gsutil cp gs://fc-secure-915a27e9-c961-4aa0-b53f-09825278579c/data/validated_proteins.tsv /home/jupyter

#set files
IN="/home/jupyter/validated_proteins.tsv"
SCRIPT="$HOME/mesa_coloc/03coloc.R"
OUT="/home/jupyter/coloc_output.txt"

#initialize output
rm -f "$OUT"

#set column names to index
PHECODE_COLNAME="aou_code"
# MESA phenotype id column (should be something like OIDxxxxx)
OID_COLNAME="OlinkID"
# MESA population id
DBPOP_COLNAME="db_pop"
# AoU population id
GWASPOP_COLNAME="gwas_pop_aou"
# -----------------------------------------------

#find column indices by header name (tab-separated)
get_col() {
  local name="$1"
  awk -v FS='\t' -v colname="$name" '
    NR==1{
      for(i=1;i<=NF;i++){
        if($i==colname){ print i; exit }
      }
      exit 1
    }' "$IN"
}

PHECODE_I=$(get_col "$PHECODE_COLNAME") || { echo "Missing column: $PHECODE_COLNAME"; exit 1; }
OID_I=$(get_col "$OID_COLNAME")         || { echo "Missing column: $OID_COLNAME"; exit 1; }
DBPOP_I=$(get_col "$DBPOP_COLNAME")     || { echo "Missing column: $DBPOP_COLNAME"; exit 1; }
GWASPOP_I=$(get_col "$GWASPOP_COLNAME") || { echo "Missing column: $GWASPOP_COLNAME"; exit 1; }

echo "Using columns:"
echo "  phecode  ($PHECODE_COLNAME) -> $PHECODE_I"
echo "  oid      ($OID_COLNAME)     -> $OID_I"
echo "  db_pop   ($DBPOP_COLNAME)   -> $DBPOP_I"
echo "  gwas_pop ($GWASPOP_COLNAME) -> $GWASPOP_I"
echo

#loop over rows (skip header)
awk -v FS='\t' -v phe_i="$PHECODE_I" -v oid_i="$OID_I" -v db_i="$DBPOP_I" -v gwas_i="$GWASPOP_I" '
  NR>1{
    print $phe_i "\t" $oid_i "\t" $db_i "\t" $gwas_i
  }' "$IN" \
| while IFS=$'\t' read -r phecode OlinkID db_pop gwas_pop; do

    # skip empty lines
    [[ -z "${phecode// }" ]] && continue

    echo "Running: phecode = $phecode, oid = $OlinkID, db_pop = $db_pop, gwas_pop = $gwas_pop"

    Rscript ~/mesa_coloc/03coloc.R \
      --phecode "$phecode" \
      --oid "$OlinkID" \
      --db_pop "$db_pop" \
      --gwas_pop "$gwas_pop"

done
