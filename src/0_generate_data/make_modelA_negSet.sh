#!/usr/bin/env bash
set -euo pipefail

# ============================================================
#   Make non-promoter windows (CLI version)
#   Usage:
#     bash make_nonpromoters_cli.sh \
#       --workdir /path/to/raw1 \
#       --species human,mouse_mm39 \
#       --windows w2k1k,w1k1k,w3k1k
# ============================================================

# ---------- Default values ----------
WORKDIR=""
OUTDIR=""
SPECIES="human,mouse_mm39"
WINDOWS="w2k1k,w1k1k,w3k1k"
# ====================================

# ---------- Parse arguments ----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --workdir)  WORKDIR="$2"; shift 2 ;;
    --outdir)   OUTDIR="$2"; shift 2 ;;
    --species)  SPECIES="$2"; shift 2 ;;
    --windows)  WINDOWS="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: bash make_nonpromoters_cli.sh --workdir <path> [--outdir <path>] [--species human,mouse_mm39] [--windows w2k1k,w1k1k,w3k1k]"
      exit 0 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# ---------- Sanity checks ----------
if [[ -z "$WORKDIR" ]]; then
  echo "ERROR: --workdir is required."
  exit 1
fi

OUTDIR="${OUTDIR:-$WORKDIR}"
mkdir -p "$WORKDIR" "$OUTDIR"
cd "$WORKDIR"

# ---------- Configuration ----------
IFS=',' read -r -a SPECIES_ARR <<< "$SPECIES"
IFS=',' read -r -a WINS <<< "$WINDOWS"

declare -A WIN_UP=( ["w2k1k"]=2000 ["w1k1k"]=1000 ["w3k1k"]=3000 )
declare -A WIN_DN=( ["w2k1k"]=1000 ["w1k1k"]=1000 ["w3k1k"]=1000 )

HG38_GENOME="$WORKDIR/hg38.genome"
MM39_GENOME="$WORKDIR/mm39.genome"

# ---------- Functions ----------
download_or_fetch_sizes () {
  local asm="$1" out="$2"
  echo ">> Preparing chrom sizes for $asm → $out"
  if command -v fetchChromSizes >/dev/null 2>&1; then
    fetchChromSizes "$asm" > "$out.tmp"
  else
    case "$asm" in
      hg38)
        curl -L -o "$out.tmp" "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes" ;;
      mm39)
        curl -L -o "$out.tmp" "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes" ;;
      *) echo "Unknown assembly: $asm" >&2; exit 1 ;;
    esac
  fi

  if [[ "$asm" == "hg38" ]]; then
    awk '$1 ~ /^chr([0-9]+|X|Y)$/' "$out.tmp" > "$out"
  else
    awk '$1 ~ /^chr([0-9]+|X|Y)$/' "$out.tmp" \
      | awk '$1 !~ /^chr([2-9][0-9]+)$/' > "$out"
  fi
  rm -f "$out.tmp"
  echo "   chroms kept → $(wc -l < "$out")"
}

[[ -s "$HG38_GENOME" ]] || download_or_fetch_sizes hg38 "$HG38_GENOME"
[[ -s "$MM39_GENOME" ]] || download_or_fetch_sizes mm39 "$MM39_GENOME"

genome_sizes_for () {
  case "$1" in
    human) echo "$HG38_GENOME" ;;
    mouse_mm39) echo "$MM39_GENOME" ;;
    *) echo "Unknown species: $1" >&2; exit 1 ;;
  esac
}

# ---------- Main logic ----------
command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found in PATH"; exit 1; }

echo "== Running make_nonpromoters_cli.sh =="
echo "WORKDIR : $WORKDIR"
echo "OUTDIR  : $OUTDIR"
echo "SPECIES : ${SPECIES_ARR[*]}"
echo "WINDOWS : ${WINS[*]}"
echo

for sp in "${SPECIES_ARR[@]}"; do
  gsizes="$(genome_sizes_for "$sp")"
  [[ -s "$gsizes" ]] || { echo "Missing genome sizes: $gsizes"; exit 1; }
  echo "== $sp (sizes: $gsizes) =="

  for win in "${WINS[@]}"; do
    up=${WIN_UP[$win]}
    dn=${WIN_DN[$win]}
    size=$(( up + dn ))

    prom="${sp}_${win}_promoter.bed"
    win6="${OUTDIR}/${sp}_${win}_windows.bed6"
    nonp="${OUTDIR}/${sp}_${win}_nonpromoter_wins.bed6"

    if [[ ! -s "$prom" ]]; then
      echo "  ! Skipping $sp/$win — promoter BED not found: $prom"
      continue
    fi

    echo "  - $sp / $win (window ${size}bp; up=${up}, down=${dn})"
    bedtools makewindows -g "$gsizes" -w "$size" \
      | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,".",0,"."}' > "$win6"

    bedtools subtract -a "$win6" -b "$prom" > "$nonp"

    echo "    windows:   $win6 ($(wc -l < "$win6"))"
    echo "    promoter:  $prom ($(wc -l < "$prom"))"
    echo "    non-prom:  $nonp ($(wc -l < "$nonp"))"
  done
done

echo "✓ Done. Outputs written to $OUTDIR"
