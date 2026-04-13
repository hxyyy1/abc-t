#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ABC_BIN="${ABC_BIN:-$ROOT_DIR/abc}"
LIB_FILE="${LIB_FILE:-${1:-$ROOT_DIR/libs/sky130.genlib}}"
BENCH_PATH="${BENCH_PATH:-${2:-$ROOT_DIR/benchmarks/EPFL}}"
OUT_DIR="${OUT_DIR:-${3:-$ROOT_DIR/results}}"

CMD1="${CMD1:-map}"
CMD2="${CMD2:-map -j -K 24}"
CMD1_LABEL="${CMD1_LABEL:-$CMD1}"
CMD2_LABEL="${CMD2_LABEL:-$CMD2}"
BENCH_GLOB="${BENCH_GLOB:-*.aig}"
ABC_TIMEOUT_SECS="${ABC_TIMEOUT_SECS:-120}"

mkdir -p "$OUT_DIR"

log_str() {
    local s="${1:-}"
    if [[ -z "$s" ]]; then
        printf 'unknown'
    else
        printf '%s' "$s"
    fi
}

slugify() {
    local s
    s="$(printf '%s' "$1" | tr '[:upper:]' '[:lower:]')"
    s="$(printf '%s' "$s" | sed -E 's/[^a-z0-9]+/_/g; s/^_+//; s/_+$//; s/_+/_/g')"
    if [[ -z "$s" ]]; then
        s="cmd"
    fi
    printf '%s' "$s"
}

detect_lib_read_cmd() {
    case "${LIB_FILE##*.}" in
        lib|liberty)
            printf 'read_lib'
            ;;
        genlib|gen|g)
            printf 'read_library'
            ;;
        *)
            printf 'read_lib'
            ;;
    esac
}

collect_benches() {
    if [[ -f "$BENCH_PATH" ]]; then
        printf '%s\n' "$BENCH_PATH"
    elif [[ -d "$BENCH_PATH" ]]; then
        find "$BENCH_PATH" -maxdepth 1 -type f -name "$BENCH_GLOB" | sort
    else
        printf 'Benchmark path not found: %s\n' "$BENCH_PATH" >&2
        return 1
    fi
}

LIB_READ_CMD="${LIB_READ_CMD:-$(detect_lib_read_cmd)}"
LIB_TAG="$(basename "$LIB_FILE")"
BENCH_TAG="$(basename "$BENCH_PATH")"
CMD1_SLUG="$(slugify "$CMD1_LABEL")"
CMD2_SLUG="$(slugify "$CMD2_LABEL")"
OUT_PREFIX="${OUT_PREFIX:-compare_${BENCH_TAG}_${LIB_TAG%.*}_${CMD1_SLUG}_vs_${CMD2_SLUG}}"
CSV_OUT="$OUT_DIR/${OUT_PREFIX}.csv"
MD_OUT="$OUT_DIR/${OUT_PREFIX}.md"
TMP_BASE="${TMPDIR:-/tmp}"
TMP_CSV="$(mktemp "$TMP_BASE/${OUT_PREFIX}.csv.XXXXXX")"
TMP_MD="$(mktemp "$TMP_BASE/${OUT_PREFIX}.md.XXXXXX")"

cleanup() {
    rm -f "$TMP_CSV" "$TMP_MD"
}
trap cleanup EXIT

run_case() {
    local bench="$1"
    local cmd="$2"
    local raw line parsed area delay status

    set +e
    if [[ "$ABC_TIMEOUT_SECS" =~ ^[0-9]+$ ]] && [[ "$ABC_TIMEOUT_SECS" -gt 0 ]]; then
        raw="$(timeout "${ABC_TIMEOUT_SECS}s" "$ABC_BIN" -c "$LIB_READ_CMD $LIB_FILE; r $bench; $cmd; ps" 2>&1)"
    else
        raw="$("$ABC_BIN" -c "$LIB_READ_CMD $LIB_FILE; r $bench; $cmd; ps" 2>&1)"
    fi
    status=$?
    set -e

    raw="$(printf '%s\n' "$raw" | sed -E 's/\x1b\[[0-9;]*m//g')"
    if [[ "$status" -ne 0 ]]; then
        printf 'NA,NA,exit_%s\n' "$status"
        return 0
    fi

    line="$(printf '%s\n' "$raw" | grep ' area =' | tail -n 1 || true)"
    if [[ -z "$line" ]]; then
        printf 'NA,NA,no_stats\n'
        return 0
    fi

    parsed="$(printf '%s\n' "$line" | sed -nE 's/.*area = *([0-9.]+).*delay = *([0-9.]+).*/\1,\2/p')"
    if [[ -z "$parsed" ]]; then
        printf 'NA,NA,parse_error\n'
        return 0
    fi

    IFS=',' read -r area delay <<<"$parsed"
    printf '%s,%s,ok\n' "$area" "$delay"
}

printf 'benchmark,%s_area,%s_delay,%s_status,%s_area,%s_delay,%s_status,area_diff,area_diff_pct,delay_diff,delay_diff_pct\n' \
    "$CMD1_SLUG" "$CMD1_SLUG" "$CMD1_SLUG" "$CMD2_SLUG" "$CMD2_SLUG" "$CMD2_SLUG" > "$TMP_CSV"

while IFS= read -r bench; do
    name="$(basename "${bench%.*}")"
    printf 'Running %s\n' "$name" >&2
    IFS=',' read -r area1 delay1 status1 <<<"$(run_case "$bench" "$CMD1")"
    IFS=',' read -r area2 delay2 status2 <<<"$(run_case "$bench" "$CMD2")"

    awk -v name="$name" -v s1="$status1" -v s2="$status2" \
        -v a1="$area1" -v d1="$delay1" \
        -v a2="$area2" -v d2="$delay2" '
        BEGIN {
            if (s1 == "ok" && s2 == "ok") {
                ad = a2 - a1;
                ap = (a1 == 0 ? 0 : ad * 100.0 / a1);
                dd = d2 - d1;
                dp = (d1 == 0 ? 0 : dd * 100.0 / d1);
                printf "%s,%.2f,%.2f,%s,%.2f,%.2f,%s,%.2f,%.2f,%.2f,%.2f\n",
                    name, a1, d1, s1, a2, d2, s2, ad, ap, dd, dp;
            } else {
                printf "%s,%s,%s,%s,%s,%s,%s,NA,NA,NA,NA\n",
                    name, a1, d1, s1, a2, d2, s2;
            }
        }' >> "$TMP_CSV"
done < <(collect_benches)

awk -F, -v cmd1="$CMD1_LABEL" -v cmd2="$CMD2_LABEL" \
    -v lib_file="$LIB_FILE" -v lib_cmd="$LIB_READ_CMD" -v bench_path="$BENCH_PATH" '
BEGIN {
    print "| benchmark | " cmd1 " area | " cmd1 " delay | " cmd1 " status | " cmd2 " area | " cmd2 " delay | " cmd2 " status | area diff | area diff % | delay diff | delay diff % |";
    print "| --- | ---: | ---: | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: |";
}
NR == 1 { next }
{
    if ($4 == "ok" && $7 == "ok") {
        sum_a1 += $2; sum_d1 += $3; sum_a2 += $5; sum_d2 += $6;
        sum_ad += $8; sum_ap += $9; sum_dd += $10; sum_dp += $11;
        n++;
    }
    printf "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n",
        $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11;
}
END {
    if (n > 0) {
        printf "| **avg (success only)** | %.2f | %.2f | ok | %.2f | %.2f | ok | %.2f | %.2f%% | %.2f | %.2f%% |\n",
            sum_a1 / n, sum_d1 / n, sum_a2 / n, sum_d2 / n,
            sum_ad / n, sum_ap / n, sum_dd / n, sum_dp / n;
    }
}' "$TMP_CSV" > "$TMP_MD"

mv "$TMP_CSV" "$CSV_OUT"
mv "$TMP_MD" "$MD_OUT"
trap - EXIT

printf 'Config: lib_read_cmd=%s lib_file=%s bench_path=%s cmd1=%s cmd2=%s\n' \
    "$(log_str "$LIB_READ_CMD")" "$(log_str "$LIB_FILE")" "$(log_str "$BENCH_PATH")" \
    "$(log_str "$CMD1")" "$(log_str "$CMD2")"
printf 'Wrote %s\n' "$CSV_OUT"
printf 'Wrote %s\n' "$MD_OUT"
