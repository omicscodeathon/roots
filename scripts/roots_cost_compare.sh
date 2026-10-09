#!/usr/bin/env bash
# roots_cost_compare.sh -- live, side-by-side Hetzner vs GPU.ai cost table for the ROOTS
# pipeline's remaining steps (MD extension/convergence, MM-PBSA, trajectory analysis).
#
# RUN THIS ON YOUR LOCAL MACHINE, not in a cloud sandbox: it needs real internet egress to
# api.hetzner.cloud (needs a read-only HCLOUD_TOKEN) and to api.gpu.ai via the `gpu` CLI (needs
# GPUAI_API_KEY, same as gpuai_resume.sh). It does NOT create, launch, or bill anything -- every
# call below is a read-only pricing/catalog lookup.
#
# WHAT IT COSTS OUT (three steps, matching gpuai_resume.sh's own defaults/comments):
#   1. MD extend-to-convergence  -- 4 systems (Sys1-4), MIN_NS=50 floor, EXTEND_NS=10 per
#      round until convergence. No one knows the real number of extension rounds up front
#      (that's what "extend until it ACTUALLY CONVERGES" means) -- so this reports
#      cost PER 10ns-EXTENSION-ROUND-ACROSS-ALL-4-SYSTEMS, not a total, and you multiply by
#      however many rounds your own RMSD trace needs. Pass ROUNDS=N to also print an N-round
#      total.
#   2. MM-PBSA on the converged frames -- CPU-bound, no GPU needed; costed on a CPU-only
#      instance class on each provider.
#   3. Trajectory analysis/visualization (MDAnalysis/PLIP/matplotlib/Py3Dmol) -- also
#      CPU-bound, no GPU needed; same CPU-only class.
#
# Both #2 and #3 are included ONLY because you said you want them priced on rented
# infrastructure for workflow consistency -- they do not need a GPU and would cost $0 run
# locally. That tradeoff is flagged in the table output, not hidden.
#
# Hours per step are ESTIMATES you can override via env vars (see ESTIMATE DEFAULTS below);
# this script does not know your actual GPU/CPU wall-clock, only the sizing already encoded
# in gpuai_resume.sh (4 systems, 50ns floor, 10ns extend rounds, per-second GPU.ai billing,
# hourly-rounded Hetzner billing).
set -euo pipefail

need_bin() { command -v "$1" >/dev/null || { echo "'$1' not found in PATH. $2" >&2; exit 1; }; }
need_bin curl "Install curl."
need_bin jq "Install it: brew install jq / sudo apt-get install -y jq"
need_bin awk "awk should ship with your OS; something is very wrong if it's missing."

# ---------------- ESTIMATE DEFAULTS (override via env) ----------------
# MD extend round: wall-clock hours for ONE 10ns EXTEND_NS round across all 4 systems combined,
# sharing one instance's GPUs. Set from your own last segment's actual wall time if you have it
# (check: tail -1 each Sys*/md.log's "Time:" line from a completed EXTEND_NS=10 run) -- this
# default (0.6h) is only a placeholder for an A100-class GPU doing ~4 systems/8 GPUs-worth of
# throughput; CORRECT IT with ROUND_HOURS=<your number> once you have one real data point.
ROUND_HOURS="${ROUND_HOURS:-0.6}"
ROUNDS="${ROUNDS:-0}"                 # 0 = don't print a multi-round total, just per-round cost

# MM-PBSA: CPU-bound, hours for all 4 systems' converged-frame sets combined.
MMPBSA_HOURS="${MMPBSA_HOURS:-4}"

# Trajectory analysis/visualization: CPU-bound, hours for all 4 systems combined.
ANALYSIS_HOURS="${ANALYSIS_HOURS:-3}"

# ---------------- GPU.ai live pricing ----------------
GPUAI_GPU_TYPE="${GPUAI_TYPES:-a100_80gb}"
# No GPUAI_CPU_TYPE default guessed -- GPU.ai's CPU-only catalog id is discovered live below
# from `gpu gpu-types list` (or whatever the equivalent "list types" subcommand is called on
# your CLI version) rather than assumed, since "cpu_small" turned out not to exist on this
# account's catalog.
GPUAI_CPU_TYPE="${GPUAI_CPU_TYPE:-}"

echo "== Querying GPU.ai live pricing (gpu CLI, needs GPUAI_API_KEY) ==" >&2

gpuai_discover_cpu_type() {
  # Tries a few plausible subcommands/field names since the exact catalog shape isn't
  # documented here; first one that returns a non-GPU type id wins. Never falls back to a
  # guessed literal -- if none of these parse, GPUAI_CPU_TYPE stays empty and the CPU rows
  # report "unknown" honestly.
  local out
  out=$(gpu gpu-types list --output json 2>/dev/null || true)
  if [[ -n "$out" ]]; then
    jq -r '[.[]? | select((.id // .name // "" ) | test("cpu"; "i"))] | .[0] | (.id // .name // empty)' <<<"$out" 2>/dev/null
    return
  fi
  out=$(gpu instance-types list --output json 2>/dev/null || true)
  if [[ -n "$out" ]]; then
    jq -r '[.[]? | select((.id // .name // "" ) | test("cpu"; "i"))] | .[0] | (.id // .name // empty)' <<<"$out" 2>/dev/null
    return
  fi
}

gpuai_rate() {  # $1 = type id -> prints USD/hr or "unknown (...)"
  local type="$1" out rate
  if ! command -v gpu >/dev/null; then echo "unknown (gpu CLI not installed)"; return; fi
  [[ -n "$type" ]] || { echo "unknown (no matching type id found -- run: gpu gpu-types list)"; return; }
  out=$(gpu pricing --gpu-type "$type" --output json 2>/dev/null || true)
  rate=$(jq -r '[.[]?.price_per_hour // .[]?.price // empty] | map(tonumber? // empty) | min // empty' <<<"$out" 2>/dev/null || true)
  [[ -n "$rate" ]] && echo "$rate" || echo "unknown (run: gpu pricing --gpu-type $type)"
}

GPUAI_GPU_RATE=$(gpuai_rate "$GPUAI_GPU_TYPE")
if [[ -z "$GPUAI_CPU_TYPE" ]]; then
  GPUAI_CPU_TYPE=$(gpuai_discover_cpu_type)
  [[ -n "$GPUAI_CPU_TYPE" && "$GPUAI_CPU_TYPE" != "null" ]] || GPUAI_CPU_TYPE=""
fi
GPUAI_CPU_RATE=$(gpuai_rate "$GPUAI_CPU_TYPE")
GPUAI_CPU_LABEL="${GPUAI_CPU_TYPE:-CPU-only (catalog id not found)}"

# ---------------- Hetzner live pricing ----------------
# Hetzner's STABLE, documented, token-authenticated Cloud API pricing endpoint. (The earlier
# version of this script used an unauthenticated configurator JSON the website itself calls --
# that URL is not part of Hetzner's documented API surface and came back unreachable, so it's
# dropped in favor of the one Hetzner actually commits to keeping working: GET /v1/pricing on
# the Cloud API. This needs a Hetzner API token -- the token only needs READ access, it is
# never used to create or bill anything here.)
#   Get one: Hetzner Cloud Console -> your project -> Security -> API Tokens -> Generate (Read-only is enough)
#   export HCLOUD_TOKEN=your_token_here
HETZNER_API="https://api.hetzner.cloud/v1/pricing"

echo "== Querying Hetzner live pricing (Cloud API /v1/pricing, needs HCLOUD_TOKEN) ==" >&2

HETZNER_JSON=""
if [[ -n "${HCLOUD_TOKEN:-}" ]]; then
  HETZNER_JSON=$(curl -fsSL -H "Authorization: Bearer $HCLOUD_TOKEN" "$HETZNER_API" 2>/dev/null || true)
else
  echo "   (HCLOUD_TOKEN not set -- Hetzner rows will show 'unknown'. Get a read-only token from" >&2
  echo "    the Cloud Console: Project > Security > API Tokens, then: export HCLOUD_TOKEN=...)" >&2
fi

hetzner_rate_hourly() {  # $1 = server type name, e.g. "cpx31" or "gex44" -> USD/hr or "unknown (...)"
  local needle="$1" hourly
  [[ -n "$HETZNER_JSON" ]] || { echo "unknown (HCLOUD_TOKEN not set or /v1/pricing unreachable)"; return; }
  hourly=$(jq -r --arg n "$needle" '
    .pricing.server_types[]? | select(.name? and (.name | test($n; "i")))
    | .prices[]?.price_hourly.net // empty
  ' <<<"$HETZNER_JSON" 2>/dev/null | head -1 || true)
  if [[ -n "$hourly" && "$hourly" != "null" ]]; then
    echo "$hourly"
  else
    echo "unknown (type '$needle' not in /v1/pricing response -- list real names with: curl -s -H \"Authorization: Bearer \$HCLOUD_TOKEN\" $HETZNER_API | jq -r '.pricing.server_types[].name')"
  fi
}

# Candidate Hetzner lines: GEX for GPU (dedicated GPU servers), CPX for shared-vCPU CPU-only.
# You said "any that is available, cost effective, and equally fast" -- so this tries the
# cheapest plausible match in each class rather than assuming a specific size. If GEX isn't
# actually orderable as on-demand hourly in your account/region, this will say so via "unknown"
# rather than inventing a number.
HETZNER_GPU_RATE=$(hetzner_rate_hourly "gex")
HETZNER_CPU_RATE=$(hetzner_rate_hourly "cpx31")

# ---------------- cost math ----------------
# GPU.ai bills per-second, no rounding (see gpuai_resume.sh header, point 2): cost = hours * rate
# Hetzner's /v1/pricing price_hourly is already a true hourly rate with no monthly conversion
# needed; on-demand hourly billing still rounds UP to the next full hour: cost = ceil(hours) * rate
is_number() { [[ "$1" =~ ^[0-9]+(\.[0-9]+)?$ ]]; }

cost_gpuai() {
  local h="$1" r="$2"
  is_number "$r" || { echo "unknown"; return; }
  awk -v h="$h" -v r="$r" 'BEGIN{printf "%.2f", h*r}'
}
cost_hetzner() {
  local h="$1" r="$2"
  is_number "$r" || { echo "unknown"; return; }
  awk -v h="$h" -v r="$r" 'BEGIN{hh=(h==int(h))?h:int(h)+1; printf "%.2f", hh*r}'
}

fmt_rate() { is_number "$1" && echo "\$$1" || echo "$1"; }
fmt_cost() { [[ "$1" == "unknown" ]] && echo "unknown" || echo "\$$1"; }

round_total_note=""
if [[ "$ROUNDS" -gt 0 ]]; then
  round_total_note=" (x$ROUNDS rounds = $(awk -v h="$ROUND_HOURS" -v n="$ROUNDS" 'BEGIN{printf "%.1f", h*n}')h total)"
fi

cat <<EOF

ROOTS pipeline -- Hetzner vs GPU.ai, per step (live rates queried just now)
=============================================================================
Rates found:
  GPU.ai  $GPUAI_GPU_TYPE (GPU) ............ $(fmt_rate "$GPUAI_GPU_RATE")/hr
  GPU.ai  $GPUAI_CPU_LABEL (CPU-only) ....... $(fmt_rate "$GPUAI_CPU_RATE")/hr
  Hetzner gex* line (GPU) .................. $(fmt_rate "$HETZNER_GPU_RATE")/hr  (per-second-rounded-to-hour on-demand billing)
  Hetzner cpx31 (CPU-only) ................ $(fmt_rate "$HETZNER_CPU_RATE")/hr  (same billing note)

EOF

printf "%-42s | %-12s | %-18s | %-12s | %-18s\n" "Step" "Hetzner \$/hr" "Hetzner cost" "GPU.ai \$/hr" "GPU.ai cost"
printf '%s\n' "---------------------------------------------------------------------------------------------------------------------"
printf "%-42s | %-12s | %-18s | %-12s | %-18s\n" \
  "1. MD extend round (10ns x4 sys)$round_total_note" \
  "$(fmt_rate "$HETZNER_GPU_RATE")" "$(fmt_cost "$(cost_hetzner "$ROUND_HOURS" "$HETZNER_GPU_RATE")")" \
  "$(fmt_rate "$GPUAI_GPU_RATE")" "$(fmt_cost "$(cost_gpuai "$ROUND_HOURS" "$GPUAI_GPU_RATE")")"
if [[ "$ROUNDS" -gt 0 ]]; then
  total_h=$(awk -v h="$ROUND_HOURS" -v n="$ROUNDS" 'BEGIN{printf "%.2f", h*n}')
  printf "%-42s | %-12s | %-18s | %-12s | %-18s\n" \
    "   -> $ROUNDS-round total (${total_h}h)" \
    "$(fmt_rate "$HETZNER_GPU_RATE")" "$(fmt_cost "$(cost_hetzner "$total_h" "$HETZNER_GPU_RATE")")" \
    "$(fmt_rate "$GPUAI_GPU_RATE")" "$(fmt_cost "$(cost_gpuai "$total_h" "$GPUAI_GPU_RATE")")"
fi
printf "%-42s | %-12s | %-18s | %-12s | %-18s\n" \
  "2. MM-PBSA, all 4 sys (CPU, ${MMPBSA_HOURS}h) [*]" \
  "$(fmt_rate "$HETZNER_CPU_RATE")" "$(fmt_cost "$(cost_hetzner "$MMPBSA_HOURS" "$HETZNER_CPU_RATE")")" \
  "$(fmt_rate "$GPUAI_CPU_RATE")" "$(fmt_cost "$(cost_gpuai "$MMPBSA_HOURS" "$GPUAI_CPU_RATE")")"
printf "%-42s | %-12s | %-18s | %-12s | %-18s\n" \
  "3. Trajectory analysis (CPU, ${ANALYSIS_HOURS}h) [*]" \
  "$(fmt_rate "$HETZNER_CPU_RATE")" "$(fmt_cost "$(cost_hetzner "$ANALYSIS_HOURS" "$HETZNER_CPU_RATE")")" \
  "$(fmt_rate "$GPUAI_CPU_RATE")" "$(fmt_cost "$(cost_gpuai "$ANALYSIS_HOURS" "$GPUAI_CPU_RATE")")"

cat <<'EOF'

[*] Steps 2 and 3 are CPU-bound and need no GPU -- they cost $0 on your own laptop. They're
    priced here only because you chose to keep them on rented infra for workflow consistency.

Notes on reading this table:
  - Row 1's hours (ROUND_HOURS) is a PLACEHOLDER until you have one real convergence-round
    wall-clock time from your own run; override with ROUND_HOURS=<hours> ./roots_cost_compare.sh.
    The number of rounds needed to reach convergence is NOT predictable up front (that's the
    whole point of "extend until it converges") -- pass ROUNDS=N only to sanity-check an
    assumed worst case, not as a forecast.
  - If every Hetzner cell shows "unknown", you likely don't have HCLOUD_TOKEN set. Get a
    read-only token from the Hetzner Cloud Console (your project -> Security -> API Tokens ->
    Generate API Token, Permissions: Read), then: export HCLOUD_TOKEN=... and re-run.
  - If GPU.ai's CPU row shows "unknown (no matching type id found...)", your account's catalog
    doesn't expose a type with "cpu" in its id/name the way this script guesses -- run
    `gpu gpu-types list --output json | jq .` yourself, find the right id, and pass it as
    GPUAI_CPU_TYPE=<id> ./roots_cost_compare.sh.
  - "unknown" anywhere means the live query failed or didn't match -- it is never silently
    replaced with a guessed number, and it never gets costed as $0.00 either (that was a bug
    in the first version of this script -- fixed: cost functions now check numeric-ness of the
    rate, not emptiness, before doing arithmetic on it).
EOF
