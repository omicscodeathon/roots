#!/usr/bin/env bash
# gpuai_resume.sh (v1 -- ported from civo_resume.sh v3) -- resume the four ROOTS production runs
# (Sys1-Sys4) on ONE GPU.ai GPU instance from the checkpoints in your downloaded results folders,
# and EXTEND each run until it ACTUALLY CONVERGES -- no more stopping at a flat 50 ns or 100 ns
# mark. Run this on YOUR LOCAL machine: it shells out to the `gpu` CLI (talks to api.gpu.ai with
# YOUR GPUAI_API_KEY) and to `ssh`/`scp` through the tunnel `gpu` resolves. Nothing runs against
# GPU.ai from anywhere but your machine.
#
# THIS IS A PORT, NOT A FIND/REPLACE. Three things about GPU.ai are genuinely different from
# Civo, and the provisioning/connection layer below was rewritten (not renamed) to match:
#
#   1. NO PUBLIC IP. GPU.ai instances sit behind a shared FRP tunnel gateway; there is no
#      "civo instance show | grep IP". `gpu instances ssh <id> --print` is the documented way
#      to resolve the real host:port:key for that instance. This script calls it ONCE per
#      command invocation (resolve_conn), parses the printed ssh command, and reuses the
#      parsed host/port/key for every ssh/scp call after that -- exactly like Civo's fixed
#      $CIVO_IP, just sourced from a different place because GPU.ai's networking model is
#      different.
#   2. PER-SECOND BILLING, NOT HOURLY-ROUNDED. Civo bills the hour, rounded up, so the old
#      spend-guard math did `hh=ceil(elapsed_hours); billed=hh*rate`. GPU.ai's CLI docs are
#      explicit: "per-second metering... no minimums, no hourly rounding." The spend-guard
#      and status math below were changed to `billed = elapsed_hours * rate`, no ceiling --
#      the Civo version would have overcharged every estimate here by up to an hour.
#   3. THE CLI RETURNS THE REAL PRICE. `gpu instances create` hands back the actual
#      price_per_hour that placement landed on (GPU.ai is a price-sorted marketplace across
#      many providers, not a fixed per-size rate card), so this script reads that back and
#      uses it everywhere, instead of asking you to hand-type a guessed HOURLY_RATE like the
#      Civo version did.
#
# Everything else -- the GROMACS resume/append logic, the RMSD convergence screen, the
# extend-until-converged loop, the budget-guard flag names -- is UNCHANGED from civo_resume.sh.
# That logic never talked to Civo; it only talked to GROMACS and the local filesystem.
#
# Subcommands:
#   types       free, read-only: live price + availability for A100 80GB, across GPU.ai's
#               whole marketplace (not per-region loops -- GPU.ai routes placement itself)
#   create      tries each candidate GPU type in GPUAI_TYPES, capped at MAX_PRICE_USD,
#               until one boots. Nothing bills until one succeeds.
#   preflight   free, local only: finds each system's checkpoint and says what can resume
#   launch      uploads, installs GROMACS (GPU), resumes all systems, attaches convergence loops
#   status      progress per system + this instance's actual elapsed/billed, from `gpu instances get`
#   balance     your GPU.ai credit balance, burn rate and runway (`gpu billing balance`) --
#               this is where your deposit shows up
#   watch       polls, desktop-notifies on any flag, exits when every system has finished
#   fetch       copies results back to $OUT_ROOT (never touches your original downloads)
#   delete      terminates the instance (`gpu instances delete`) -- GPU.ai bills until this runs
#
# CONVERGENCE (per system, evaluated at the end of every segment, first at the end of the
# existing 50 ns run) -- IDENTICAL to civo_resume.sh, not touched by this port:
#   - RMSD of the SOLU group (protein+ligand; Sys3/Sys4's index.ndx has no Backbone group).
#   - Looks at the last WIN_NS (10 ns) of trace against the 10 ns before it:
#       |slope of the last window| <= SLOPE_MAX      (default 0.01 nm/ns)
#       |mean(last window) - mean(previous)| <= MEAN_DELTA_MAX   (default 0.05 nm)
#     The mean-difference test catches slow drift that a slope threshold alone lets through.
#   - Converged -> stops, writes <Sys>_converged_frame_ns.txt = start (ns) of the stable window,
#     for MM-PBSA to use frames from there to the end.
#   - Not converged -> extends by EXTEND_NS (10) via gmx convert-tpr -extend and resumes from
#     the final checkpoint, repeating FOR AS LONG AS IT TAKES. MAX_NS defaults to effectively
#     no ceiling (set MAX_NS=50 to reproduce a plain 50 ns run with no extension).
#   - The real stop condition besides convergence is a SPENDING GUARD: before every extension,
#     each system checks elapsed_hours x ACTUAL_RATE (per-second, no rounding -- see point 2
#     above) against MAX_SPEND_USD. Over budget -> stops and flags
#     ALERT_<Sys>_BUDGET_GUARD_STOPPED instead of extending blindly.
#   These are SCREENING heuristics, not a validated criterion. Always look at the RMSD plots.
#
# PREREQUISITES (one-time):
#   curl -fsSL https://gpu.ai/install | sh                         # installs the `gpu` CLI
#   export GPUAI_API_KEY=gpuai_live_xxxxxxxxxxxxxxxxxxxxxxxx        # private beta: ask GPU.ai
#                                                                    # for a key; put this in
#                                                                    # your shell rc, not just
#                                                                    # this terminal
#   gpu keys create --name laptop --public-key-file ~/.ssh/id_ed25519.pub   # if you haven't
#   command -v jq >/dev/null || (brew install jq || sudo apt-get install -y jq)
#
# BEFORE `launch`, you need a running instance. If you don't have one:
#   ./gpuai_resume.sh types                 # see live price/availability for A100 80GB
#   ./gpuai_resume.sh create                # tries GPUAI_TYPES until one boots, capped at
#                                            # MAX_PRICE_USD (default 1.23, your locked-in rate)
# then, once you have an instance (create prints the id, also saved to ~/.gpuai_create_env):
#   GPUAI_ID=<id> ./gpuai_resume.sh launch
#
# Lessons carried over from civo_resume.sh: state values written with printf %q; resume uses
# the same output names (-deffnm) and appends, falling back to -noappend if the .xtc/.edr/.log
# are missing; -cpt 15; .tpr/.cpt hashes verified after upload; jobs spread across however many
# GPUs the instance has (system i -> GPU i mod #GPUs), threads per job = min(vCPUs/#jobs, MAX_THREADS=8).
set -euo pipefail

LOCAL_ROOT="${LOCAL_ROOT:-$HOME/roots_pipeline1}"
RESULT_ROOTS=("$LOCAL_ROOT/Sys1+Sys2_results" "$LOCAL_ROOT/Sys3+Sys4_results")
SYSTEMS=(Sys1 Sys2 Sys3 Sys4)
OUT_ROOT="${OUT_ROOT:-$LOCAL_ROOT/GPUai_resumed}"
STATE_FILE="${STATE_FILE:-$HOME/.gpuai_md_state}"
SSH_KEY="${SSH_KEY:-$HOME/.ssh/id_ed25519}"          # fallback identity if --print doesn't name one
MAX_SPEND_USD="${MAX_SPEND_USD:-50}"                 # spending guard: convergence loops stop extending past this
POLL_SEC="${POLL_SEC:-300}"
MAX_THREADS="${MAX_THREADS:-8}"                      # cap OpenMP threads per job
MIN_NS="${MIN_NS:-50}"; MAX_NS="${MAX_NS:-100000}"; EXTEND_NS="${EXTEND_NS:-10}"
WIN_NS="${WIN_NS:-10}"; SLOPE_MAX="${SLOPE_MAX:-0.01}"; MEAN_DELTA_MAX="${MEAN_DELTA_MAX:-0.05}"
SSH_OPTS=(-o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null -o ServerAliveInterval=30 \
          -o ServerAliveCountMax=6 -o ConnectTimeout=20)
# ^ StrictHostKeyChecking=no / UserKnownHostsFile=/dev/null is intentional, not a weakening:
#   GPU.ai's own `gpu instances ssh` does the same thing, because the FRP gateway is a shared
#   host whose key isn't meaningfully "the instance's" key the way a Civo public IP's was.
LAUNCH_EPOCH=""

# --- only used by `types` / `create` (instance provisioning), not by the MD side ---
GPUAI_NAME="${GPUAI_NAME:-roots-md}"
GPUAI_ENVIRONMENT="${GPUAI_ENVIRONMENT:-raw-vm}"     # root VM access, like Civo's cuda12 image; NOT "certified" (may be container-native / no root)
GPUAI_SSH_KEY_ID="${GPUAI_SSH_KEY_ID:-}"             # set this, or let ensure_ssh_key() find/create one
GPUAI_KEY_NAME="${GPUAI_KEY_NAME:-roots-laptop}"     # display name used if we have to `gpu keys create` for you
GPUAI_REGION="${GPUAI_REGION:-}"                     # empty = GPU.ai places anywhere (it already price-sorts across providers)
# Candidate GPU type IDs for "A100 80GB", tried in order. CONFIRMED live against this
# account's own `gpu gpu-types list`: the catalog has no pcie/sxm split at all -- it's just
# VRAM-based: a100_80gb and a100_40gb. (The docs' own examples used "a100_pcie", which turned
# out not to exist on this catalog -- `gpu pricing --gpu-type a100_pcie` returns a header and
# zero rows. Fixed here; a100_80gb is the real id.)
GPUAI_TYPES="${GPUAI_TYPES:-a100_80gb}"
MAX_PRICE_USD="${MAX_PRICE_USD:-1.23}"               # hard price cap on create -- your locked-in rate
RETRY_PASS_SLEEP="${RETRY_PASS_SLEEP:-120}"          # seconds between full passes over GPUAI_TYPES
MAX_WAIT_MIN="${MAX_WAIT_MIN:-0}"                    # 0 = keep retrying until it works (Ctrl+C to give up); else give up after N minutes

declare -A SYS_DIR SYS_FNM SYS_NDX SYS_MODE SYS_NOTE SYS_NS
READY=()
CONN_HOST=""; CONN_PORT=""; CONN_KEY=""; CONN_USER=""

need_bin() { command -v "$1" >/dev/null || { echo "'$1' not found in PATH. $2" >&2; exit 1; }; }

ensure_cli() {
  need_bin gpu "Install it first: curl -fsSL https://gpu.ai/install | sh"
  need_bin jq "Install it first (brew install jq / sudo apt-get install -y jq) -- this script parses the gpu CLI's JSON output."
  [[ -n "${GPUAI_API_KEY:-}" ]] || { echo "GPUAI_API_KEY is not set. export GPUAI_API_KEY=gpuai_live_... (ask GPU.ai for a key during private beta)." >&2; exit 1; }
}

ensure_ssh_key() {
  [[ -n "$GPUAI_SSH_KEY_ID" ]] && return 0
  local existing
  existing=$(gpu keys list --output json 2>/dev/null | jq -r --arg n "$GPUAI_KEY_NAME" '.[]? | select(.name==$n) | .id' | head -1 || true)
  if [[ -n "$existing" ]]; then GPUAI_SSH_KEY_ID="$existing"; return 0; fi
  [[ -f "${SSH_KEY}.pub" ]] || { echo "No SSH key id set (GPUAI_SSH_KEY_ID) and ${SSH_KEY}.pub not found to upload. Run: ssh-keygen -t ed25519 -f $SSH_KEY" >&2; exit 1; }
  echo "No '$GPUAI_KEY_NAME' key on GPU.ai yet -- uploading ${SSH_KEY}.pub ..."
  GPUAI_SSH_KEY_ID=$(gpu keys create --name "$GPUAI_KEY_NAME" --public-key-file "${SSH_KEY}.pub" --output json | jq -r '.id')
  [[ -n "$GPUAI_SSH_KEY_ID" && "$GPUAI_SSH_KEY_ID" != "null" ]] || { echo "Failed to register SSH key with GPU.ai." >&2; exit 1; }
  echo "Registered as key id $GPUAI_SSH_KEY_ID"
}

# Resolves CONN_HOST/CONN_PORT/CONN_KEY/CONN_USER for $GPUAI_ID by parsing
# `gpu instances ssh --print`, the documented way to get the FRP tunnel's real host:port
# (GPU.ai instances have no public IP of their own -- see header comment, point 1).
resolve_conn() {
  [[ -n "${GPUAI_ID:-}" ]] || { echo "GPUAI_ID is not set." >&2; return 1; }
  local out
  out=$(gpu instances ssh "$GPUAI_ID" --print --wait=true 2>&1) || { echo "could not resolve SSH tunnel for $GPUAI_ID:" >&2; echo "$out" >&2; return 1; }
  CONN_PORT=$(grep -oE -- '-p[[:space:]]+[0-9]+' <<<"$out" | grep -oE '[0-9]+' | head -1)
  CONN_KEY=$(grep -oE -- '-i[[:space:]]+[^ ]+' <<<"$out" | awk '{print $2}' | head -1)
  local userhost; userhost=$(grep -oE '[A-Za-z0-9_.-]+@[A-Za-z0-9_.-]+' <<<"$out" | head -1)
  CONN_USER="${userhost%@*}"; CONN_HOST="${userhost#*@}"
  CONN_KEY="${CONN_KEY:-$SSH_KEY}"
  [[ -n "$CONN_PORT" && -n "$userhost" ]] || { echo "unexpected 'gpu instances ssh --print' output, could not parse host/port:" >&2; echo "$out" >&2; return 1; }
}

ssh_run() { ssh "${SSH_OPTS[@]}" -p "$CONN_PORT" -i "$CONN_KEY" "${CONN_USER}@${CONN_HOST}" "$@"; }
scp_up()  { scp "${SSH_OPTS[@]}" -P "$CONN_PORT" -i "$CONN_KEY" -q "$1" "${CONN_USER}@${CONN_HOST}:$2"; }
tar_up() {
  local src="$1"; shift
  local dest="${@: -1}"
  local files=("${@:1:$(($#-1))}")
  tar -C "$src" -cf - "${files[@]}" \
    | ssh "${SSH_OPTS[@]}" -p "$CONN_PORT" -i "$CONN_KEY" "${CONN_USER}@${CONN_HOST}" "mkdir -p $dest && tar -xf - -C $dest"
}

last_ns() {
  awk '/^ *Step *Time$/{getline; t=$2; f=1} END{if(f) printf "%.2f", t/1000; else printf "0.00"}' "$1" 2>/dev/null || printf "0.00"
}

cmd_types() {
  ensure_cli
  echo "Live A100 80GB price/availability on GPU.ai (read-only, no cost -- GPU.ai aggregates every provider, you don't pick a region loop like Civo):"
  local t
  for t in $GPUAI_TYPES; do
    echo; echo "--- $t ---"
    gpu pricing --gpu-type "$t" 2>&1 || true
  done
  echo
  echo "If none of these look like 'A100 80GB', run 'gpu gpu-types list' to find the current catalog ID and set GPUAI_TYPES."
}

cmd_create() {
  ensure_cli
  ensure_ssh_key
  local start_epoch; start_epoch=$(date +%s)
  local t out rc id price pass=1
  while :; do
    echo "=== create: pass $pass -- types [$GPUAI_TYPES], capped at \$${MAX_PRICE_USD}/hr ==="
    for t in $GPUAI_TYPES; do
      printf '[%s] trying %-16s ... ' "$(date +%H:%M:%S)" "$t"
      local errfile; errfile=$(mktemp)
      set +e
      # IMPORTANT: capture stdout and stderr SEPARATELY. The CLI writes human-readable
      # progress ("Submitted: operation ... polling...") to stderr and the structured
      # result to stdout; merging them with 2>&1 (the old bug) fed jq a blob of prose
      # plus JSON and silently lost real, already-billing instance ids on success.
      out=$(gpu instances create -t "$t" --tier on_demand --environment "$GPUAI_ENVIRONMENT" \
            --ssh-key-id "$GPUAI_SSH_KEY_ID" --max-price "$MAX_PRICE_USD" --name "$GPUAI_NAME" \
            ${GPUAI_REGION:+--region "$GPUAI_REGION"} --output json 2>"$errfile")
      rc=$?
      local err; err=$(cat "$errfile"); rm -f "$errfile"
      set -e
      if [[ $rc -eq 0 ]]; then
        id=$(jq -r '.id' <<<"$out" 2>/dev/null || true)
        price=$(jq -r '.price_per_hour' <<<"$out" 2>/dev/null || true)
        if [[ -z "$id" || "$id" == "null" ]]; then
          echo "FAILED (couldn't parse instance id from output -- but exit code was 0, so check"
          echo "'gpu instances list' before retrying, an instance may still have been created):"
          echo "$out" | tail -5; echo "$err" | tail -5; continue
        fi
        echo "OK"
        {
          printf 'GPUAI_ID=%q\n' "$id"
          printf 'GPUAI_TYPE=%q\n' "$t"
          printf 'GPUAI_RATE=%q\n' "${price:-$MAX_PRICE_USD}"
        } > "$HOME/.gpuai_create_env"
        echo
        echo "Instance is up: $id ($t), actual rate \$${price:-?}/hr (capped at \$${MAX_PRICE_USD})."
        echo "Saved to $HOME/.gpuai_create_env. Next:"
        echo "  GPUAI_ID=$id ./gpuai_resume.sh launch"
        return 0
      fi
      if grep -qiE '422|no capacity|unavailable|price_changed|insufficient|price (for this offering )?moved|above the \$|approved at' <<<"$err$out"; then
        echo "no stock / over price cap"
        grep -iE 'moved to \$[0-9.]+' <<<"$err" | head -1
      else
        echo "FAILED (not a stock/price issue)"
        echo "$err" | tail -5
        echo "    -- skipping, but that error is worth reading: likely env/ssh-key/type-id, not availability."
      fi
    done
    if [[ "$MAX_WAIT_MIN" -gt 0 ]]; then
      local elapsed_min=$(( ($(date +%s) - start_epoch) / 60 ))
      if [[ "$elapsed_min" -ge "$MAX_WAIT_MIN" ]]; then
        echo "Gave up after ${elapsed_min} min (MAX_WAIT_MIN=$MAX_WAIT_MIN). Nothing was left running/billing."
        return 1
      fi
    fi
    echo "Nothing free this pass at \$${MAX_PRICE_USD}/hr or under. Sleeping ${RETRY_PASS_SLEEP}s before pass $((pass+1)). Ctrl+C to stop -- nothing is billing yet."
    sleep "$RETRY_PASS_SLEEP"
    pass=$((pass+1))
  done
}

discover() {
  local sys ov cpt root dir fnm ndx
  READY=()
  for sys in "${SYSTEMS[@]}"; do
    ov="${sys^^}_DIR"
    cpt=""
    if [[ -n "${!ov:-}" ]]; then
      cpt=$(find "${!ov}" -maxdepth 4 -type f -name 'prod*.cpt' ! -name '*_prev.cpt' 2>/dev/null | head -1 || true)
    else
      for root in "${RESULT_ROOTS[@]}"; do
        [[ -d "$root" ]] || continue
        cpt=$(find "$root" -type f -path "*/${sys}/*" -name 'prod*.cpt' ! -name '*_prev.cpt' 2>/dev/null | head -1 || true)
        if [[ -n "$cpt" ]]; then break; fi
      done
    fi
    if [[ -z "$cpt" ]]; then
      SYS_NOTE[$sys]="BLOCKED: no prod*.cpt checkpoint found (cannot resume continuously)"
      continue
    fi
    dir=$(dirname "$cpt"); fnm=$(basename "$cpt" .cpt)
    ndx=""
    if [[ -f "$dir/index.ndx" ]]; then ndx="$dir/index.ndx"
    else ndx=$(find "$dir/.." -maxdepth 2 -type f -name index.ndx 2>/dev/null | head -1 || true); fi
    SYS_DIR[$sys]="$dir"; SYS_FNM[$sys]="$fnm"; SYS_NDX[$sys]="$ndx"
    if [[ -f "$dir/$fnm.log" ]]; then SYS_NS[$sys]=$(last_ns "$dir/$fnm.log"); else SYS_NS[$sys]="?"; fi
    if [[ ! -f "$dir/$fnm.tpr" ]]; then
      SYS_NOTE[$sys]="BLOCKED: $fnm.tpr missing next to the checkpoint"; continue
    fi
    if [[ -f "$dir/$fnm.xtc" && -f "$dir/$fnm.edr" && -f "$dir/$fnm.log" ]]; then
      SYS_MODE[$sys]="append"; SYS_NOTE[$sys]="OK, will append to existing outputs"
    else
      SYS_MODE[$sys]="noappend"
      SYS_NOTE[$sys]="OK, but .xtc/.edr/.log incomplete -> -noappend (new part files, NO auto-extension); progress unknown if log missing (counted as ${MIN_NS} ns remaining)"
    fi
    [[ -n "$ndx" ]] || SYS_NOTE[$sys]+="; no index.ndx (RMSD check will fail)"
    READY+=("$sys")
  done
}

cmd_preflight() {
  discover
  printf '\n%-6s %-12s %9s  %s\n' System deffnm "ns so far" Status
  local sys total_rem=0 n=0 c
  for sys in "${SYSTEMS[@]}"; do
    printf '%-6s %-12s %9s  %s\n' "$sys" "${SYS_FNM[$sys]:--}" "${SYS_NS[$sys]:--}" "${SYS_NOTE[$sys]}"
    if [[ -n "${SYS_MODE[$sys]:-}" ]]; then
      n=$((n+1)); c="${SYS_NS[$sys]//[!0-9.]/}"
      total_rem=$(awk -v a="$total_rem" -v t="$MIN_NS" -v c="${c:-0}" 'BEGIN{printf "%.2f", a + (t>c ? t-c : 0)}')
    fi
  done
  echo
  echo "Resumable systems: $n of ${#SYSTEMS[@]}   |   remaining to the first convergence check (${MIN_NS} ns): ${total_rem} ns"
  if [[ "$n" -gt 0 ]]; then
    awk -v r="$total_rem" -v rate="$MAX_PRICE_USD" -v mn="$MIN_NS" -v cap="$MAX_SPEND_USD" -v ext="$EXTEND_NS" 'BEGIN{
      lo=r/620*24; hi=r/480*24;
      printf "Wall-clock to reach the first convergence check (%s ns each): %.1f-%.1f h (+~0.5 h setup) on an A100-class GPU.\n", mn, lo, hi;
      printf "At $%.2f/h, billed per-second (no hourly rounding): about $%.0f-$%.0f just to get there.\n", rate, lo*rate, hi*rate;
      print "(Assumes 480-620 ns/day aggregate; NOT measured for 4 concurrent jobs on an A100. Adjust for your actual card.)";
      printf "Past that: each system extends %s ns at a time until its RMSD screen converges -- no fixed 100 ns stop anymore.\n", ext;
      printf "Spending guard: MAX_SPEND_USD=$%.0f -- once estimated billing hits that, a system stops extending (flags BUDGET_GUARD) even if not converged.\n", cap }'
  fi
  echo
  echo "Nothing was uploaded or started. If a system says BLOCKED, paste this table and we decide how to handle it."
}

save_state() {
  {
    printf 'GPUAI_ID=%q\n' "$GPUAI_ID"
    printf 'GPUAI_RATE=%q\n' "$GPUAI_RATE"
    printf 'LAUNCH_EPOCH=%q\n' "$LAUNCH_EPOCH"
    printf 'MAX_SPEND_USD=%q\n' "$MAX_SPEND_USD"
    printf 'READY_LIST=%q\n' "${READY[*]}"
    local s
    for s in "${READY[@]}"; do printf 'FNM_%s=%q\n' "$s" "${SYS_FNM[$s]}"; done
  } > "$STATE_FILE"
  chmod 600 "$STATE_FILE"
}

load_state() {
  [[ -f "$STATE_FILE" ]] || { echo "No state file at $STATE_FILE -- run 'launch' first." >&2; exit 1; }
  # shellcheck disable=SC1090
  source "$STATE_FILE"
  read -ra READY <<< "$READY_LIST"
}

cmd_launch() {
  ensure_cli
  [[ -n "${GPUAI_ID:-}" ]] || { echo "Set GPUAI_ID=<instance id> (from 'create', or 'gpu instances list')." >&2; exit 1; }
  LAUNCH_EPOCH=$(date +%s)
  discover
  [[ ${#READY[@]} -gt 0 ]] || { echo "No resumable systems. Run './gpuai_resume.sh preflight'." >&2; exit 1; }

  echo "=== 0. Resolve instance + confirm actual billed rate ==="
  local info; info=$(gpu instances get "$GPUAI_ID" --output json)
  GPUAI_RATE=$(jq -r '.price_per_hour // empty' <<<"$info")
  [[ -n "$GPUAI_RATE" ]] || { echo "Couldn't read price_per_hour from 'gpu instances get $GPUAI_ID' -- is the instance id right?" >&2; echo "$info" >&2; exit 1; }
  echo "Instance $GPUAI_ID is billing at \$${GPUAI_RATE}/hr (per-second, no rounding)."
  resolve_conn || exit 1
  echo "Tunnel resolved: ${CONN_USER}@${CONN_HOST}:${CONN_PORT}"

  echo; echo "=== 1. Reachability + GPU check ==="
  ssh_run 'echo "login ok as $(whoami)"; echo "vCPUs: $(nproc)"; command -v nvidia-smi >/dev/null && nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv || echo "NO nvidia-smi"' \
    || { echo "SSH failed through the tunnel. Try: gpu instances ssh $GPUAI_ID  (interactively) to see the raw error." >&2; exit 1; }
  if ! ssh_run 'command -v nvidia-smi >/dev/null'; then
    echo "This instance has no NVIDIA driver visible. STOP: wrong type/environment. Delete it ('./gpuai_resume.sh delete') and re-create." >&2
    exit 1
  fi
  read -rp "Does the GPU line above show the card you intended (A100 80GB)? Type YES to continue: " OKGPU
  [[ "$OKGPU" == "YES" ]] || { echo "Aborted. Delete the instance if it is wrong ('./gpuai_resume.sh delete') -- it bills until deleted."; exit 1; }

  echo; echo "=== 2. Upload checkpoints + outputs (tar over ssh, through the FRP tunnel) ==="
  local sys dir fnm files f
  for sys in "${READY[@]}"; do
    dir="${SYS_DIR[$sys]}"; fnm="${SYS_FNM[$sys]}"; files=()
    for f in "$fnm.tpr" "$fnm.cpt" "${fnm}_prev.cpt" "$fnm.xtc" "$fnm.edr" "$fnm.log"; do
      [[ -f "$dir/$f" ]] && files+=("$f")
    done
    echo "--- $sys: ${files[*]}  ($(du -ch "${files[@]/#/$dir/}" 2>/dev/null | tail -1 | cut -f1))"
    ssh_run "mkdir -p ~/md/$sys"
    tar_up "$dir" "${files[@]}" "~/md/$sys"
    [[ -n "${SYS_NDX[$sys]}" ]] && scp_up "${SYS_NDX[$sys]}" "md/$sys/index.ndx"
    for f in "$fnm.tpr" "$fnm.cpt"; do
      local a b
      a=$(sha256sum "$dir/$f" | cut -d' ' -f1)
      b=$(ssh_run "sha256sum ~/md/$sys/$f" | cut -d' ' -f1)
      [[ "$a" == "$b" ]] || { echo "HASH MISMATCH on $sys/$f -- upload corrupted. Re-run launch." >&2; exit 1; }
    done
    echo "    hashes match for $fnm.tpr and $fnm.cpt"
  done

  echo; echo "=== 3. Install GPU GROMACS on the instance (skipped if already present) ==="
  ssh_run bash -s <<'REMOTE_INSTALL'
set -eo pipefail
if [ ! -d "$HOME/miniforge3" ]; then
  URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
  if command -v curl >/dev/null; then curl -fsSL "$URL" -o /tmp/mf.sh; else wget -q "$URL" -O /tmp/mf.sh; fi
  bash /tmp/mf.sh -b -p "$HOME/miniforge3"
fi
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda env list | grep -q gmx_gpu || conda create -y -n gmx_gpu -c conda-forge "gromacs=*=*cuda*"
conda activate gmx_gpu
GMXRC=$(find "$HOME/miniforge3/envs/gmx_gpu" -name GMXRC.bash | head -1)
source "$GMXRC"
gmx --version 2>&1 | grep -E 'GROMACS version|GPU support|CUDA driver|CUDA runtime'
REMOTE_INSTALL

  echo; echo "=== 4. Launch resumed production for all systems ==="
  local TMP; TMP=$(mktemp -d)
  cat > "$TMP/launch_all.sh" <<'REMOTE_LAUNCH'
#!/usr/bin/env bash
# args: SYS:FNM:MODE ...   Writes run.env per system so the convergence loop can relaunch identically.
set -eo pipefail
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate gmx_gpu
GMXRC=$(find "$HOME/miniforge3/envs/gmx_gpu" -name GMXRC.bash | head -1)
source "$GMXRC"
N=$#; NPROC=$(nproc); NT=$(( NPROC / N )); [ "$NT" -lt 1 ] && NT=1
[ "$NT" -gt "${MAX_THREADS:-8}" ] && NT="${MAX_THREADS:-8}"
NGPU=$(nvidia-smi -L 2>/dev/null | wc -l); [ "$NGPU" -lt 1 ] && NGPU=1
OFF=0; I=0
echo "cores=$NPROC gpus=$NGPU jobs=$N threads/job=$NT"
for entry in "$@"; do
  IFS=: read -r SYS FNM MODE <<< "$entry"
  cd "$HOME/md/$SYS"
  GPU=$(( I % NGPU ))
  printf 'FNM=%s\nMODE=%s\nNT=%s\nOFF=%s\nGPU=%s\n' "$FNM" "$MODE" "$NT" "$OFF" "$GPU" > run.env
  APP=""; [ "$MODE" = "noappend" ] && APP="-noappend"
  nohup gmx mdrun -deffnm "$FNM" -s "$FNM.tpr" -cpi "$FNM.cpt" -cpt 15 $APP \
      -nb gpu -pme gpu -bonded gpu -update gpu \
      -ntmpi 1 -ntomp "$NT" -pin on -pinoffset "$OFF" -pinstride 1 -gpu_id "$GPU" \
      > resume_seg1.stdout 2>&1 &
  echo $! > mdrun.pid
  disown
  echo "$SYS launched, pid $(cat mdrun.pid), gpu $GPU, pin offset $OFF"
  OFF=$(( OFF + NT )); I=$(( I + 1 ))
done
sleep 75
for entry in "$@"; do
  IFS=: read -r SYS FNM MODE <<< "$entry"
  cd "$HOME/md/$SYS"
  if kill -0 "$(cat mdrun.pid)" 2>/dev/null; then echo "[$SYS] ALIVE"; else echo "[$SYS] DEAD -- see resume_seg1.stdout"; fi
  grep -iE 'checkpoint|fatal|error' resume_seg1.stdout | head -4 || true
  tail -n 2 "$FNM.log" || true
done
REMOTE_LAUNCH

  cat > "$TMP/converge.sh" <<'REMOTE_CONVERGE'
#!/usr/bin/env bash
# usage: converge.sh WORKDIR SYS FNM
# Waits for each production segment, screens RMSD convergence, extends or stops.
WD="$1"; SYS="$2"; FNM="$3"
MIN_NS="${MIN_NS:-50}"; MAX_NS="${MAX_NS:-100000}"; EXTEND_NS="${EXTEND_NS:-10}"
WIN_NS="${WIN_NS:-10}"; SLOPE_MAX="${SLOPE_MAX:-0.01}"; MEAN_DELTA_MAX="${MEAN_DELTA_MAX:-0.05}"
EARLY_NS="${EARLY_NS:-25}"; GROUP="${RMSD_GROUP:-SOLU}"; INTERVAL="${INTERVAL_SEC:-300}"
cd "$WD" || exit 1
[ -f run.env ] && . ./run.env          # FNM MODE NT OFF GPU written by the launcher
if [ -z "${SKIP_ENV:-}" ]; then
  source "$HOME/miniforge3/etc/profile.d/conda.sh"; conda activate gmx_gpu
  GMXRC=$(find "$HOME/miniforge3/envs/gmx_gpu" -name GMXRC.bash | head -1); source "$GMXRC"
fi
log() { echo "[$(date -u +%FT%TZ)] $*" >> watchdog.log; }
ge()  { awk -v a="$1" -v b="$2" 'BEGIN{exit !(a>=b)}'; }
cur_ns() { awk '/^ *Step *Time$/{getline; t=$2; f=1} END{if(f) printf "%.3f", t/1000; else print 0}' "$FNM.log" 2>/dev/null; }
seg_end_ns() {   # end time of the current .tpr in ns, from nsteps*delta_t; fall back to MIN_NS
  gmx dump -s "$FNM.tpr" 2>/dev/null | awk '/^ *nsteps *=/{n=$3} /^ *delta_t *=/{d=$3} END{if(n>0&&d>0) printf "%.3f", n*d/1000; else printf "%s", ENVIRON["MIN_NS"]}'
}
rms_trace() {    # $1 = tag; prints the xvg path (empty if it failed)
  printf '%s\n%s\n' "$GROUP" "$GROUP" | gmx rms -s "$FNM.tpr" -f "$FNM.xtc" -o "rmsd_$1.xvg" -tu ns -n index.ndx >/dev/null 2>"rmsd_$1.err" || true
  [ -s "rmsd_$1.xvg" ] && echo "rmsd_$1.xvg" || echo ""
}
evaluate() {     # $1 = xvg ; prints: T slope delta VERDICT
  [ -n "$1" ] || { echo "NA NA NA FAILED"; return; }
  grep -v '^[#@]' "$1" | awk -v W="$WIN_NS" -v SM="$SLOPE_MAX" -v DM="$MEAN_DELTA_MAX" -v MIN="$MIN_NS" '
    {x[NR]=$1; y[NR]=$2; n=NR}
    END{
      if(n<4){print "NA NA NA FAILED"; exit}
      T=x[n]; s1=T-W; s0=T-2*W
      if(s0 < x[1]-1e-9){printf "%.2f NA NA INSUFFICIENT\n", T; exit}
      c=0;sx=0;sy=0;sxy=0;sxx=0;m1=0;k1=0;m0=0;k0=0
      for(i=1;i<=n;i++){
        if(x[i]>=s1-1e-9){c++;sx+=x[i];sy+=y[i];sxy+=x[i]*y[i];sxx+=x[i]*x[i];m1+=y[i];k1++}
        else if(x[i]>=s0-1e-9){m0+=y[i];k0++}
      }
      d=c*sxx-sx*sx
      if(c<3||k0<3||d==0){print "NA NA NA FAILED"; exit}
      sl=(c*sxy-sx*sy)/d; if(sl<0)sl=-sl
      dm=(m1/k1)-(m0/k0); if(dm<0)dm=-dm
      v=(sl<=SM && dm<=DM && T>=MIN-1e-9) ? "CONVERGED" : "NOT_CONVERGED"
      printf "%.2f %.5f %.4f %s\n", T, sl, dm, v }'
}
finish() { touch "DONE_${SYS}.flag"; log "DONE: $1"; }

SEG_END=$(seg_end_ns); SEG=1; EARLY_DONE=0
log "convergence loop start: sys=$SYS fnm=$FNM segment_end=${SEG_END}ns min=${MIN_NS} max=${MAX_NS} extend=${EXTEND_NS} window=${WIN_NS} slope_max=${SLOPE_MAX} mean_delta_max=${MEAN_DELTA_MAX} group=${GROUP} mode=${MODE:-append} spend_cap=\$${MAX_SPEND_USD:-none} rate=\$${HOURLY_RATE:-?}/hr(per-second)"
while :; do
  PID=$(cat mdrun.pid 2>/dev/null || echo 0)
  while kill -0 "$PID" 2>/dev/null; do
    sleep "$INTERVAL"; NS=$(cur_ns); log "progress ${NS} ns (segment ${SEG}, ends at ${SEG_END})"
    if [ $EARLY_DONE -eq 0 ] && ge "$NS" "$EARLY_NS"; then
      EARLY_DONE=1; R=$(evaluate "$(rms_trace early)")
      log "EARLY LOOK (informational only, never stops or extends the run) at ~${NS} ns: ${R}"
    fi
  done
  NS=$(cur_ns); NEED=$(awk -v a="$SEG_END" 'BEGIN{printf "%.3f", a-0.05}')
  if ! ge "$NS" "$NEED"; then
    log "mdrun EXITED EARLY at ${NS} ns (segment ${SEG} should end at ${SEG_END}) -- see resume_seg*.stdout"
    touch "STOPPED_EARLY_${SYS}.flag"; exit 1
  fi
  R=$(evaluate "$(rms_trace "seg${SEG}")")
  echo "segment=${SEG} ${R}" >> convergence_summary.txt
  read -r T SL DM V <<< "$R"
  log "Segment ${SEG} ended at ${NS} ns. Convergence screen: T=${T} slope=${SL} nm/ns mean_delta=${DM} nm -> ${V}"
  if [ "$V" = "CONVERGED" ]; then
    START=$(awk -v t="$T" -v w="$WIN_NS" 'BEGIN{printf "%.1f", t-w}')
    echo "$START" > "${SYS}_converged_frame_ns.txt"; touch "CONVERGED_${SYS}.flag"
    log "CONVERGED at ${T} ns; stable window starts at ${START} ns (wrote ${SYS}_converged_frame_ns.txt)"
    finish "converged"; exit 0
  fi
  if [ "$V" = "FAILED" ]; then
    log "RMSD screen FAILED (see rmsd_seg${SEG}.err) -- cannot judge convergence, NOT extending"
    touch "ALERT_${SYS}_RMSD_FAILED.flag"; finish "rmsd failed"; exit 1
  fi
  if [ -n "${MAX_SPEND_USD:-}" ] && [ -n "${LAUNCH_EPOCH:-}" ] && [ -n "${HOURLY_RATE:-}" ]; then
    NOWSEC=$(date +%s)
    ELAPSED_H=$(awk -v n="$NOWSEC" -v s="$LAUNCH_EPOCH" 'BEGIN{h=(n-s)/3600; if(h<0)h=0; printf "%.4f", h}')
    # GPU.ai bills per second, no hourly rounding -- straight elapsed_hours x rate (the Civo
    # version did hh=ceil(elapsed_hours); this is the fix described in the script header).
    BILLED=$(awk -v h="$ELAPSED_H" -v r="$HOURLY_RATE" 'BEGIN{printf "%.2f", h*r}')
    if awk -v b="$BILLED" -v cap="$MAX_SPEND_USD" 'BEGIN{exit !(b>=cap)}'; then
      log "SPENDING GUARD: ~${ELAPSED_H}h elapsed, ~\$${BILLED} billed (per-second) >= cap \$${MAX_SPEND_USD}. Not extending further."
      touch "ALERT_${SYS}_BUDGET_GUARD_STOPPED.flag"
      finish "budget guard stopped at ${T} ns, not converged"; exit 0
    fi
  fi
  if [ "${MODE:-append}" = "noappend" ]; then
    log "noappend mode: output is split into part files, automatic extension disabled. Not converged by this screen."
    touch "ALERT_${SYS}_NOT_CONVERGED.flag"; finish "noappend, not extended"; exit 0
  fi
  if ge "$SEG_END" "$MAX_NS"; then
    log "!!! NOT CONVERGED at the ${MAX_NS} ns ceiling (screen: ${V})"
    touch "ALERT_${SYS}_NOT_CONVERGED_AT_CEILING.flag"; finish "ceiling reached"; exit 0
  fi
  EXT=$(awk -v e="$EXTEND_NS" -v s="$SEG_END" -v m="$MAX_NS" 'BEGIN{r=m-s; printf "%d", (e<r?e:r)}')
  log "Not converged -> extending by ${EXT} ns (segment ${SEG} -> $((SEG+1)), new end $(awk -v s="$SEG_END" -v e="$EXT" 'BEGIN{print s+e}') ns)"
  cp "$FNM.tpr" "${FNM}_seg${SEG}.tpr"
  if ! gmx convert-tpr -s "$FNM.tpr" -extend "$((EXT*1000))" -o "${FNM}_ext.tpr" > "extend_seg${SEG}.log" 2>&1; then
    log "convert-tpr FAILED -- see extend_seg${SEG}.log"; touch "STOPPED_EARLY_${SYS}.flag"; exit 1
  fi
  mv "${FNM}_ext.tpr" "$FNM.tpr"
  SEG=$((SEG+1)); SEG_END=$(awk -v s="$SEG_END" -v e="$EXT" 'BEGIN{print s+e}')
  nohup gmx mdrun -deffnm "$FNM" -s "$FNM.tpr" -cpi "$FNM.cpt" -cpt 15 \
      -nb gpu -pme gpu -bonded gpu -update gpu \
      -ntmpi 1 -ntomp "${NT:-4}" -pin on -pinoffset "${OFF:-0}" -pinstride 1 -gpu_id "${GPU:-0}" \
      > "resume_seg${SEG}.stdout" 2>&1 &
  echo $! > mdrun.pid; disown
  sleep "${RESTART_CHECK_SEC:-20}"
  if ! kill -0 "$(cat mdrun.pid)" 2>/dev/null && ! ge "$(cur_ns)" "$(awk -v a="$SEG_END" 'BEGIN{printf "%.3f", a-0.05}')"; then
    log "extended run failed to start -- see resume_seg${SEG}.stdout"; touch "STOPPED_EARLY_${SYS}.flag"; exit 1
  fi
done
REMOTE_CONVERGE

  scp "${SSH_OPTS[@]}" -P "$CONN_PORT" -i "$CONN_KEY" -q "$TMP/launch_all.sh" "$TMP/converge.sh" "${CONN_USER}@${CONN_HOST}:md/"
  local entries=()
  for sys in "${READY[@]}"; do entries+=("${sys}:${SYS_FNM[$sys]}:${SYS_MODE[$sys]}"); done
  ssh_run "MAX_THREADS=$MAX_THREADS bash ~/md/launch_all.sh ${entries[*]}"

  echo; echo "=== 5. Attach convergence loops (screen + automatic extension; first check at ${MIN_NS} ns, extends ${EXTEND_NS} ns at a time until converged or spend cap \$${MAX_SPEND_USD} is hit) ==="
  local envp="MIN_NS=$MIN_NS MAX_NS=$MAX_NS EXTEND_NS=$EXTEND_NS WIN_NS=$WIN_NS SLOPE_MAX=$SLOPE_MAX MEAN_DELTA_MAX=$MEAN_DELTA_MAX MAX_SPEND_USD=$MAX_SPEND_USD HOURLY_RATE=$GPUAI_RATE LAUNCH_EPOCH=$LAUNCH_EPOCH"
  for sys in "${READY[@]}"; do
    ssh_run "cd ~/md/$sys && ( $envp nohup bash ~/md/converge.sh \$HOME/md/$sys $sys ${SYS_FNM[$sys]} >/dev/null 2>&1 & ); sleep 2; if pgrep -f 'converge.sh .*$sys ' >/dev/null; then echo '$sys convergence loop ALIVE'; else echo '$sys convergence loop NOT RUNNING'; fi"
  done
  save_state
  rm -rf "$TMP"
  echo
  echo "Launched ${READY[*]}. The instance bills per-second (no rounding) until you DELETE it ('./gpuai_resume.sh delete')."
  echo "Each system extends ${EXTEND_NS} ns at a time past ${MIN_NS} ns until its RMSD screen converges -- no fixed ns ceiling anymore."
  echo "Spending guard: stops extending (flags *_BUDGET_GUARD_STOPPED) once estimated billing hits \$${MAX_SPEND_USD}. 'watch' tells you when every system is finished."
  echo "Next:  ./gpuai_resume.sh status     |     ./gpuai_resume.sh watch     |     ./gpuai_resume.sh fetch"
}

cmd_status() {
  load_state
  resolve_conn || exit 1
  local info status rate last_reach now hrs billed
  info=$(gpu instances get "$GPUAI_ID" --output json 2>/dev/null || echo '{}')
  status=$(jq -r '.status // "unknown"' <<<"$info")
  rate=$(jq -r '.price_per_hour // empty' <<<"$info"); rate="${rate:-$GPUAI_RATE}"
  last_reach=$(jq -r '.last_reachable_at // "not probed yet"' <<<"$info")
  now=$(date +%s); hrs=$(awk -v n="$now" -v s="$LAUNCH_EPOCH" 'BEGIN{printf "%.2f", (n-s)/3600}')
  billed=$(awk -v h="$hrs" -v r="$rate" 'BEGIN{printf "%.2f", h*r}')   # per-second, no rounding
  echo "Instance ${GPUAI_ID}: status=${status}  ~${hrs} h since launch, ~\$${billed} billed so far at \$${rate}/hr (cap \$${MAX_SPEND_USD}; per-second billing until deleted)."
  echo "Last reachable (GPU.ai's own SSH probe): ${last_reach}"
  [[ "$status" == "unreachable" ]] && echo "!! GPU.ai reports the tunnel as UNREACHABLE -- the instance is still provisioned/metered but SSH isn't responding. See 'gpu instances get $GPUAI_ID' for details."
  echo
  for sys in "${READY[@]}"; do
    fnm_var="FNM_${sys}"; fnm="${!fnm_var}"
    ssh_run "cd ~/md/$sys 2>/dev/null && { if kill -0 \$(cat mdrun.pid) 2>/dev/null; then A=RUNNING; else A=STOPPED; fi; NS=\$(awk '/^ *Step *Time\$/{getline; t=\$2; f=1} END{if(f) printf \"%.2f\", t/1000; else print 0}' $fnm.log); echo \"$sys  \$A  \${NS} ns (first convergence check at ${MIN_NS})\"; ls *.flag *_converged_frame_ns.txt 2>/dev/null | sed 's/^/    /'; tail -n 2 watchdog.log 2>/dev/null | sed 's/^/    /'; }" || echo "$sys: could not read status"
  done
}

cmd_balance() {
  ensure_cli
  # CONFIRMED on this CLI build: `gpu billing balance` isn't implemented -- it falls through
  # to the `billing` parent command's help text and still exits 0, so it can't be detected
  # reliably by exit code. Don't call it; go straight to what actually works.
  echo "=== GPU.ai account: deposits (CLI) ==="
  gpu billing deposits || true
  echo
  echo "Live balance/burn/runway isn't exposed by 'gpu billing' on this CLI build -- check"
  echo "https://gpu.ai/cloud -> Billing for the current number. With Auto-Pay off, a running"
  echo "instance gets STOPPED (not gracefully paused) the moment balance hits \$0 -- make sure"
  echo "MAX_SPEND_USD leaves real headroom under your actual balance."
}

cmd_watch() {
  load_state
  resolve_conn || exit 1
  declare -A seen
  echo "Polling every ${POLL_SEC}s. Ctrl+C stops WATCHING only; the jobs keep running on the instance."
  while true; do
    local all_done=1 sys flags f
    for sys in "${READY[@]}"; do
      flags=$(ssh_run "ls ~/md/$sys/*.flag 2>/dev/null" || true)
      for f in $flags; do
        if [[ -z "${seen[$f]:-}" ]]; then
          seen[$f]=1
          echo "[$(date +%H:%M)] $sys: $(basename "$f")"
          command -v notify-send >/dev/null && notify-send "ROOTS $sys" "$(basename "$f")" || printf '\a'
        fi
      done
      [[ "$flags" == *"DONE_"* || "$flags" == *"STOPPED_EARLY_"* ]] || all_done=0
    done
    if [[ "$all_done" -eq 1 ]]; then echo "All systems finished (or stopped). Run: ./gpuai_resume.sh fetch"; break; fi
    sleep "$POLL_SEC"
  done
}

cmd_fetch() {
  load_state
  resolve_conn || exit 1
  local sys fnm out ns
  mkdir -p "$OUT_ROOT"
  for sys in "${READY[@]}"; do
    fnm_var="FNM_${sys}"; fnm="${!fnm_var}"; out="$OUT_ROOT/$sys"; mkdir -p "$out"
    echo "--- $sys -> $out"
    ssh_run "cd ~/md/$sys && tar -cf - \$(ls ${fnm}.* ${fnm}_prev.cpt ${fnm}_seg*.tpr watchdog.log convergence_summary.txt rmsd_*.xvg *.flag *_converged_frame_ns.txt resume_seg*.stdout extend_seg*.log run.env 2>/dev/null)" | tar -xf - -C "$out"
    ns=$(last_ns "$out/$fnm.log")
    echo "    local log reaches ${ns} ns"
    [[ -f "$out/${sys}_converged_frame_ns.txt" ]] && echo "    converged window starts at $(cat "$out/${sys}_converged_frame_ns.txt") ns"
    ls "$out"/ALERT_*.flag 2>/dev/null | sed 's/^/    ALERT: /' || true
  done
  echo
  echo "Copied to $OUT_ROOT. Check each system's ns above and any ALERT flags, THEN delete the instance."
  echo "Only run 'delete' after confirming each system's .xtc, .edr, .log and .cpt are present locally."
}

cmd_delete() {
  ensure_cli
  load_state
  echo "Terminating $GPUAI_ID ..."
  gpu instances delete "$GPUAI_ID"
  echo "Deleted. Billing for this instance stops now. State file left at $STATE_FILE for reference."
}

case "${1:-help}" in
  types)     cmd_types ;;
  create)    cmd_create ;;
  preflight) cmd_preflight ;;
  launch)    cmd_launch ;;
  status)    cmd_status ;;
  balance)   cmd_balance ;;
  watch)     cmd_watch ;;
  fetch)     cmd_fetch ;;
  delete)    cmd_delete ;;
  *) sed -n '2,70p' "$0" ;;
esac
