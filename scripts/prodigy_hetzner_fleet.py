#!/usr/bin/env python3
"""
ROOTS -- PRODIGY-LIG cross-validation on a disposable Hetzner CPU fleet.

Workflow:
    plan -> provision -> deploy -> progress -> pull -> delete

The fleet runs prodigy_cross_validation.py against already-docked poses.
Each instance receives a round-robin assignment of the original docking-result
archives, preserving the roots-docking-*/results_*.tar.gz directory layout.

IMPORTANT:
The estimated scoring rate is a local sandbox estimate, not a verified
measurement on Hetzner hardware. Check the plan and account limits before
provisioning. Provisioned servers incur charges until deleted.
"""

import argparse
import json
import os
import subprocess
import sys
import tarfile
import time
from pathlib import Path

PRODIGY_DIR = Path.home() / "roots_pipeline1" / "prodigy_fleet"
IDS_FILE = PRODIGY_DIR / "prodigy_server_ids.json"
IPS_FILE = PRODIGY_DIR / "prodigy_server_ips.json"
CHUNKS_DIR = PRODIGY_DIR / "chunks"

LOCAL_MEASURED_RATE = 10.0  # Local sandbox estimate only: poses/second/worker.

SSH_OPTS_BASE = [
    "-o", "StrictHostKeyChecking=no",
    "-o", "UserKnownHostsFile=/dev/null",
    "-o", "GlobalKnownHostsFile=/dev/null",
    "-o", "ConnectTimeout=15",
    "-o", "BatchMode=yes",
]

# A location-specific server-type availability error is treated as a reason
# to try the next location, not as a reason to abandon provisioning entirely.
RETRYABLE_LOCATION_ERRORS = {
    "resource_unavailable",
    "server_type_not_available_in_location",
}

DEPLOY_SETUP = r"""
set -e
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq
apt-get install -y -qq openbabel python3-pip > /dev/null
pip install --break-system-packages -q biopython prodigy-lig rdkit meeko scipy gemmi
echo "SETUP_OK"
"""

# Keep the launch as one inline command so that backgrounding and disowning
# happen in the same remote shell.
def run_script_cmd(workers):
    return (
        "cd /root && "
        "nohup python3 -u prodigy_cross_validation.py "
        "--receptor-pdbqt receptor.pdbqt "
        "--results-dir fleet_input "
        "--out prodigy_scores.csv "
        f"--workers {workers} "
        "> run.log 2>&1 < /dev/null & "
        "disown; sleep 1; "
        "pgrep -f '[p]rodigy_cross_validation.py' >/dev/null "
        "&& echo LAUNCH=RUNNING || echo LAUNCH=FAILED"
    )


PROGRESS_CHECK = r"""
cd /root
DONE=$(wc -l < prodigy_scores.csv 2>/dev/null || echo 0)
if [ "$DONE" -gt 0 ]; then DONE=$((DONE - 1)); fi
echo "DONE=$DONE"
pgrep -f '[p]rodigy_cross_validation.py' >/dev/null && echo "RUNNING=yes" || echo "RUNNING=no"
tail -n 3 run.log 2>/dev/null | sed 's/^/LOG: /'
"""

PACK_RESULTS = r"""
set -e
cd /root
test -f prodigy_scores.csv
tar -czf prodigy_results.tar.gz prodigy_scores.csv run.log
sha256sum prodigy_results.tar.gz > SHA256SUMS
echo "SIZE=$(stat -c %s prodigy_results.tar.gz)"
"""


def ssh(ip, key_path, script=None, cmd="bash -s", timeout=900):
    """Run a command or send a shell script to a remote server."""
    return subprocess.run(
        [
            "ssh", "-i", str(key_path),
            *SSH_OPTS_BASE,
            f"root@{ip}",
            cmd,
        ],
        input=script,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def kv(text):
    """Parse simple KEY=VALUE status lines."""
    return {
        line.split("=", 1)[0]: line.split("=", 1)[1]
        for line in text.splitlines()
        if "=" in line and not line.startswith("LOG:")
    }


def upload(ip, key_path, local, remote, tries=4):
    """Upload a local file using rsync, retrying transient SSH failures."""
    ssh_command = f"ssh -i {key_path} " + " ".join(SSH_OPTS_BASE)
    result = None

    for attempt in range(tries):
        result = subprocess.run(
            [
                "rsync", "-a", "--partial", "-e", ssh_command,
                str(local), f"root@{ip}:{remote}",
            ],
            capture_output=True,
            text=True,
        )

        if result.returncode == 0:
            return

        if attempt < tries - 1:
            time.sleep(3 * (attempt + 1))

    raise RuntimeError(
        f"Upload to {ip}:{remote} failed: "
        f"{(result.stderr if result else '')[-300:]}"
    )


def download(ip, key_path, remote, local, tries=4):
    """Download a remote file using rsync, retrying transient SSH failures."""
    local = Path(local)
    local.parent.mkdir(parents=True, exist_ok=True)

    ssh_command = f"ssh -i {key_path} " + " ".join(SSH_OPTS_BASE)
    result = None

    for attempt in range(tries):
        result = subprocess.run(
            [
                "rsync", "-a", "--partial", "-e", ssh_command,
                f"root@{ip}:{remote}", str(local),
            ],
            capture_output=True,
            text=True,
        )

        if result.returncode == 0:
            return

        if attempt < tries - 1:
            time.sleep(3 * (attempt + 1))

    raise RuntimeError(
        f"Download from {ip}:{remote} failed: "
        f"{(result.stderr if result else '')[-300:]}"
    )

def api_request(token, method, url, payload=None, timeout=60):
    """
    Make a Hetzner API request and parse its response.

    Empty responses are accepted for endpoints such as successful DELETE
    requests. HTTP error responses remain inspectable as JSON where provided.
    """
    command = [
        "curl", "-sS",
        "--max-time", str(timeout),
        "-X", method,
        "-H", f"Authorization: Bearer {token}",
    ]

    if payload is not None:
        command.extend([
            "-H", "Content-Type: application/json",
            "-d", json.dumps(payload),
        ])

    command.append(url)

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=timeout + 10,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"Hetzner API request failed: {result.stderr.strip()[-300:]}"
        )

    if not result.stdout.strip():
        return {}

    try:
        return json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError(
            f"Hetzner API returned invalid JSON: {result.stdout[-300:]!r}"
        ) from exc


def count_total_poses(results_dir):
    """Count docked poses from the pulled docking-fleet archives."""
    total = 0

    for server_dir in sorted(Path(results_dir).glob("roots-docking-*")):
        archive = next(server_dir.glob("results_*.tar.gz"), None)

        if archive is None:
            continue

        with tarfile.open(archive, "r:gz") as tar:
            total += sum(
                1 for member in tar.getmembers()
                if member.name.endswith("_docked.pdbqt")
            )

    return total


def cmd_plan(args):
    results_dir = Path(args.results_dir)

    if not results_dir.exists():
        sys.exit(f"Results directory not found: {results_dir}")

    total = count_total_poses(results_dir)

    print(f"Total docked poses found across pulled fleet results: {total:,}")

    if total == 0:
        sys.exit(
            "No poses found. Check that --results-dir points to the "
            "pull_results.py output directory."
        )

    print("\n=== ESTIMATE: LOCAL SANDBOX RATE, NOT HETZNER MEASURED ===")
    print(f"Estimated rate: ~{LOCAL_MEASURED_RATE:.0f} poses/s/worker")
    print(
        "The actual rate on Hetzner CX CPUs has not been verified. "
        "Treat the ETA and cost as planning estimates."
    )

    total_workers = args.n_instances * args.cores_per_instance

    if total_workers <= 0:
        sys.exit("The number of instances and cores per instance must be positive.")

    seconds = total / (LOCAL_MEASURED_RATE * total_workers)
    hours = seconds / 3600

    print(
        f"\n{args.n_instances} instances x {args.cores_per_instance} cores "
        f"= {total_workers} workers"
    )
    print(f"Estimated time: {hours:.1f} hours ({hours / 24:.1f} days)")

    # This is the configured historical estimate, not a live price quote.
    if args.server_type == "cx33":
        hourly_price = 0.016
        cost = args.n_instances * hours * hourly_price
        print(
            f"Estimated compute cost at the configured EUR {hourly_price:.3f}/h "
            f"per instance: EUR {cost:.2f}"
        )
        print("Verify the current price and billing terms in your Hetzner account.")


def cmd_provision(args):
    if not args.ssh_key_id:
        sys.exit("--ssh-key-id is required for provisioning.")

    if args.n_instances < 1:
        sys.exit("--n-instances must be at least 1.")

    if args.provision_retries < 1:
        sys.exit("--provision-retries must be at least 1.")

    token = os.environ.get("HETZNER_TOKEN")

    if not token:
        sys.exit("HETZNER_TOKEN environment variable is not set.")

    PRODIGY_DIR.mkdir(parents=True, exist_ok=True)

    try:
        account_data = api_request(
            token,
            "GET",
            "https://api.hetzner.cloud/v1/servers",
        )
    except (RuntimeError, subprocess.SubprocessError) as exc:
        sys.exit(f"Could not inspect existing Hetzner servers: {exc}")

    existing = sum(
        server.get("server_type", {}).get("cores", 0)
        for server in account_data.get("servers", [])
        if server.get("server_type", {}).get("cpu_type") == "shared"
    )

    default_cores = {"cx23": 2, "cx33": 4, "cx43": 8}
    cores_per = default_cores.get(args.server_type, args.cores_per_instance)
    requested = cores_per * args.n_instances

    print(f"Shared cores currently in use: {existing} | this request: {requested}")

    if existing + requested > 20:
        print(
            "WARNING: this request exceeds the previously observed 20-shared-core "
            "account ceiling. Provisioning may fail. Check your current account limits."
        )

    # Prefer the requested location, then try the other configured locations.
    # Some server types are not offered in every location; that error is handled
    # as a location fallback rather than a fatal authentication/quota error.
    locations = list(dict.fromkeys(
        [args.location] + [
            location
            for location in ("fsn1", "hel1", "nbg1", "ash", "hil")
            if location != args.location
        ]
    ))

    def create_one(name):
        for attempt in range(args.provision_retries):
            for location in locations:
                payload = {
                    "name": name,
                    "server_type": args.server_type,
                    "image": "ubuntu-24.04",
                    "location": location,
                    "ssh_keys": [args.ssh_key_id],
                }

                try:
                    response = api_request(
                        token,
                        "POST",
                        "https://api.hetzner.cloud/v1/servers",
                        payload=payload,
                    )
                except (RuntimeError, subprocess.SubprocessError) as exc:
                    print(f"FAILED for {name} in {location}: {exc}")
                    return None

                if "server" in response:
                    if location != args.location:
                        print(
                            f"  Fallback placement: {location} "
                            f"(preferred: {args.location})"
                        )
                    return response["server"]["id"]

                error = response.get("error", {})
                error_code = error.get("code", "unknown_error")

                if error_code not in RETRYABLE_LOCATION_ERRORS:
                    print(
                        f"FAILED for {name} in {location}: "
                        f"{error or response}"
                    )
                    return None

                print(
                    f"  {name}: {error_code} in {location}; "
                    "trying another location."
                )

            if attempt < args.provision_retries - 1:
                wait_seconds = min(5 * (2 ** attempt), 60)
                print(
                    f"  No placement available for {name}; retrying in "
                    f"{wait_seconds}s "
                    f"(attempt {attempt + 1}/{args.provision_retries})."
                )
                time.sleep(wait_seconds)

        print(
            f"FAILED for {name}: no capacity after "
            f"{args.provision_retries} attempts across {locations}"
        )
        return None

    server_ids = []

    for index in range(args.n_instances):
        name = f"roots-prodigy-{index:02d}"
        server_id = create_one(name)

        if server_id:
            server_ids.append(server_id)
            print(f"Created {name}: id={server_id}")
        else:
            print(f"Could not create {name}.")

    # Save IDs even if provisioning was partial, so created servers are recorded.
    IDS_FILE.write_text(json.dumps(server_ids, indent=2))

    if not server_ids:
        IPS_FILE.write_text("{}")
        sys.exit(
            "\nNo servers were created. No boot wait was performed. "
            "Check the API errors, server-type availability, SSH key ID, "
            "account limits, and locations before retrying."
        )

    print(f"\nCreated {len(server_ids)}/{args.n_instances} requested servers.")
    print("Waiting 30 seconds for initial boot...")
    time.sleep(30)

    server_ips = {}

    for server_id in server_ids:
        try:
            response = api_request(
                token,
                "GET",
                f"https://api.hetzner.cloud/v1/servers/{server_id}",
            )
        except (RuntimeError, subprocess.SubprocessError) as exc:
            print(f"WARNING: could not inspect server {server_id}: {exc}")
            continue

        server = response.get("server")

        if not server:
            print(f"WARNING: no server object returned for ID {server_id}.")
            continue

        name = server["name"]
        public_net = server.get("public_net", {})
        ipv4 = public_net.get("ipv4", {}).get("ip")

        if not ipv4:
            print(f"WARNING: {name} has no public IPv4 address yet.")
            continue

        server_ips[name] = ipv4
        print(f"{name}: {ipv4} ({server.get('status', 'unknown status')})")

    IPS_FILE.write_text(json.dumps(server_ips, indent=2))

    print(f"\nSaved server IDs: {IDS_FILE}")
    print(f"Saved server IPs: {IPS_FILE}")

    if len(server_ips) != len(server_ids):
        print(
            "WARNING: Some created servers have no recorded IP address. "
            "Inspect the Hetzner console before continuing."
        )

    if len(server_ids) != args.n_instances:
        print(
            "WARNING: Provisioning was partial. The successfully created "
            "servers are recorded above; check your intended fleet size "
            "before deployment."
        )


def cmd_deploy(args):
    if not IPS_FILE.exists():
        sys.exit(f"{IPS_FILE} not found. Run --action provision first.")

    if not args.results_dir:
        sys.exit("--results-dir is required for deployment.")

    if not args.receptor_pdbqt:
        sys.exit("--receptor-pdbqt is required for deployment.")

    if args.workers < 1:
        sys.exit("--workers must be at least 1.")

    receptor = Path(args.receptor_pdbqt).expanduser()
    script_path = Path(args.script_path).expanduser()
    results_dir = Path(args.results_dir).expanduser()

    if not receptor.is_file():
        sys.exit(f"Receptor not found: {receptor}")

    if not script_path.is_file():
        sys.exit(
            f"prodigy_cross_validation.py not found at {script_path}. "
            "Pass --script-path."
        )

    if not results_dir.is_dir():
        sys.exit(f"Results directory not found: {results_dir}")

    server_ips = json.loads(IPS_FILE.read_text())
    servers = sorted(server_ips.items())
    n_servers = len(servers)

    if n_servers == 0:
        sys.exit(
            "No provisioned server IPs are recorded. "
            "Provision servers successfully before deploying."
        )

    source_dirs = sorted(results_dir.glob("roots-docking-*"))

    if not source_dirs:
        sys.exit(
            f"No roots-docking-* directories found in {results_dir}."
        )

    if not any(
        next(directory.glob("results_*.tar.gz"), None)
        for directory in source_dirs
    ):
        sys.exit(
            f"No results_*.tar.gz archives found in {results_dir}/roots-docking-*."
        )

    print(
        f"Splitting {len(source_dirs)} pulled result directories "
        f"across {n_servers} instances."
    )

    if len(source_dirs) != n_servers:
        print(
            "NOTE: The number of result directories differs from the number "
            "of instances. Archives will be assigned round-robin."
        )

    CHUNKS_DIR.mkdir(parents=True, exist_ok=True)

    for index, (name, ip) in enumerate(servers):
        assigned = [
            source_dirs[j]
            for j in range(len(source_dirs))
            if j % n_servers == index
        ]

        print(f"\n=== {name} ({ip}) <- {[directory.name for directory in assigned]} ===")

        # Install dependencies on this instance.
        try:
            result = ssh(ip, args.ssh_key_path, DEPLOY_SETUP, timeout=900)
        except (subprocess.SubprocessError, OSError) as exc:
            print(f"  SETUP FAILED: {exc}")
            continue

        if result.returncode != 0 or "SETUP_OK" not in result.stdout:
            print(
                f"  SETUP FAILED: "
                f"{(result.stderr or result.stdout)[-500:]}"
            )
            continue

        print("  Dependencies installed.")

        try:
            upload(ip, args.ssh_key_path, receptor, "/root/receptor.pdbqt")
            upload(
                ip,
                args.ssh_key_path,
                script_path,
                "/root/prodigy_cross_validation.py",
            )

            # Create every remote directory BEFORE uploading its archive.
            result = ssh(
                ip,
                args.ssh_key_path,
                cmd="mkdir -p /root/fleet_input",
            )

            if result.returncode != 0:
                raise RuntimeError(
                    f"Could not create /root/fleet_input: "
                    f"{result.stderr[-300:]}"
                )

            for directory in assigned:
                remote_dir = f"/root/fleet_input/{directory.name}"

                result = ssh(
                    ip,
                    args.ssh_key_path,
                    cmd=f"mkdir -p {remote_dir}",
                )

                if result.returncode != 0:
                    raise RuntimeError(
                        f"Could not create remote directory {remote_dir}: "
                        f"{result.stderr[-300:]}"
                    )

                archive = next(directory.glob("results_*.tar.gz"), None)

                if archive is None:
                    print(
                        f"  WARNING: no results archive found in {directory}; "
                        "this directory will be skipped."
                    )
                    continue

                upload(
                    ip,
                    args.ssh_key_path,
                    archive,
                    f"{remote_dir}/{archive.name}",
                )

            print("  Receptor, scoring script, and assigned archives uploaded.")

        except (RuntimeError, subprocess.SubprocessError, OSError) as exc:
            print(f"  DEPLOYMENT FAILED: {exc}")
            continue

        try:
            result = ssh(
                ip,
                args.ssh_key_path,
                cmd=run_script_cmd(args.workers),
            )
        except (subprocess.SubprocessError, OSError) as exc:
            print(f"  LAUNCH FAILED: {exc}")
            continue

        status = kv(result.stdout)

        if result.returncode != 0:
            print(f"  LAUNCH COMMAND FAILED: {result.stderr[-300:]}")
        else:
            print(f"  {status.get('LAUNCH', 'UNKNOWN')}")


def cmd_progress(args):
    if not IPS_FILE.exists():
        sys.exit(f"{IPS_FILE} not found. Run --action provision first.")

    server_ips = json.loads(IPS_FILE.read_text())

    if not server_ips:
        sys.exit("No server IPs are recorded. Check provisioning first.")

    total_done = 0

    for name, ip in sorted(server_ips.items()):
        try:
            result = ssh(ip, args.ssh_key_path, PROGRESS_CHECK, timeout=60)
        except (subprocess.SubprocessError, OSError) as exc:
            print(f"{name} ({ip}): STATUS CHECK FAILED: {exc}")
            continue

        status = kv(result.stdout)

        try:
            done = int(status.get("DONE", 0))
        except ValueError:
            done = 0

        total_done += done

        print(
            f"{name} ({ip}): {done:,} scored | "
            f"running={status.get('RUNNING', '?')}"
        )

        if result.returncode != 0:
            print(f"  SSH error: {result.stderr[-200:]}")

        for line in result.stdout.splitlines():
            if line.startswith("LOG:"):
                print(f"  {line}")

    print(f"\nTotal scored so far across fleet: {total_done:,}")


def cmd_pull(args):
    if not IPS_FILE.exists():
        sys.exit(f"{IPS_FILE} not found. Run --action provision first.")

    server_ips = json.loads(IPS_FILE.read_text())

    if not server_ips:
        sys.exit("No server IPs are recorded. Check provisioning first.")

    out_dir = Path(args.local_out).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)

    all_ok = True

    for name, ip in sorted(server_ips.items()):
        server_out = out_dir / name
        server_out.mkdir(parents=True, exist_ok=True)

        try:
            result = ssh(ip, args.ssh_key_path, PACK_RESULTS, timeout=300)
        except (subprocess.SubprocessError, OSError) as exc:
            print(f"{name}: FAILED to pack results: {exc}")
            all_ok = False
            continue

        status = kv(result.stdout)

        if result.returncode != 0 or "SIZE" not in status:
            print(
                f"{name}: FAILED to pack results: "
                f"{(result.stderr or result.stdout)[-300:]}"
            )
            all_ok = False
            continue

        print(
            f"{name}: packed ({int(status['SIZE']) / 1e6:.1f} MB), "
            "downloading..."
        )

        archive_path = server_out / "prodigy_results.tar.gz"
        checksum_path = server_out / "SHA256SUMS"

        try:
            download(
                ip, args.ssh_key_path,
                "/root/prodigy_results.tar.gz",
                archive_path,
            )
            download(
                ip, args.ssh_key_path,
                "/root/SHA256SUMS",
                checksum_path,
            )
        except (RuntimeError, subprocess.SubprocessError, OSError) as exc:
            print(f"  DOWNLOAD FAILED: {exc}")
            all_ok = False
            continue

        try:
            expected_hash = checksum_path.read_text().split()[0]
            checksum_result = subprocess.run(
                ["sha256sum", str(archive_path)],
                capture_output=True,
                text=True,
                check=True,
            )
            actual_hash = checksum_result.stdout.split()[0]
        except (OSError, IndexError, subprocess.SubprocessError) as exc:
            print(f"  CHECKSUM FAILED: {exc}")
            all_ok = False
            continue

        if actual_hash != expected_hash:
            print("  VERIFICATION FAILED: checksum mismatch.")
            all_ok = False
            continue

        try:
            with tarfile.open(archive_path, "r:gz") as tar:
                # Only extract the two expected files from our own result archive.
                allowed_names = {"prodigy_scores.csv", "run.log"}
                members = [
                    member for member in tar.getmembers()
                    if member.name in allowed_names and member.isfile()
                ]

                if not any(m.name == "prodigy_scores.csv" for m in members):
                    raise RuntimeError(
                        "Archive does not contain prodigy_scores.csv."
                    )

                tar.extractall(server_out, members=members)

            (server_out / "VERIFIED.json").write_text(
                json.dumps(
                    {
                        "server": name,
                        "ip": ip,
                        "sha256": actual_hash,
                    },
                    indent=2,
                )
            )

        except (OSError, tarfile.TarError, RuntimeError) as exc:
            print(f"  EXTRACTION/VERIFICATION FAILED: {exc}")
            all_ok = False
            continue

        print("  VERIFIED")

    if all_ok:
        print(
            "\nAll recorded servers passed download verification. "
            "Review the downloaded CSVs and confirm that scoring has finished "
            "before running --action delete."
        )
    else:
        print(
            "\nSome servers were not verified. Do not delete the fleet "
            "until the missing results have been recovered or you explicitly "
            "decide to discard them."
        )


def cmd_delete(args):
    if not IDS_FILE.exists() or not IPS_FILE.exists():
        sys.exit("No provisioned fleet found on record.")

    token = os.environ.get("HETZNER_TOKEN")

    if not token:
        sys.exit("HETZNER_TOKEN environment variable is not set.")

    server_ips = json.loads(IPS_FILE.read_text())
    pull_dir = Path(args.local_out).expanduser() if args.local_out else None

    for name, ip in sorted(server_ips.items()):
        if pull_dir and not (pull_dir / name / "VERIFIED.json").exists():
            print(
                f"SKIPPING {name}: no VERIFIED.json found in "
                f"{pull_dir / name}."
            )

            if not args.force:
                print(
                    "Pass --force only if you intentionally accept the risk "
                    "of losing results."
                )
                continue

        try:
            response = api_request(
                token,
                "GET",
                "https://api.hetzner.cloud/v1/servers",
            )
        except (RuntimeError, subprocess.SubprocessError) as exc:
            print(f"{name}: could not query servers: {exc}")
            continue

        servers = [
            server for server in response.get("servers", [])
            if server.get("name") == name
        ]

        if not servers:
            print(f"{name}: not found (possibly already deleted).")
            continue

        server_id = servers[0]["id"]

        try:
            delete_response = api_request(
                token,
                "DELETE",
                f"https://api.hetzner.cloud/v1/servers/{server_id}",
            )
        except (RuntimeError, subprocess.SubprocessError) as exc:
            print(f"{name} (id {server_id}): deletion request failed: {exc}")
            continue

        # Successful Hetzner DELETE requests commonly return an empty body.
        # api_request expects JSON, so the API wrapper may reject that response.
        # This branch is retained for API responses that include JSON.
        print(f"{name} (id {server_id}): delete response: {delete_response}")


def main():
    parser = argparse.ArgumentParser(
        description="ROOTS PRODIGY-LIG Hetzner fleet orchestration."
    )

    parser.add_argument(
        "--action",
        required=True,
        choices=["plan", "provision", "deploy", "progress", "pull", "delete"],
    )
    parser.add_argument(
        "--results-dir",
        help="Directory containing roots-docking-*/ subfolders from pull_results.py.",
    )
    parser.add_argument("--receptor-pdbqt")
    parser.add_argument(
        "--script-path",
        default=str(Path.home() / "prodigy_cross_validation.py"),
    )
    parser.add_argument(
        "--ssh-key-path",
        default=str(Path.home() / ".ssh" / "hetzner_step4"),
    )
    parser.add_argument("--ssh-key-id", type=int)
    parser.add_argument(
        "--location",
        default="nbg1",
        help="Preferred Hetzner location; fallback locations will be tried.",
    )
    parser.add_argument(
        "--provision-retries",
        type=int,
        default=4,
        help="Number of capacity-placement passes per instance.",
    )
    parser.add_argument("--n-instances", type=int, default=5)
    parser.add_argument("--server-type", default="cx33")
    parser.add_argument("--cores-per-instance", type=int, default=4)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--local-out",
        default=str(Path.home() / "roots_pipeline1" / "prodigy_results"),
    )
    parser.add_argument("--force", action="store_true")

    args = parser.parse_args()

    actions = {
        "plan": cmd_plan,
        "provision": cmd_provision,
        "deploy": cmd_deploy,
        "progress": cmd_progress,
        "pull": cmd_pull,
        "delete": cmd_delete,
    }

    try:
        actions[args.action](args)
    except KeyboardInterrupt:
        sys.exit("\nInterrupted by user.")
    except (OSError, subprocess.SubprocessError, tarfile.TarError) as exc:
        sys.exit(f"Operation failed: {exc}")


if __name__ == "__main__":
    main()