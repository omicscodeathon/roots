#!/usr/bin/env python3
"""
ROOTS -- Step 6: High-Throughput Docking Fleet Orchestration (Layer 2)

Splits the compound library into N chunks, provisions N Hetzner CX-tier
cloud instances via the API, deploys the QuickVina2 worker and each
instance's chunk, and launches docking in the background on all of them
simultaneously.

This is a genuinely parallel workload by design: each compound-receptor
docking run is fully independent and stateless, so splitting the library
across many small, cheap instances running concurrently is both correct
and, per Hetzner's shared-core account limits, the only way to achieve
real wall-clock speedup without a large, continuously-billed dedicated
instance.

Prerequisites:
    - HETZNER_TOKEN environment variable set (a fresh API token -- never
      commit one to source control or paste it into a shared session)
    - An SSH key already registered on the Hetzner account; its numeric ID
      set in SSH_KEY_ID below
    - step6_dock_worker.sh in the same directory as this script
    - A prepared ALDH3B2 receptor PDBQT and the full ligand PDBQT library
      available locally to distribute to each instance

Usage:
    export HETZNER_TOKEN="your_token_here"
    python3 step6_orchestrate_fleet.py --n-instances 5 --server-type cx33
"""

import argparse
import json
import os
import subprocess
import time
from pathlib import Path


def check_shared_core_limit(token, requested_cores):
    """Hetzner accounts have a shared-vCPU ceiling across all running
    instances combined (observed limit: 20 shared cores on a standard
    account). Check current usage before requesting more and hitting a
    partial, confusing failure partway through provisioning."""
    result = subprocess.run(
        ["curl", "-s", "-H", f"Authorization: Bearer {token}",
         "https://api.hetzner.cloud/v1/servers"],
        capture_output=True, text=True,
    )
    data = json.loads(result.stdout)
    current_shared_cores = sum(
        s["server_type"]["cores"] for s in data.get("servers", [])
        if s["server_type"]["cpu_type"] == "shared"
    )
    print(f"Current shared cores in use: {current_shared_cores}")
    print(f"This request would add: {requested_cores}")
    if current_shared_cores + requested_cores > 20:
        print(
            "WARNING: this may exceed the account's shared-core limit "
            "(observed ceiling: 20). Provisioning will likely fail partway "
            "through with 'resource_limit_exceeded' -- reduce --n-instances "
            "or --server-type, or delete unused instances first."
        )
    return current_shared_cores


def provision_fleet(token, n_instances, server_type, ssh_key_id, location="nbg1"):
    server_ids = []
    for i in range(n_instances):
        payload = {
            "name": f"roots-docking-{i:02d}",
            "server_type": server_type,
            "image": "ubuntu-24.04",
            "location": location,
            "ssh_keys": [ssh_key_id],
        }
        result = subprocess.run(
            ["curl", "-s", "-X", "POST",
             "-H", f"Authorization: Bearer {token}",
             "-H", "Content-Type: application/json",
             "-d", json.dumps(payload),
             "https://api.hetzner.cloud/v1/servers"],
            capture_output=True, text=True,
        )
        resp = json.loads(result.stdout)
        if "server" in resp:
            sid = resp["server"]["id"]
            server_ids.append(sid)
            print(f"Created roots-docking-{i:02d}: id={sid}")
        else:
            print(f"FAILED for instance {i}: {resp.get('error', resp)}")
    return server_ids


def get_server_ips(token, server_ids, boot_wait_s=30):
    print(f"Waiting {boot_wait_s}s for boot...")
    time.sleep(boot_wait_s)
    server_ips = {}
    for sid in server_ids:
        result = subprocess.run(
            ["curl", "-s", "-H", f"Authorization: Bearer {token}",
             f"https://api.hetzner.cloud/v1/servers/{sid}"],
            capture_output=True, text=True,
        )
        resp = json.loads(result.stdout)
        name = resp["server"]["name"]
        ip = resp["server"]["public_net"]["ipv4"]["ip"]
        status = resp["server"]["status"]
        server_ips[name] = ip
        print(f"{name}: {ip} ({status})")
    return server_ips


def split_library(local_ligands_dir, n_chunks, chunks_dir):
    chunks_dir.mkdir(parents=True, exist_ok=True)
    ligand_files = sorted(os.listdir(local_ligands_dir))
    total = len(ligand_files)
    chunk_size = (total + n_chunks - 1) // n_chunks

    print(f"Total ligands: {total:,}")
    print(f"Chunks: {n_chunks}, ~{chunk_size:,} compounds each")

    for i in range(n_chunks):
        start = i * chunk_size
        end = min(start + chunk_size, total)
        chunk_files = ligand_files[start:end]
        chunk_path = chunks_dir / f"chunk_{i:02d}.txt"
        chunk_path.write_text("\n".join(chunk_files))
        print(f"  chunk_{i:02d}.txt: {len(chunk_files):,} compounds")
    return total


def deploy_and_launch(ssh_key_path, server_ips, worker_script, receptor_pdbqt,
                       local_ligands_dir, chunks_dir):
    for i, (name, ip) in enumerate(sorted(server_ips.items())):
        chunk_id = f"{i:02d}"
        print(f"\n=== Deploying to {name} ({ip}) ===")

        ssh_base = ["ssh", "-i", ssh_key_path, "-o", "StrictHostKeyChecking=no", f"root@{ip}"]
        scp_base = ["scp", "-i", ssh_key_path, "-o", "StrictHostKeyChecking=no"]

        subprocess.run(ssh_base + [
            "mkdir -p /root/step5_docking/receptors /root/docking_results /root/docking_logs"
        ], check=True)

        subprocess.run(ssh_base + [
            "wget -q https://github.com/QVina/qvina/raw/master/bin/qvina2.1 "
            "-O /usr/local/bin/qvina2.1 && chmod +x /usr/local/bin/qvina2.1"
        ], check=True)
        print("  QuickVina2 installed")

        subprocess.run(scp_base + [str(receptor_pdbqt), f"root@{ip}:/root/step5_docking/receptors/"], check=True)
        print("  Receptor copied")

        subprocess.run(scp_base + [str(worker_script), f"root@{ip}:/root/dock_worker.sh"], check=True)
        subprocess.run(ssh_base + ["chmod +x /root/dock_worker.sh"], check=True)
        print("  Worker script copied")

        chunk_file = chunks_dir / f"chunk_{chunk_id}.txt"
        subprocess.run(scp_base + [str(chunk_file), f"root@{ip}:/root/chunk_{chunk_id}.txt"], check=True)

        compound_files = chunk_file.read_text().strip().split("\n")
        tarball = chunks_dir / f"chunk_{chunk_id}_ligands.tar.gz"
        filelist = chunks_dir / f"chunk_{chunk_id}_filelist.txt"
        filelist.write_text("\n".join(compound_files))
        subprocess.run(f"cd {local_ligands_dir} && tar -czf {tarball} -T {filelist}", shell=True, check=True)
        print(f"  Tarball built: {tarball.stat().st_size / 1e6:.1f} MB")

        subprocess.run(scp_base + [str(tarball), f"root@{ip}:/root/chunk_ligands.tar.gz"], check=True)
        subprocess.run(ssh_base + [
            "mkdir -p /root/ligands_final && tar -xzf /root/chunk_ligands.tar.gz -C /root/ligands_final"
        ], check=True)
        print("  Extracted on remote")

        # Launched individually per instance, NOT in a batched loop with shared
        # state: a prior orchestration attempt using one Python subprocess loop
        # for all launches silently only started the first instance, leaving
        # the rest never launched with no error raised. One explicit,
        # independently-verified ssh call per instance avoids that failure mode.
        launch_cmd = f"cd /root && nohup bash dock_worker.sh {chunk_id} > /root/worker_stdout.log 2>&1 < /dev/null & disown"
        subprocess.run(ssh_base + [launch_cmd])
        print(f"  Launched chunk {chunk_id}")


def verify_fleet(ssh_key_path, server_ips):
    print("\n=== Fleet verification ===")
    all_running = True
    for name, ip in sorted(server_ips.items()):
        ssh_base = ["ssh", "-i", ssh_key_path, "-o", "StrictHostKeyChecking=no", f"root@{ip}"]
        proc_check = subprocess.run(ssh_base + ["pgrep -f dock_worker.sh"], capture_output=True, text=True)
        status = "RUNNING" if proc_check.stdout.strip() else "NOT RUNNING"
        if status == "NOT RUNNING":
            all_running = False
        count_check = subprocess.run(
            ssh_base + ["find /root/docking_results -name '*.pdbqt' 2>/dev/null | wc -l"],
            capture_output=True, text=True,
        )
        count = count_check.stdout.strip()
        print(f"{name} ({ip}): {status} | {count} results so far")
    if not all_running:
        print("\nWARNING: at least one instance is NOT RUNNING. Re-launch it individually:")
        print('  ssh -i <key> root@<ip> "cd /root && nohup bash dock_worker.sh <NN> '
              '> /root/worker_stdout.log 2>&1 < /dev/null & disown"')
    return all_running


def check_fleet_progress(ssh_key_path, server_ips, total_compounds):
    total_done = 0
    print("=== HETZNER FLEET STATUS ===\n")
    for name, ip in sorted(server_ips.items()):
        ssh_base = ["ssh", "-i", ssh_key_path, "-o", "StrictHostKeyChecking=no", f"root@{ip}"]
        count_result = subprocess.run(
            ssh_base + ["find /root/docking_results -name '*.pdbqt' | wc -l"],
            capture_output=True, text=True,
        )
        count = int(count_result.stdout.strip() or 0)
        total_done += count
        complete_check = subprocess.run(
            ssh_base + ["test -f /root/DOCKING_COMPLETE && echo DONE || echo IN_PROGRESS"],
            capture_output=True, text=True,
        )
        completion = complete_check.stdout.strip()
        print(f"{name} ({ip}): {completion} | {count:,} results")
    pct = 100 * total_done / total_compounds if total_compounds else 0
    print(f"\n=== TOTAL: {total_done:,} / {total_compounds:,} ({pct:.1f}%) ===")
    return total_done


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Orchestrate a ROOTS Step 6 QuickVina2 docking fleet on Hetzner")
    parser.add_argument("--n-instances", type=int, default=5)
    parser.add_argument("--server-type", default="cx33")
    parser.add_argument("--ssh-key-id", type=int, required=True)
    parser.add_argument("--ssh-key-path", required=True)
    parser.add_argument("--local-ligands-dir", required=True)
    parser.add_argument("--receptor-pdbqt", required=True)
    parser.add_argument("--worker-script", default="step6_dock_worker.sh")
    parser.add_argument("--chunks-dir", default="./chunks")
    parser.add_argument("--server-ids-file", default="server_ids.json")
    parser.add_argument("--server-ips-file", default="server_ips.json")
    parser.add_argument("--action", choices=["provision", "deploy", "verify", "progress"], required=True,
                         help="provision: create instances. deploy: split+push+launch. "
                              "verify: confirm all instances running. progress: check completion.")
    args = parser.parse_args()

    token = os.environ.get("HETZNER_TOKEN")
    if not token:
        raise SystemExit("HETZNER_TOKEN environment variable not set.")

    chunks_dir = Path(args.chunks_dir)

    if args.action == "provision":
        cores_per_instance = {"cx23": 2, "cx33": 4, "cx43": 8, "cx53": 16}.get(args.server_type, 4)
        check_shared_core_limit(token, cores_per_instance * args.n_instances)
        server_ids = provision_fleet(token, args.n_instances, args.server_type, args.ssh_key_id)
        Path(args.server_ids_file).write_text(json.dumps(server_ids, indent=2))
        server_ips = get_server_ips(token, server_ids)
        Path(args.server_ips_file).write_text(json.dumps(server_ips, indent=2))

    elif args.action == "deploy":
        server_ips = json.loads(Path(args.server_ips_file).read_text())
        total = split_library(args.local_ligands_dir, len(server_ips), chunks_dir)
        deploy_and_launch(args.ssh_key_path, server_ips, args.worker_script,
                           args.receptor_pdbqt, args.local_ligands_dir, chunks_dir)
        time.sleep(15)
        verify_fleet(args.ssh_key_path, server_ips)

    elif args.action == "verify":
        server_ips = json.loads(Path(args.server_ips_file).read_text())
        verify_fleet(args.ssh_key_path, server_ips)

    elif args.action == "progress":
        server_ips = json.loads(Path(args.server_ips_file).read_text())
        total = len(list(Path(args.local_ligands_dir).glob("*.pdbqt"))) if Path(args.local_ligands_dir).exists() else 0
        check_fleet_progress(args.ssh_key_path, server_ips, total)
