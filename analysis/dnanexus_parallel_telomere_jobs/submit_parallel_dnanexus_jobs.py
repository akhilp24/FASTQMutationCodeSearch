#!/usr/bin/env python3
"""Submit telomere CRAM analyses in parallel on DNANexus.

Input: a text file containing one participant ID per line.
Behavior:
1) Resolve each participant's CRAM in the configured UKB folder.
2) Submit one Swiss Army Knife job per participant (parallel fan-out).
3) Wait for all jobs.
4) Download per-sample CSVs and create a combined CSV.
5) Upload combined CSV back to the run folder.
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import List

import pandas as pd


def _run(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    try:
        return subprocess.run(cmd, check=check, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        stderr = (exc.stderr or "").strip()
        stdout = (exc.stdout or "").strip()
        cmd_text = " ".join(cmd)
        msg = f"Command failed (exit {exc.returncode}): {cmd_text}"
        if stdout:
            msg += f"\nstdout:\n{stdout}"
        if stderr:
            msg += f"\nstderr:\n{stderr}"
        raise RuntimeError(msg) from exc


def read_participant_ids(path: str) -> List[str]:
    ids: List[str] = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            ids.append(line)
    if not ids:
        raise ValueError("Participant ID file is empty.")
    return ids


def resolve_cram_file_id(project: str, cram_folder: str, participant_id: str) -> str:
    pattern = f"{participant_id}_*.cram"
    result = _run([
        "dx", "find", "data",
        "--project", project,
        "--folder", cram_folder,
        "--name", pattern,
        "--class", "file",
        "--brief",
    ])
    ids = [x.strip() for x in result.stdout.splitlines() if x.strip()]
    if not ids:
        raise RuntimeError(f"No CRAM found for participant {participant_id} in {cram_folder}.")
    if len(ids) > 1:
        print(f"[warn] Multiple CRAM files found for {participant_id}; using first: {ids[0]}")
    return ids[0]


def upload_worker_script(project: str, folder: str, worker_script: str) -> str:
    res = _run(["dx", "upload", worker_script, "--project", project, "--path", folder, "--brief"])
    file_id = res.stdout.strip()
    if not file_id:
        raise RuntimeError("Failed to upload worker script to DNANexus.")
    return file_id


def ensure_remote_folder(project: str, folder: str) -> None:
    folder = folder if folder.startswith("/") else f"/{folder}"
    _run(["dx", "mkdir", "-p", f"{project}:{folder}"])


def submit_one_job(
    destination: str,
    cram_file_id: str,
    worker_file_id: str,
    reference_fasta_file: str,
    metadata_csv_file: str,
    participant_id: str,
    per_sample_output_folder: str,
    k: int,
    instance_type: str,
) -> str:
    out_name = f"{participant_id}_telomere.csv"
    job_cmd = " && ".join([
        "set -euo pipefail",
        "python3 -m pip install --quiet pysam pandas numpy",
        "CRAM_FILE=$(ls *.cram | head -n 1)",
        "REF_FILE=$(ls *.fa* | head -n 1)",
        "AGE_FILE=$(ls *.csv | head -n 1)",
        "python3 run_single_cram_analysis.py "
        + "--cram \"$CRAM_FILE\" "
        + f"--participant-id {shlex.quote(participant_id)} "
        + "--reference-fasta \"$REF_FILE\" "
        + "--metadata-csv \"$AGE_FILE\" "
        + f"--output-csv {shlex.quote(out_name)} "
        + f"--k {k}",
        f"dx upload {shlex.quote(out_name)} --path {shlex.quote(per_sample_output_folder)} --brief",
    ])

    job = _run([
        "dx", "run", "app-swiss-army-knife",
        "-iin=" + cram_file_id,
        "-iin=" + worker_file_id,
        "-iin=" + reference_fasta_file,
        "-iin=" + metadata_csv_file,
        "-icmd=" + job_cmd,
        "--instance-type", instance_type,
        "--name", f"telomere-{participant_id}",
        "--destination", destination,
        "--brief",
        "-y",
    ])
    job_id = job.stdout.strip()
    if not job_id.startswith("job-"):
        raise RuntimeError(f"Failed to submit job for {participant_id}. Output: {job.stdout} {job.stderr}")
    return job_id


def wait_for_jobs(job_ids: List[str], poll_seconds: int = 30) -> None:
    pending = set(job_ids)
    terminal = {"done", "failed", "terminated", "partially_failed"}
    while pending:
        finished = set()
        for job_id in list(pending):
            desc = _run(["dx", "describe", job_id, "--json"])
            state = json.loads(desc.stdout).get("state", "unknown")
            if state in terminal:
                print(f"{job_id}: {state}")
                finished.add(job_id)
                if state != "done":
                    raise RuntimeError(f"Job {job_id} finished with non-success state: {state}")
        pending -= finished
        if pending:
            print(f"Waiting on {len(pending)} jobs...")
            time.sleep(poll_seconds)


def combine_per_sample_csvs(project: str, per_sample_folder: str, output_csv_path: str) -> None:
    find_res = _run([
        "dx", "find", "data",
        "--project", project,
        "--folder", per_sample_folder,
        "--name", "*.csv",
        "--class", "file",
        "--brief",
    ])
    file_ids = [line.strip() for line in find_res.stdout.splitlines() if line.strip()]
    if not file_ids:
        raise RuntimeError(f"No per-sample CSV files found in {per_sample_folder}.")

    with tempfile.TemporaryDirectory(prefix="telomere_csvs_") as td:
        frames = []
        for file_id in file_ids:
            _run(["dx", "download", file_id, "-o", td])
        for p in Path(td).glob("*.csv"):
            frames.append(pd.read_csv(p))
        combined = pd.concat(frames, ignore_index=True)
        combined.to_csv(output_csv_path, index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Submit UKB telomere CRAM jobs in parallel and combine CSV.")
    parser.add_argument("--participant-ids", required=True, help="Text file with one participant ID per line")
    parser.add_argument("--project", required=True, help="DNANexus project ID (e.g., project-XXXX)")
    parser.add_argument(
        "--cram-folder",
        default="/Bulk/GATK and GraphTyper WGS/Whole genome GATK CRAM files and indices [500k release]/10/",
        help="Project folder containing CRAM files",
    )
    parser.add_argument(
        "--reference-fasta-file",
        default="file-G572Pj8JykJZ52BP4z6GqB21",
        help="DNANexus file ID/path for CRAM reference FASTA",
    )
    parser.add_argument(
        "--metadata-csv-file",
        default="file-J4PfY50J7pfKqx2KzZpzy5zb",
        help="DNANexus file ID/path for participant metadata CSV",
    )
    parser.add_argument(
        "--output-base-folder",
        default="/Telomere/parallel_ukb_runs",
        help="Base DNANexus folder where run outputs are written",
    )
    parser.add_argument("--instance-type", default="mem2_hdd2_x4", help="DNANexus instance type for each participant job")
    parser.add_argument("--k", type=int, default=3, help="k-repeat size (default 3)")
    parser.add_argument("--combined-local-csv", default="combined_telomere_results.csv", help="Local path for final combined CSV")
    args = parser.parse_args()

    participant_ids = read_participant_ids(args.participant_ids)
    run_tag = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_folder = f"{args.output_base_folder}/run_{run_tag}"
    jobs_folder = f"{run_folder}/jobs"
    per_sample_folder = f"{run_folder}/per_sample_csv"
    support_folder = f"{run_folder}/support"
    combined_folder = f"{run_folder}/combined"

    ensure_remote_folder(args.project, jobs_folder)
    ensure_remote_folder(args.project, per_sample_folder)
    ensure_remote_folder(args.project, support_folder)
    ensure_remote_folder(args.project, combined_folder)

    worker_script = str(Path(__file__).with_name("run_single_cram_analysis.py"))
    worker_file_id = upload_worker_script(args.project, support_folder, worker_script)

    job_ids: List[str] = []
    print(f"Submitting {len(participant_ids)} participant jobs in parallel...")
    for pid in participant_ids:
        cram_id = resolve_cram_file_id(args.project, args.cram_folder, pid)
        job_id = submit_one_job(
            destination=jobs_folder,
            cram_file_id=cram_id,
            worker_file_id=worker_file_id,
            reference_fasta_file=args.reference_fasta_file,
            metadata_csv_file=args.metadata_csv_file,
            participant_id=pid,
            per_sample_output_folder=per_sample_folder,
            k=args.k,
            instance_type=args.instance_type,
        )
        job_ids.append(job_id)
        print(f"  {pid}: {job_id}")

    print("All jobs submitted. Waiting for completion...")
    wait_for_jobs(job_ids)
    print("All participant jobs completed successfully.")

    combine_per_sample_csvs(args.project, per_sample_folder, args.combined_local_csv)
    print(f"Combined CSV written locally: {args.combined_local_csv}")

    _run([
        "dx", "upload", args.combined_local_csv,
        "--project", args.project,
        "--path", combined_folder,
        "--brief",
    ])
    print(f"Uploaded combined CSV to {args.project}:{combined_folder}")


if __name__ == "__main__":
    main()
