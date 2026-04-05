#!/bin/bash
# Wait for the Snakemake pipeline to finish, then stop the VM.
# Usage: run AFTER starting snakemake in the same tmux session, e.g.:
#   snakemake --cores $(nproc) --use-conda 2>&1 | tee pipeline.log ; bash auto_stop.sh

ZONE=us-central1-a
INSTANCE=splice-pipeline
LOG=/home/jin-holee/splice-neoepitope-pipeline/auto_stop.log

echo "[$(date)] Waiting for snakemake to finish..." >> "$LOG"

while pgrep -x snakemake > /dev/null; do
    sleep 60
done

echo "[$(date)] Snakemake finished. Stopping VM in 60 seconds (Ctrl+C to cancel)..." >> "$LOG"
echo "Snakemake finished. Stopping VM in 60 seconds (Ctrl+C to cancel)..."
sleep 60

echo "[$(date)] Shutting down OS (GCP will stop the instance automatically)..." >> "$LOG"
sudo shutdown -h now >> "$LOG" 2>&1
echo "[$(date)] Done." >> "$LOG"
