#!/bin/bash

# Run the pipeline inside the official Snakemake container.
# Use -v to mount current dir, -e to pass environment variables, and place the image
# name before the command to execute inside the container.
LOG_FILE="$(pwd)/data/output/logs/pipeline_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(pwd)/data/output/logs"

echo "Log salvat în: $LOG_FILE"

docker run --rm \
  --name vegas \
  --platform linux/amd64 \
  -v "$(pwd)":/vegas \
  -w /vegas \
  -e HOST_GENOME=/vegas/data/host/ \
  -e REFERENCE_GENOME=/vegas/data/reference/ \
  snakemake/snakemake \
  python src/main.py \
    -d /vegas/data/input/ \
    -o /vegas/data/output/ \
    -r /vegas/data/reference/ \
    -t /vegas/data/host/ \
    -c 1 -cc 1 --overwrite \
  2>&1 | tee "$LOG_FILE"

# Get the absolute path on the host system to show the user
HOST_OUTPUT_DIR="$(pwd)/data/output"

# Start the IGV HTTP server automatically
echo -e "\n========================================"
echo "Pornire server HTTP pentru vizualizarea IGV..."
bash scripts/start_igv_server.sh 8000 data/output

echo -e "\n========================================"
echo "Raportul IGV HTML este gata și serverul rulează!"
echo "Îl poți accesa din browser la adresa de mai jos:"
echo "http://localhost:8000/igv_report.html"
echo "(Apasă Cmd + Click pe link pentru a-l deschide)"
echo -e "========================================\n"
