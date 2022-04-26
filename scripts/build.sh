set -e

# Build the Docker image from the root
docker build -t phylodm_perf:latest -f scripts/Dockerfile .

# Export as singularity image
mkdir -p /tmp/singularity
docker run -v /var/run/docker.sock:/var/run/docker.sock \
  -v /tmp/singularity:/output \
  --privileged -t --rm \
  quay.io/singularity/docker2singularity:v3.8.4 \
  --name phylodm_perf phylodm_perf

scp /private/tmp/singularity/phylodm_perf.sif uqamussi@cook.ace.uq.edu.au:/srv/home/uqamussi/phylodm/phylodm.sif

# singularity run --bind /srv/home/uqamussi/phylodm/data:/io:rw phylodm.sif /io 10,50,100,500,1000,2500,5000,7500,10000,15000,20000,30000 5 1
# nice singularity run --bind /srv/home/uqamussi/phylodm/data:/io:rw phylodm.sif /io 8000 100 40