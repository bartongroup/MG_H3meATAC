#!/bin/env bash
 snakemake  --rerun-incomplete --latency-wait 120 --jobs 100 --keep-going --cluster-config config/cluster_config.json \
  --drmaa " -cwd -j y -o cluster_logs/{rule}.log -mods l_hard mfree {cluster.mem_free} -adds l_hard local_free {cluster.local_free}  -pe smp {threads}"
