#!/usr/bin/bash

DATA_DIR="/work/vetmed_data/jj/projects/juanJovel/courses/scRNAseq_2024/Parent_NGSC3_DI_PBMC_fastqs"
OUTDIR="/work/vetmed_data/jj/projects/juanJovel/courses/scRNAseq_2024/scVelo"

salmon alevin -i "${OUTDIR}/hs_GRCh38_splici_fl146_idx" -p 16 -l ISR --chromium --sketch \
-1 "${DATA_DIR}/Parent_NGSC3_DI_PBMC_R1.fq.gz" \
-2 "${DATA_DIR}/Parent_NGSC3_DI_PBMC_R2.fq.gz" \
-o "${OUTDIR}/pbmc_splici_map"
