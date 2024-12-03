#!/usr/bin/bash


alevin-fry generate-permit-list -d fw -k -i pbmc_splici_map -o pbmc_splici_quant

alevin-fry collate -t 16 -i pbmc_splici_quant -r pbmc_splici_map

alevin-fry quant -t 16 -i pbmc_splici_quant -o pbmc_splici_quant_res --tg-map hs_GRCh38_113_splici_python/splici_fl146_t2g_3col.tsv --resolution cr-like --use-mtx

