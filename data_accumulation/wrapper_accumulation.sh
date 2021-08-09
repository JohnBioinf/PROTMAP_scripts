#!/bin/bash

set -e

echo "Count peptides per species"
python ./data_accumulation/count_pep_species.py

echo "Count transcriptomic reads per species"
python ./data_accumulation/count_reads_per_species.py

echo "Calculateumber of deteceted proteins"
python ./data_accumulation/detected_protein_k_psm.py

echo "Get scores of manual inspected proteins"
python ./data_accumulation/get_scores_after_manual_inspection.py

echo "Predict orf interval"
python ./data_accumulation/simple_orf_prediction.py

echo "Calculate number of ambiguous PSMs"
python ./data_accumulation/ambiguous_mapping.py

echo "Collect data for protein table"
python ./data_accumulation/detected_protein_table.py

echo "Calculate sets for venn diagramms"
python ./data_accumulation/venn_data.py
