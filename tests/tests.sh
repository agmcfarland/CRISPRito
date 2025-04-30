conda activate CRISPRito

cd /data/CRISPRito

python -m pytest --disable-warnings -sv tests

/home/ubuntu/miniconda3/envs/CRISPRito/bin/python -m pytest -sv tests/unit/test_misc.py

python -m pytest --disable-warnings -sv /data/CRISPRito/tests/unit/test_RunParameters.py

python -m pytest -sv tests/unit/test_ProcessGroup.py

python -m pytest -sv tests/unit/test_StandardCuts.py

python -m pytest -sv tests/unit/test_SampleManager.py

python -m pytest -sv tests/unit/test_greedy_clustering.py

python -m pytest -sv tests/unit/test_extract_sequences.py

python -m pytest -sv tests/unit/test_various_annotations_pyranges.py

pytest --disable-warnings -sv tests/unit/test_standard_cut_alignments.py

pytest -sv tests/unit/test_StandardCuts.py::test_multiprocessing_cut_site_alignment

pytest -sv tests/unit/test_StandardCuts.py::test_extract_in_refseq_feature


# pytest -sv tests/unit/test_standard_cut_alignments.py::test_alignment6


pytest --disable-warnings -sv tests/unit/test_standard_cut_alignments.py

pytest -sv tests/unit/test_standard_cut_alignments.py::test_alignment1
pytest -sv tests/unit/test_standard_cut_alignments.py::test_alignment2
pytest -sv tests/unit/test_standard_cut_alignments.py::test_alignment3
pytest -sv tests/unit/test_standard_cut_alignments.py::test_alignment4
pytest -sv tests/unit/test_standard_cut_alignments.py::test_alignment5

pytest -sv tests/unit/test_standard_cut_alignments.py::test_alignment6

pytest --disable-warnings -sv tests/unit/test_standard_cut_alignments.py::test_alignment7


setup_run -h


nextflow run main.nf \
  --sample_sheet /data/CRISPRito/tests/data/input_samplesheet_ptprc_reduced.csv \
  --genome /data/CRISPRito/tests/temp/new \
  --features /data/GenomicTrackRepository/data/processed/hg38/hg38.fasta.gz \
  --genes /data/GenomicTrackRepository/data/processed/hg38/ncbiRefSeqCurated_expanded.csv \
  --output_dir /data/GenomicTrackRepository/data/external/hg38/gene_names.csv \
  --overwrite_output_dir


