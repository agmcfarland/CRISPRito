conda activate CRISPRito

cd /data/CRISPRito

python -m pytest -sv tests

python -m pytest -sv tests/unit/test_SampleManager.py

python -m pytest -sv tests/unit/test_greedy_clustering.py

python -m pytest -sv tests/unit/test_extract_sequences.py

python -m pytest -sv tests/unit/test_StandardCuts.py

python -m pytest -sv tests/unit/test_various_annotations_pyranges.py


pytest -sv tests/unit/test_StandardCuts.py::test_multiprocessing_cut_site_alignment

pytest -sv tests/unit/test_StandardCuts.py::test_extract_in_refseq_feature