conda activate CRISPRito

cd /data/CRISPRito


python -m pytest -sv tests/unit/test_SampleManager.py

python -m pytest -sv tests/unit/test_greedy_clustering.py


