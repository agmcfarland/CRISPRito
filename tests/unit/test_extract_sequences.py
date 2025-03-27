import pytest
import gzip
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
from CRISPRito.Utils import genome_to_dict_memoryview, retrieve_genome_slices_memoryview

@pytest.fixture
def sample_fasta_gz():
    """Creates a temporary gzipped FASTA file for testing."""
    fasta_content = '>chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>chr2\nTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA'
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fa.gz") as tmp_fasta:
        with gzip.open(tmp_fasta.name, "wt") as f:
            f.write(fasta_content)
        return tmp_fasta.name  # Return the file path

def test_genome_to_dict_memoryview(sample_fasta_gz):
    """
    Test genome_to_dict_memoryview loads sequences as memoryviews.
    pytest -sv tests/unit/test_extract_sequences.py::test_genome_to_dict_memoryview
    """
    print(sample_fasta_gz)
    genome_dict = genome_to_dict_memoryview(sample_fasta_gz)
    assert isinstance(genome_dict, dict)
    assert "chr1" in genome_dict
    assert "chr2" in genome_dict
    assert isinstance(genome_dict["chr1"], memoryview)
    
    # Convert memoryview back to string and check contents
    assert genome_dict["chr1"].tobytes().decode().startswith("ACGTACGT")
    assert genome_dict["chr2"].tobytes().decode().startswith("TTTTGGGG")

def test_retrieve_genome_slices_memoryview():
    """
    Test retrieve_genome_slices_memoryview correctly extracts sequence slices.
    pytest -sv tests/unit/test_extract_sequences.py::test_retrieve_genome_slices_memoryview
    """
    sequence = memoryview(b"ACGTACGTACGTACGTACGTACGTACGTACGT")  # 32 bp sequence
    positions = [10, 20]  # Extract slices around these positions
    flank_size = 5  # 5 bp upstream and downstream

    slices = retrieve_genome_slices_memoryview(sequence, positions, flank_size)

    assert isinstance(slices, dict)
    assert 10 in slices
    assert 20 in slices
    print(slices[10])
    # Check actual extracted sequences
    assert slices[10] == "CGTACGTACG"
    # assert slices[20] == "ACGTACGTG"

