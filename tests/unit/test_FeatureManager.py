import pytest
import pandas as pd
from unittest.mock import patch, MagicMock
from CRISPRito.FeatureManager import FeatureManager  # Replace with actual import

# Sample mocked feature table
valid_feature_table = pd.DataFrame({
    'feature': ['gene', 'exon'],
    'type': ['annotation', 'presence_absence'],
    'file_path': ['gene_names.csv', 'exon_ranges.csv']
})

# Mocked data returned by pd.read_csv
valid_gene_df = pd.DataFrame({
    'chromosome': ['chr1'],
    'strand': ['+'],
    'start': [1000],
    'end': [3000],
    'annotation': ['TP53']
})

valid_exon_df = pd.DataFrame({
    'chromosome': ['chr2'],
    'strand': ['-'],
    'start': [2000],
    'end': [3000]
})


@patch('os.path.exists', return_value=True)
@patch('pandas.read_csv')
def test_feature_manager_valid(mock_read_csv, mock_exists):
    """
    pytest -sv tests/unit/test_FeatureManager.py::test_feature_manager_valid
    """
    # Order of calls matters here: gene first, then exon
    mock_read_csv.side_effect = [valid_gene_df, valid_exon_df]

    fm = FeatureManager(valid_feature_table)
    registry = fm.registry

    assert 'gene' in registry
    assert 'exon' in registry
    assert registry['gene']['type'] == 'annotation'
    assert registry['exon']['type'] == 'presence_absence'


@patch('os.path.exists', return_value=False)
def test_feature_manager_missing_file(mock_exists):
    """
    pytest -sv tests/unit/test_FeatureManager.py::test_feature_manager_missing_file
    """
    with pytest.raises(FileNotFoundError):
        FeatureManager(valid_feature_table)


def test_feature_manager_invalid_type():
    """
    pytest -sv tests/unit/test_FeatureManager.py::test_feature_manager_invalid_type
    """
    table = valid_feature_table.copy()
    table.loc[0, 'type'] = 'invalid_type'
    with patch('os.path.exists', return_value=True), \
         patch('pandas.read_csv', return_value=valid_gene_df):
        with pytest.raises(ValueError, match="Invalid type"):
            FeatureManager(table)


def test_feature_manager_duplicate_features():
    """
    pytest -sv tests/unit/test_FeatureManager.py::test_feature_manager_duplicate_features
    """
    table = valid_feature_table.copy()
    table.loc[1, 'feature'] = 'gene'  # duplicate 'gene'
    with patch('os.path.exists', return_value=True), \
         patch('pandas.read_csv', return_value=valid_gene_df):
        with pytest.raises(ValueError, match="Duplicate feature"):
            FeatureManager(table)


def test_feature_manager_missing_annotation_column():
    """
    pytest -sv tests/unit/test_FeatureManager.py::test_feature_manager_missing_annotation_column
    """
    table = valid_feature_table.copy()
    table.loc[0, 'type'] = 'annotation'

    bad_df = pd.DataFrame({
        'chromosome': ['chr1'],
        'strand': ['+'],
        'start': [1000],
        'end' : [2000]
    })

    with patch('os.path.exists', return_value=True), \
         patch('pandas.read_csv', return_value=bad_df):
        with pytest.raises(ValueError, match="Missing 'annotation' column"):
            FeatureManager(table)

def test_feature_manager_missing_standard_column():
    """
    pytest -sv tests/unit/test_FeatureManager.py::test_feature_manager_missing_standard_column
    """
    table = valid_feature_table.copy()
    table.loc[0, 'type'] = 'annotation'

    bad_df = pd.DataFrame({
        'chromosome': ['chr1'],
        'strand': ['+'],
        'start': [1000],
        # 'end' : [2000]
    })

    with patch('os.path.exists', return_value=True), \
         patch('pandas.read_csv', return_value=bad_df):
        with pytest.raises(ValueError, match="Missing required column 'end' in gene"):
            FeatureManager(table)