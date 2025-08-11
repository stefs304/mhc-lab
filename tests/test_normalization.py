import unittest
import pandas as pd
import tempfile
import os
from mhc_lab.iedb.parsers import IedbMhcDataParser, IedbMhcNamesParser, IedbMhcName


class MockMhcNamesParser(IedbMhcNamesParser):
    """Mock MHC names parser for testing normalization"""
    def __init__(self):
        super().__init__()
        # Add test MHC names with aliases
        self.names = {
            'HLA-A*02:01': IedbMhcName(
                name='HLA-A*02:01',
                mhc_class='I',
                chain_1_name='HLA-A*02:01',
                chain_2_name=None,
                aliases=['HLA-A2', 'A*02:01', 'A2'],
                organism='Human',
                organism_id='9606',
                restriction_level='Allele'
            ),
            'H2-Kb': IedbMhcName(
                name='H2-Kb',
                mhc_class='I',
                chain_1_name='H2-Kb',
                chain_2_name=None,
                aliases=['Kb', 'H-2Kb'],
                organism='Mouse',
                organism_id='10090',
                restriction_level='Allele'
            )
        }


class TestNormalization(unittest.TestCase):
    """Test cases for MHC name normalization functionality"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.test_data = """Assay,Assay,Assay,Assay,Assay,Epitope,MHC Restriction
"Qualitative Measurement","Quantitative measurement","Response measured","Units","IRI",Name,Name
Positive-High,100.5,"qualitative binding",nM,http://example.com/1,KLEDLERDL,HLA-A*02:01
Positive,75.2,"quantitative binding",μM,http://example.com/2,LITGRLQSL,HLA-A2
Negative,,,"",http://example.com/3,TRVAFAGL,H-2Kb
Positive,50.0,"binding",nM,http://example.com/4,UNKNOWN,UnknownMhc"""
        
        # Create temporary file
        self.test_file = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
        self.test_file.write(self.test_data)
        self.test_file.close()
        
        # Create parser instances
        self.names_parser = MockMhcNamesParser()
        self.data_parser = IedbMhcDataParser(self.names_parser)
        
    def tearDown(self):
        """Clean up test fixtures"""
        if os.path.exists(self.test_file.name):
            os.unlink(self.test_file.name)
            
    def test_mhc_name_normalization_direct_match(self):
        """Test normalization of MHC names that match directly"""
        self.data_parser.parse(self.test_file.name)
        
        # Check direct match: HLA-A*02:01 should remain as is
        direct_match_rows = self.data_parser.data[
            self.data_parser.data['epitope'] == 'KLEDLERDL'
        ]
        self.assertEqual(len(direct_match_rows), 1)
        row = direct_match_rows.iloc[0]
        self.assertEqual(row['mhc_name'], 'HLA-A*02:01')
        self.assertEqual(row['mhc_class'], 'I')
        self.assertEqual(row['species'], 'Human')
        
    def test_mhc_name_normalization_alias_match(self):
        """Test normalization of MHC names through aliases"""
        self.data_parser.parse(self.test_file.name)
        
        # Check alias match: HLA-A2 should be normalized to HLA-A*02:01
        alias_match_rows = self.data_parser.data[
            self.data_parser.data['epitope'] == 'LITGRLQSL'
        ]
        self.assertEqual(len(alias_match_rows), 1)
        row = alias_match_rows.iloc[0]
        self.assertEqual(row['mhc_name'], 'HLA-A*02:01')
        self.assertEqual(row['mhc_class'], 'I')
        self.assertEqual(row['species'], 'Human')
        
    def test_mhc_name_normalization_mouse_alias(self):
        """Test normalization of mouse MHC names through aliases"""
        self.data_parser.parse(self.test_file.name)
        
        # Check mouse alias: H-2Kb should be normalized to H2-Kb
        mouse_rows = self.data_parser.data[
            self.data_parser.data['epitope'] == 'TRVAFAGL'
        ]
        self.assertEqual(len(mouse_rows), 1)
        row = mouse_rows.iloc[0]
        self.assertEqual(row['mhc_name'], 'H2-Kb')
        self.assertEqual(row['mhc_class'], 'I')
        self.assertEqual(row['species'], 'Mouse')
        
    def test_mhc_name_normalization_no_match(self):
        """Test handling of unknown MHC names"""
        self.data_parser.parse(self.test_file.name)
        
        # Check unknown MHC: UnknownMhc should remain unchanged with None values
        unknown_rows = self.data_parser.data[
            self.data_parser.data['epitope'] == 'UNKNOWN'
        ]
        self.assertEqual(len(unknown_rows), 1)
        row = unknown_rows.iloc[0]
        self.assertEqual(row['mhc_name'], 'UnknownMhc')
        self.assertIsNone(row['mhc_class'])
        self.assertIsNone(row['species'])
        
    def test_normalization_data_integrity(self):
        """Test that normalization preserves other data integrity"""
        self.data_parser.parse(self.test_file.name)
        
        # Verify all rows are present
        self.assertEqual(len(self.data_parser.data), 4)
        
        # Verify epitope data is preserved
        expected_epitopes = ['KLEDLERDL', 'LITGRLQSL', 'TRVAFAGL', 'UNKNOWN']
        actual_epitopes = self.data_parser.data['epitope'].tolist()
        self.assertEqual(sorted(actual_epitopes), sorted(expected_epitopes))
        
        # Verify qualitative data is preserved
        qual_data = self.data_parser.data['qualitative_data'].tolist()
        expected_qual = ['Positive-High', 'Positive', 'Negative', 'Positive']
        self.assertEqual(sorted(qual_data), sorted(expected_qual))
        
    def test_mock_names_parser_setup(self):
        """Test that the mock names parser is properly configured"""
        # Test direct name lookup
        self.assertIn('HLA-A*02:01', self.names_parser.names)
        self.assertIn('H2-Kb', self.names_parser.names)
        
        # Test alias lookup
        hla_name = self.names_parser.names['HLA-A*02:01']
        self.assertIn('HLA-A2', hla_name.aliases)
        self.assertIn('A*02:01', hla_name.aliases)
        
        h2_name = self.names_parser.names['H2-Kb']
        self.assertIn('H-2Kb', h2_name.aliases)


if __name__ == '__main__':
    unittest.main()