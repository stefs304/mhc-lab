import unittest
import pandas as pd
import tempfile
import os
from mhc_lab.iedb.parsers import IedbMhcDataParser, IedbMhcNamesParser


class TestParser(unittest.TestCase):
    """Test cases for IEDB data parser functionality"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.test_data = """Assay,Assay,Assay,Assay,Assay,Epitope,MHC Restriction
"Qualitative Measurement","Quantitative measurement","Response measured","Units","IRI",Name,Name
Positive-High,100.5,"qualitative binding",nM,http://example.com/1,KLEDLERDL,HLA-A*02:01
Positive,75.2,"quantitative binding",μM,http://example.com/2,LITGRLQSL,HLA-A2
Negative,,,"",http://example.com/3,TRVAFAGL,H2-Kb"""
        
        # Create temporary file
        self.test_file = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
        self.test_file.write(self.test_data)
        self.test_file.close()
        
        # Create parser instances
        self.names_parser = IedbMhcNamesParser()
        self.data_parser = IedbMhcDataParser(self.names_parser)
        
    def tearDown(self):
        """Clean up test fixtures"""
        if os.path.exists(self.test_file.name):
            os.unlink(self.test_file.name)
            
    def test_parse_function_execution(self):
        """Test that the parse function executes without errors"""
        try:
            self.data_parser.parse(self.test_file.name)
        except Exception as e:
            self.fail(f"Parse function raised an exception: {e}")
            
    def test_parsed_data_shape(self):
        """Test that parsed data has expected shape"""
        self.data_parser.parse(self.test_file.name)
        
        # Should have 3 rows from test data
        self.assertEqual(len(self.data_parser.data), 3)
        self.assertGreater(len(self.data_parser.data.columns), 0)
        
    def test_parsed_data_columns(self):
        """Test that parsed data has expected column structure"""
        self.data_parser.parse(self.test_file.name)
        
        expected_columns = ['mhc_name', 'mhc_class', 'qualitative_data', 
                          'quantitative_data', 'assay_response', 'species', 'epitope']
        actual_columns = list(self.data_parser.data.columns)
        
        self.assertEqual(actual_columns, expected_columns)
        
    def test_parsed_data_content(self):
        """Test that parsed data contains expected content"""
        self.data_parser.parse(self.test_file.name)
        
        # Test epitope data preservation
        expected_epitopes = ['KLEDLERDL', 'LITGRLQSL', 'TRVAFAGL']
        actual_epitopes = self.data_parser.data['epitope'].tolist()
        self.assertEqual(actual_epitopes, expected_epitopes)
        
        # Test qualitative data preservation
        expected_qual = ['Positive-High', 'Positive', 'Negative']
        actual_qual = self.data_parser.data['qualitative_data'].tolist()
        self.assertEqual(actual_qual, expected_qual)
        
    def test_parsed_data_types(self):
        """Test that parsed data has appropriate data types"""
        self.data_parser.parse(self.test_file.name)
        
        # Check that we have a DataFrame
        self.assertIsInstance(self.data_parser.data, pd.DataFrame)
        
        # Check that string columns are object type
        string_columns = ['mhc_name', 'qualitative_data', 'assay_response', 'species', 'epitope']
        for col in string_columns:
            if col in self.data_parser.data.columns:
                self.assertEqual(self.data_parser.data[col].dtype, 'object')
                
    def test_quantitative_data_handling(self):
        """Test handling of quantitative data column"""
        self.data_parser.parse(self.test_file.name)
        
        quant_data = self.data_parser.data['quantitative_data']
        
        # Should have 2 non-null numeric values and 1 null value
        non_null_count = quant_data.notna().sum()
        self.assertEqual(non_null_count, 2)
        
        # Check specific values
        expected_values = ['100.5', '75.2']  # Values as strings from CSV
        actual_values = [str(val) for val in quant_data.dropna().tolist()]
        self.assertEqual(sorted(actual_values), sorted(expected_values))
        
    def test_mhc_name_handling(self):
        """Test handling of MHC name column"""
        self.data_parser.parse(self.test_file.name)
        
        mhc_names = self.data_parser.data['mhc_name'].tolist()
        
        # Should preserve original MHC names from input
        # Note: Exact names depend on normalization logic in parser
        self.assertEqual(len(mhc_names), 3)
        for name in mhc_names:
            self.assertIsInstance(name, str)
            self.assertNotEqual(name.strip(), '')
            
    def test_empty_values_handling(self):
        """Test handling of empty/null values in data"""
        self.data_parser.parse(self.test_file.name)
        
        # Check that the parser handles empty quantitative values gracefully
        # Row 3 has empty quantitative measurement
        row_3 = self.data_parser.data.iloc[2]
        self.assertEqual(row_3['qualitative_data'], 'Negative')
        
        # The quantitative value should be handled appropriately (null/NaN)
        self.assertTrue(pd.isna(row_3['quantitative_data']) or row_3['quantitative_data'] == '')


if __name__ == '__main__':
    unittest.main()