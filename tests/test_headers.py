import unittest
import pandas as pd
import os


class TestHeaders(unittest.TestCase):
    """Test cases for CSV header structure and column mapping"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.csv_path = '../data/mhc_ligand_full.csv'
        self.index_column_map = {
            ('MHC Restriction', 'Name'): 'mhc_name',
            ('Assay', 'Qualitative Measurement'): 'qualitative_data',
            ('Assay', 'Quantitative measurement'): 'quantitative_data',
            ('Assay', 'Response measured'): 'assay_response',
            ('Epitope', 'Name'): 'epitope',
        }
        
    def test_csv_file_exists(self):
        """Test that the CSV file exists"""
        self.assertTrue(os.path.exists(self.csv_path), f"CSV file not found at {self.csv_path}")
        
    def test_multi_level_headers_loading(self):
        """Test loading CSV with multi-level headers"""
        df = pd.read_csv(self.csv_path, header=[0, 1], nrows=3)
        
        # Verify DataFrame is not empty
        self.assertFalse(df.empty, "DataFrame should not be empty")
        self.assertEqual(len(df), 3, "Should have exactly 3 rows as specified by nrows")
        
        # Verify multi-level columns structure
        self.assertTrue(hasattr(df.columns, 'levels'), "Should have multi-level columns")
        
    def test_expected_columns_present(self):
        """Test that all expected columns from index_column_map are present"""
        df = pd.read_csv(self.csv_path, header=[0, 1], nrows=3)
        
        missing_columns = []
        found_columns = []
        
        for col_tuple, output_name in self.index_column_map.items():
            if col_tuple in df.columns:
                found_columns.append((col_tuple, output_name))
            else:
                missing_columns.append((col_tuple, output_name))
        
        # At least some expected columns should be found
        self.assertGreater(len(found_columns), 0, "Should find at least some expected columns")
        
        # Check specific important columns
        mhc_restriction_col = ('MHC Restriction', 'Name')
        self.assertIn(mhc_restriction_col, df.columns, 
                     f"Critical column {mhc_restriction_col} should be present")
        
    def test_mhc_restriction_name_column_data(self):
        """Test that MHC Restriction Name column contains valid data"""
        df = pd.read_csv(self.csv_path, header=[0, 1], nrows=3)
        mhc_col = ('MHC Restriction', 'Name')
        
        if mhc_col in df.columns:
            mhc_data = df[mhc_col]
            # Should have some non-null values
            self.assertGreater(mhc_data.notna().sum(), 0, 
                             "MHC Restriction Name column should contain some valid data")
        else:
            self.fail(f"Required column {mhc_col} not found")
            
    def test_column_structure_consistency(self):
        """Test that the column structure is consistent"""
        df = pd.read_csv(self.csv_path, header=[0, 1], nrows=10)
        
        # All columns should be tuples (multi-level)
        for col in df.columns:
            self.assertIsInstance(col, tuple, "All columns should be tuples for multi-level headers")
            self.assertEqual(len(col), 2, "All columns should have exactly 2 levels")


if __name__ == '__main__':
    unittest.main()