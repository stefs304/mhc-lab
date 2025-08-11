import unittest
import tempfile
import os
import pandas as pd
from mhc_lab.iedb.parsers import IedbMhcNamesParser, IedbMhcDataParser, IedbFilter


class TestIedbMhcNamesParser(unittest.TestCase):
    
    def setUp(self):
        """Set up test data"""
        self.test_xml_content = """<?xml version="1.0" encoding="UTF-8"?>
<MhcAlleleNameList xsi:noNamespaceSchemaLocation="http://www.iedb.org/schema/MhcAlleleNameList.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
    <MhcAlleleName>
        <MhcAlleleRestrictionId>1</MhcAlleleRestrictionId>
        <DisplayedRestriction>HLA-C*07:01</DisplayedRestriction>
        <Synonyms>HLA-C07:01|HLA-C0701|HLA-C*0701|HLA-Cw*0701|HLA-Cw*07011</Synonyms>
        <Includes>HLA-Cw*0701|HLA-Cw*07011</Includes>
        <RestrictionLevel>complete molecule</RestrictionLevel>
        <Organism>Homo sapiens (human)</Organism>
        <OrganismNcbiTaxId>9606</OrganismNcbiTaxId>
        <Class>I</Class>
        <Locus>C</Locus>
        <Molecule>HLA-C*07:01</Molecule>
        <Chain1Name>HLA-C*07:01</Chain1Name>
        <Chain2Name>Beta-2-microglobulin</Chain2Name>
    </MhcAlleleName>
    <MhcAlleleName>
        <MhcAlleleRestrictionId>2</MhcAlleleRestrictionId>
        <DisplayedRestriction>HLA-B*27:05</DisplayedRestriction>
        <Synonyms>HLA-B27:05|HLA-B2705|HLA-B*2705</Synonyms>
        <RestrictionLevel>complete molecule</RestrictionLevel>
        <Organism>Homo sapiens (human)</Organism>
        <OrganismNcbiTaxId>9606</OrganismNcbiTaxId>
        <Class>I</Class>
        <Locus>B</Locus>
        <Molecule>HLA-B*27:05</Molecule>
        <Chain1Name>HLA-B*27:05</Chain1Name>
        <Chain2Name>Beta-2-microglobulin</Chain2Name>
    </MhcAlleleName>
    <MhcAlleleName>
        <MhcAlleleRestrictionId>3</MhcAlleleRestrictionId>
        <DisplayedRestriction>H-2-Kb</DisplayedRestriction>
        <RestrictionLevel>complete molecule</RestrictionLevel>
        <Organism>Mus musculus (house mouse)</Organism>
        <OrganismNcbiTaxId>10090</OrganismNcbiTaxId>
        <Class>I</Class>
        <Locus>K</Locus>
        <Molecule>H-2-Kb</Molecule>
        <Chain1Name>H-2-Kb</Chain1Name>
        <Chain2Name>Beta-2-microglobulin</Chain2Name>
    </MhcAlleleName>
</MhcAlleleNameList>"""
        
        # Create temporary XML file
        self.temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False)
        self.temp_file.write(self.test_xml_content)
        self.temp_file.close()
        
    def tearDown(self):
        """Clean up temporary file"""
        os.unlink(self.temp_file.name)
        
    def test_parse_xml_file(self):
        """Test that parse() method correctly parses XML file"""
        parser = IedbMhcNamesParser()
        parser.parse(self.temp_file.name)
        
        # Check that names dictionary is populated
        self.assertEqual(len(parser.names), 3)
        self.assertIn('1', parser.names)
        self.assertIn('2', parser.names)
        self.assertIn('3', parser.names)
        
    def test_parse_creates_correct_objects(self):
        """Test that parsed objects have correct attributes"""
        parser = IedbMhcNamesParser()
        parser.parse(self.temp_file.name)
        
        # Test first allele
        allele1 = parser.names['1']
        self.assertEqual(allele1.name, 'HLA-C*07:01')
        self.assertEqual(allele1.mhc_class, 'I')
        self.assertEqual(allele1.chain_1_name, 'HLA-C*07:01')
        self.assertEqual(allele1.chain_2_name, 'Beta-2-microglobulin')
        self.assertEqual(allele1.organism, 'Homo sapiens (human)')
        self.assertEqual(allele1.organism_id, '9606')
        self.assertEqual(allele1.restriction_level, 'complete molecule')
        self.assertEqual(allele1.aliases, ['HLA-C07:01', 'HLA-C0701', 'HLA-C*0701', 'HLA-Cw*0701', 'HLA-Cw*07011'])
        
        # Test third allele (no synonyms)
        allele3 = parser.names['3']
        self.assertEqual(allele3.name, 'H-2-Kb')
        self.assertEqual(allele3.mhc_class, 'I')
        self.assertEqual(allele3.organism, 'Mus musculus (house mouse)')
        self.assertEqual(allele3.organism_id, '10090')
        self.assertEqual(allele3.aliases, [])  # No synonyms in this entry
        
    def test_parse_file_not_found(self):
        """Test that FileNotFoundError is raised for non-existent file"""
        parser = IedbMhcNamesParser()
        with self.assertRaises(FileNotFoundError):
            parser.parse('/non/existent/file.xml')
            
    def test_find_one_by_name(self):
        """Test finding MHC name by exact name match"""
        parser = IedbMhcNamesParser()
        parser.parse(self.temp_file.name)
        
        filter_ = IedbFilter(mhc_name='HLA-C*07:01')
        result = parser.find_one(filter_)
        
        self.assertIsNotNone(result)
        self.assertEqual(result.name, 'HLA-C*07:01')
        self.assertEqual(result.mhc_class, 'I')
        
    def test_find_one_by_alias(self):
        """Test finding MHC name by alias"""
        parser = IedbMhcNamesParser()
        parser.parse(self.temp_file.name)
        
        filter_ = IedbFilter(mhc_name='HLA-C07:01')  # This is an alias
        result = parser.find_one(filter_)
        
        self.assertIsNotNone(result)
        self.assertEqual(result.name, 'HLA-C*07:01')  # Should return main name
        
    def test_find_one_by_class(self):
        """Test finding MHC name by class"""
        parser = IedbMhcNamesParser()
        parser.parse(self.temp_file.name)
        
        filter_ = IedbFilter(mhc_class='I')
        # This should raise ValueError because multiple entries have class I
        with self.assertRaises(ValueError):
            parser.find_one(filter_)
            
    def test_find_one_by_organism(self):
        """Test finding MHC name by organism"""
        parser = IedbMhcNamesParser()
        parser.parse(self.temp_file.name)
        
        filter_ = IedbFilter(organism='Mus musculus (house mouse)')
        result = parser.find_one(filter_)
        
        self.assertIsNotNone(result)
        self.assertEqual(result.name, 'H-2-Kb')
        
    def test_find_one_not_found(self):
        """Test finding MHC name that doesn't exist"""
        parser = IedbMhcNamesParser()
        parser.parse(self.temp_file.name)
        
        filter_ = IedbFilter(mhc_name='NonExistentAllele')
        result = parser.find_one(filter_)
        
        self.assertIsNone(result)
        
    def test_find_all(self):
        """Test finding all MHC names matching a filter"""
        parser = IedbMhcNamesParser()
        parser.parse(self.temp_file.name)
        
        filter_ = IedbFilter(mhc_class='I')
        results = list(parser.find_all(filter_))
        
        self.assertEqual(len(results), 3)  # All three entries have class I
        
    def test_find_all_no_filter(self):
        """Test finding all MHC names with no filter"""
        parser = IedbMhcNamesParser()
        parser.parse(self.temp_file.name)
        
        results = list(parser.find_all(None))
        
        self.assertEqual(len(results), 3)


class TestIedbMhcDataParser(unittest.TestCase):
    
    def setUp(self):
        """Set up test data for IedbMhcDataParser tests"""
        # Create a mock IedbMhcNamesParser with test data
        self.names_parser = IedbMhcNamesParser()
        
        # Create test XML content for names parser
        test_xml_content = """<?xml version="1.0" encoding="UTF-8"?>
<MhcAlleleNameList xsi:noNamespaceSchemaLocation="http://www.iedb.org/schema/MhcAlleleNameList.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
    <MhcAlleleName>
        <MhcAlleleRestrictionId>1</MhcAlleleRestrictionId>
        <DisplayedRestriction>HLA-A*02:01</DisplayedRestriction>
        <RestrictionLevel>complete molecule</RestrictionLevel>
        <Organism>Homo sapiens (human)</Organism>
        <OrganismNcbiTaxId>9606</OrganismNcbiTaxId>
        <Class>I</Class>
        <Chain1Name>HLA-A*02:01</Chain1Name>
        <Chain2Name>Beta-2-microglobulin</Chain2Name>
    </MhcAlleleName>
    <MhcAlleleName>
        <MhcAlleleRestrictionId>2</MhcAlleleRestrictionId>
        <DisplayedRestriction>HLA-DRB1*01:01</DisplayedRestriction>
        <RestrictionLevel>complete molecule</RestrictionLevel>
        <Organism>Homo sapiens (human)</Organism>
        <OrganismNcbiTaxId>9606</OrganismNcbiTaxId>
        <Class>II</Class>
        <Chain1Name>HLA-DRA*01:01</Chain1Name>
        <Chain2Name>HLA-DRB1*01:01</Chain2Name>
    </MhcAlleleName>
    <MhcAlleleName>
        <MhcAlleleRestrictionId>3</MhcAlleleRestrictionId>
        <DisplayedRestriction>H-2-Kb</DisplayedRestriction>
        <RestrictionLevel>complete molecule</RestrictionLevel>
        <Organism>Mus musculus (house mouse)</Organism>
        <OrganismNcbiTaxId>10090</OrganismNcbiTaxId>
        <Class>I</Class>
        <Chain1Name>H-2-Kb</Chain1Name>
        <Chain2Name>Beta-2-microglobulin</Chain2Name>
    </MhcAlleleName>
</MhcAlleleNameList>"""
        
        # Create temporary XML file and parse it
        self.temp_xml_file = tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False)
        self.temp_xml_file.write(test_xml_content)
        self.temp_xml_file.close()
        self.names_parser.parse(self.temp_xml_file.name)
        
        # Create data parser with the names parser
        self.data_parser = IedbMhcDataParser(self.names_parser)
        
        # Create test DataFrame directly (simulating parsed data)
        test_data = {
            'mhc_name': ['HLA-A*02:01', 'HLA-DRB1*01:01', 'H-2-Kb', 'HLA-A*02:01', 'Unknown-MHC'],
            'mhc_class': ['I', 'II', 'I', 'I', None],
            'qualitative_data': ['Positive', None, 'Negative', '', 'Positive'],
            'quantitative_data': [None, '50.5', '', '125.0', None],
            'assay_response': ['T cell', 'B cell', 'T cell', 'T cell', 'T cell'],
            'species': ['Homo sapiens (human)', 'Homo sapiens (human)', 'Mus musculus (house mouse)', 'Homo sapiens (human)', None],
            'epitope': ['PEPTIDE1', 'PEPTIDE2', 'PEPTIDE3', 'PEPTIDE4', 'PEPTIDE5']
        }
        self.data_parser.data = pd.DataFrame(test_data)
        
    def tearDown(self):
        """Clean up temporary files"""
        os.unlink(self.temp_xml_file.name)
        
    def test_filter_empty_data(self):
        """Test filtering when data is None or empty"""
        empty_parser = IedbMhcDataParser(self.names_parser)
        filter_ = IedbFilter(mhc_name='HLA-A*02:01')
        result = empty_parser.filter(filter_)
        
        self.assertTrue(result.empty)
        self.assertEqual(list(result.columns), self.data_parser.output_data_columns)
        
    def test_filter_by_mhc_name(self):
        """Test filtering by MHC name"""
        filter_ = IedbFilter(mhc_name='HLA-A*02:01')
        result = self.data_parser.filter(filter_)
        
        self.assertEqual(len(result), 2)  # Two rows with HLA-A*02:01
        self.assertTrue(all(result['mhc_name'] == 'HLA-A*02:01'))
        
    def test_filter_by_mhc_name_case_insensitive(self):
        """Test filtering by MHC name is case-insensitive"""
        filter_ = IedbFilter(mhc_name='hla-a*02:01')
        result = self.data_parser.filter(filter_)
        
        self.assertEqual(len(result), 2)  # Two rows with HLA-A*02:01
        self.assertTrue(all(result['mhc_name'] == 'HLA-A*02:01'))
        
    def test_filter_by_mhc_class(self):
        """Test filtering by MHC class"""
        filter_ = IedbFilter(mhc_class='I')
        result = self.data_parser.filter(filter_)
        
        self.assertEqual(len(result), 3)  # Three rows with class I
        self.assertTrue(all(result['mhc_class'] == 'I'))
        
    def test_filter_by_mhc_class_case_insensitive(self):
        """Test filtering by MHC class is case-insensitive"""
        filter_ = IedbFilter(mhc_class='ii')
        result = self.data_parser.filter(filter_)
        
        self.assertEqual(len(result), 1)  # One row with class II
        self.assertTrue(all(result['mhc_class'] == 'II'))
        
    def test_filter_by_organism(self):
        """Test filtering by organism"""
        filter_ = IedbFilter(organism='Homo sapiens (human)')
        result = self.data_parser.filter(filter_)
        
        self.assertEqual(len(result), 3)  # Three human entries
        self.assertTrue(all(result['species'] == 'Homo sapiens (human)'))
        
    def test_filter_by_organism_case_insensitive(self):
        """Test filtering by organism is case-insensitive"""
        filter_ = IedbFilter(organism='mus musculus (house mouse)')
        result = self.data_parser.filter(filter_)
        
        self.assertEqual(len(result), 1)  # One mouse entry
        self.assertTrue(all(result['species'] == 'Mus musculus (house mouse)'))
        
    def test_filter_by_qualitative_data(self):
        """Test filtering by qualitative data presence"""
        filter_ = IedbFilter(qualitative_data=True)
        result = self.data_parser.filter(filter_)
        
        # Should return rows where qualitative_data is not None, empty string, or 'null'
        self.assertEqual(len(result), 3)  # 'Positive', 'Negative', 'Positive'
        valid_qualitative = result['qualitative_data'].tolist()
        self.assertIn('Positive', valid_qualitative)
        self.assertIn('Negative', valid_qualitative)
        
    def test_filter_by_quantitative_data(self):
        """Test filtering by quantitative data presence"""
        filter_ = IedbFilter(quantitative_data=True)
        result = self.data_parser.filter(filter_)
        
        # Should return rows where quantitative_data is not None, empty string, or 'null'
        self.assertEqual(len(result), 2)  # '50.5', '125.0'
        valid_quantitative = result['quantitative_data'].tolist()
        self.assertIn('50.5', valid_quantitative)
        self.assertIn('125.0', valid_quantitative)
        
    def test_filter_combined_criteria(self):
        """Test filtering with multiple criteria"""
        filter_ = IedbFilter(mhc_class='I', organism='Homo sapiens (human)')
        result = self.data_parser.filter(filter_)
        
        self.assertEqual(len(result), 2)  # Two human class I entries
        self.assertTrue(all(result['mhc_class'] == 'I'))
        self.assertTrue(all(result['species'] == 'Homo sapiens (human)'))
        
    def test_filter_chain_filtering(self):
        """Test filtering by chain information"""
        filter_ = IedbFilter(chain_1='HLA-A*02:01')
        result = self.data_parser.filter(filter_)
        
        # Should return rows where MHC names have the specified chain_1
        # Based on our test XML, HLA-A*02:01 has Chain1Name='HLA-A*02:01'
        self.assertEqual(len(result), 2)  # Two HLA-A*02:01 entries
        self.assertTrue(all(result['mhc_name'] == 'HLA-A*02:01'))
        
    def test_filter_chain_any_filtering(self):
        """Test filtering by any chain information"""
        filter_ = IedbFilter(chain_any='Beta-2-microglobulin')
        result = self.data_parser.filter(filter_)
        
        # Should return rows where MHC names have Beta-2-microglobulin as any chain
        # Both HLA-A*02:01 and H-2-Kb have this as chain_2
        expected_mhc_names = ['HLA-A*02:01', 'H-2-Kb']
        result_mhc_names = result['mhc_name'].unique().tolist()
        self.assertEqual(len(result), 3)  # Two HLA-A*02:01 + one H-2-Kb
        for name in result_mhc_names:
            self.assertIn(name, expected_mhc_names)
            
    def test_filter_no_matches(self):
        """Test filtering that returns no matches"""
        filter_ = IedbFilter(mhc_name='NonExistentAllele')
        result = self.data_parser.filter(filter_)
        
        self.assertTrue(result.empty)
        self.assertEqual(list(result.columns), self.data_parser.output_data_columns)
        
    def test_filter_returns_correct_columns(self):
        """Test that filter returns DataFrame with correct columns"""
        filter_ = IedbFilter(mhc_class='I')
        result = self.data_parser.filter(filter_)
        
        self.assertEqual(list(result.columns), self.data_parser.output_data_columns)
        expected_columns = ['mhc_name', 'mhc_class', 'qualitative_data', 'quantitative_data', 'assay_response', 'species', 'epitope']
        self.assertEqual(list(result.columns), expected_columns)


if __name__ == '__main__':
    unittest.main()
