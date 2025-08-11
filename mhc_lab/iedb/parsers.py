import dataclasses
import time
from typing import Dict, Iterable, List, Optional, Tuple, Union
import xml.etree.ElementTree as ET
import os

import pandas as pd


@dataclasses.dataclass
class IedbFilter:
    mhc_name: str = None
    mhc_class: str = None
    chain_1: str = None
    chain_2: str = None
    chain_any: str = None
    restriction_level: str = None
    organism_id: str = None
    organism: str = None

    # exclusively data filters
    qualitative_data: bool = False
    quantitative_data: bool = False


class IedbMhcNamesParser:

    def __init__(self):
        self.names: Dict[str, 'IedbMhcName'] = {}

    def parse(self, path: str) -> None:
        """
        Parse the MHC allele names XML file and store IedbMhcName objects in memory.
        
        Args:
            path (str): Path to the XML file containing MHC allele names.
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"MHC allele names file not found: {path}")

        tree = ET.parse(path)
        root = tree.getroot()

        for allele_elem in root.findall('MhcAlleleName'):
            # Extract required fields
            restriction_id = allele_elem.find('MhcAlleleRestrictionId').text
            name = allele_elem.find('DisplayedRestriction').text
            
            # Extract optional fields with None defaults
            mhc_class = allele_elem.find('Class').text if allele_elem.find('Class') is not None else None
            chain_1_name = allele_elem.find('Chain1Name').text if allele_elem.find('Chain1Name') is not None else None
            chain_2_name = allele_elem.find('Chain2Name').text if allele_elem.find('Chain2Name') is not None else None
            organism = allele_elem.find('Organism').text if allele_elem.find('Organism') is not None else None
            organism_id = allele_elem.find('OrganismNcbiTaxId').text if allele_elem.find('OrganismNcbiTaxId') is not None else None
            restriction_level = allele_elem.find('RestrictionLevel').text if allele_elem.find('RestrictionLevel') is not None else None
            
            # Process synonyms (aliases)
            aliases = []
            synonyms_elem = allele_elem.find('Synonyms')
            if synonyms_elem is not None and synonyms_elem.text:
                aliases = synonyms_elem.text.split('|')
            
            # Create IedbMhcName object
            mhc_name_obj = IedbMhcName(
                name=name,
                mhc_class=mhc_class,
                chain_1_name=chain_1_name,
                chain_2_name=chain_2_name,
                aliases=aliases,
                organism=organism,
                organism_id=organism_id,
                restriction_level=restriction_level
            )
            
            # Store in names dictionary using restriction_id as key
            self.names[restriction_id] = mhc_name_obj

    def find_one(self, filter_: IedbFilter) -> Optional['IedbMhcName']:
        """
        Returns the single matching IedbMhcName, None if no matches, and raises ValueError if more than one match is found.

        The filtering logic applies all non-None attributes from the provided IedbFilter:
        - mhc_name: matches if equal to the IedbMhcName.name or any alias (case-insensitive).
        - mhc_class: equals IedbMhcName.mhc_class (case-insensitive).
        - chain_1: equals IedbMhcName.chain_1_name (case-insensitive).
        - chain_2: equals IedbMhcName.chain_2_name (case-insensitive).
        - chain_any: equals either chain_1_name or chain_2_name (case-insensitive).
        - restriction_level: equals IedbMhcName.restriction_level (case-insensitive).
        - organism_id: equals IedbMhcName.organism_id (string-compare).
        - organism: equals IedbMhcName.organism (case-insensitive).

        """
        matches = list(self.find_all(filter_))
        if not matches:
            return None
        if len(matches) > 1:
            raise ValueError("Multiple MHC names found for filter")
        return matches[0]

    def find_all(self, filter_: IedbFilter) -> Iterable['IedbMhcName']:
        """Finds all instances of IedbMhcName from a filter."""

        def _norm(v: Optional[str]) -> Optional[str]:
            return v.lower() if isinstance(v, str) else v

        candidates = list(self.names.values())
        if filter_ is None:
            # With no filter, return all available candidates
            return candidates

        f_mhc_name = _norm(filter_.mhc_name)
        f_class = _norm(filter_.mhc_class)
        f_chain1 = _norm(filter_.chain_1)
        f_chain2 = _norm(filter_.chain_2)
        f_chain_any = _norm(filter_.chain_any)
        f_rl = _norm(filter_.restriction_level)
        f_org = _norm(filter_.organism)
        f_org_id = filter_.organism_id  # keep as-is

        def match(item: 'IedbMhcName') -> bool:
            name = _norm(item.name)
            mhc_class = _norm(item.mhc_class)
            c1 = _norm(item.chain_1_name)
            c2 = _norm(item.chain_2_name)
            rl = _norm(item.restriction_level)
            org = _norm(item.organism)
            org_id = item.organism_id

            if f_mhc_name is not None:
                aliases = [_norm(a) for a in (item.aliases or []) if a is not None]
                if not ((name == f_mhc_name) or (f_mhc_name in aliases)):
                    return False

            if f_class is not None and mhc_class != f_class:
                return False
            if f_chain1 is not None and c1 != f_chain1:
                return False
            if f_chain2 is not None and c2 != f_chain2:
                return False
            if f_chain_any is not None and not (c1 == f_chain_any or c2 == f_chain_any):
                return False
            if f_rl is not None and rl != f_rl:
                return False
            if f_org is not None and org != f_org:
                return False
            if f_org_id is not None and org_id != f_org_id:
                return False
            return True

        return [x for x in candidates if match(x)]


class IedbMhcName:
    """
    Class representing an Iedb MHC Name.
    """

    def __init__(self, name, mhc_class, chain_1_name, chain_2_name, aliases, organism, organism_id, restriction_level):
        self.name = name
        self.mhc_class = mhc_class
        self.chain_1_name = chain_1_name
        self.chain_2_name = chain_2_name
        self.aliases = aliases
        self.organism = organism
        self.organism_id = organism_id
        self.restriction_level = restriction_level

    def __str__(self):
        return f'{self.name} - class {self.mhc_class} MHC molecule'


class IedbMhcDataParser:
    """
    Parses MHC ligand data from IEDB. Provides methods for filtering data.

    Attributes:
        index_field_map: maps data column names to filter fields.
        names: Initialized instance of IedbMhcNamesParser.
        data: MHC ligand data. Has two column indexes, represented by tuple of keys in index_field_map.

    Methods:
        parse: used to load and initialize data from file
        filter: used for filtering data
    """

    index_column_map: Dict[str, str] = {
        ('MHC Restriction', 'Name'): 'mhc_name',
        ('Assay', 'Qualitative Measurement'): 'qualitative_data',
        ('Assay', 'Quantitative measurement'): 'quantitative_data',
        ('Assay', 'Response measured'): 'assay_response',
        ('Epitope', 'Name'): 'epitope',
    }
    output_data_columns: List[str] = [
        'mhc_name',
        'mhc_class',
        'qualitative_data',
        'quantitative_data',
        'assay_response',
        'species',
        'epitope'
    ]

    @property
    def field_index_map(self) -> dict:
        """Reverse of index_field_map. Maps filter fields to data column names."""
        return {v: k for k, v in self.index_column_map.items()}

    def __init__(self, names: IedbMhcNamesParser):
        """
        Creates an instance of IedbMhcDataParser.
        :param names: Initialized and loaded instance of IedbMhcNamesParser.
        """
        self.names = names
        self.data: pd.DataFrame = None

    def parse(self, path: str) -> None:
        """
        :param path:
        :return: None. Data is stored in self.data, in a pandas dataframe.
        """
        print(f"Parsing MHC ligand data from {path}...")
        
        # Read CSV in chunks due to large size (8M+ rows)
        chunk_size = 50000
        processed_chunks = []
        
        # Read with multi-level headers (2 header rows as indicated by index_column_map)
        csv_reader = pd.read_csv(path, header=[0, 1], chunksize=chunk_size, low_memory=False)
        
        chunk_count = 0
        for chunk in csv_reader:
            chunk_count += 1
            print(f"Processing chunk {chunk_count} with {len(chunk)} rows...")
            
            # Select only the columns we need based on index_column_map
            selected_data = {}
            
            for col_tuple, output_name in self.index_column_map.items():
                if col_tuple in chunk.columns:
                    selected_data[output_name] = chunk[col_tuple]
                else:
                    # If column is missing, create empty series
                    selected_data[output_name] = pd.Series([None] * len(chunk))
            
            # Create DataFrame with renamed columns
            processed_chunk = pd.DataFrame(selected_data)
            
            # Infer MHC class and species from MHC names using the names parser
            # Also normalize MHC names by translating aliases to main names
            mhc_class_list = []
            species_list = []
            normalized_mhc_names = []
            
            for mhc_name in processed_chunk['mhc_name']:
                if pd.isna(mhc_name) or mhc_name is None:
                    mhc_class_list.append(None)
                    species_list.append(None)
                    normalized_mhc_names.append(mhc_name)  # Keep original None/NaN
                    continue
                
                try:
                    # Create filter to find MHC name object
                    filter_ = IedbFilter(mhc_name=str(mhc_name))
                    mhc_obj = self.names.find_one(filter_)
                    
                    if mhc_obj:
                        mhc_class_list.append(mhc_obj.mhc_class)
                        species_list.append(mhc_obj.organism)
                        # Normalize MHC name to the main name (translate aliases)
                        normalized_mhc_names.append(mhc_obj.name)
                    else:
                        # If no MHC object found, set to None and keep original name
                        mhc_class_list.append(None)
                        species_list.append(None)
                        normalized_mhc_names.append(mhc_name)
                        
                except (ValueError, Exception):
                    # If multiple matches or other error, set to None and keep original name
                    mhc_class_list.append(None)
                    species_list.append(None)
                    normalized_mhc_names.append(mhc_name)
            
            # Replace original MHC names with normalized ones
            processed_chunk['mhc_name'] = normalized_mhc_names
            # Add inferred MHC class and species columns
            processed_chunk['mhc_class'] = mhc_class_list
            processed_chunk['species'] = species_list
            
            # Ensure all output columns are present
            for col in self.output_data_columns:
                if col not in processed_chunk.columns:
                    processed_chunk[col] = None
            
            # Select only the output columns in the correct order
            processed_chunk = processed_chunk[self.output_data_columns]
            
            processed_chunks.append(processed_chunk)
        
        # Concatenate all processed chunks
        if processed_chunks:
            self.data = pd.concat(processed_chunks, ignore_index=True)
            print(f"Finished parsing. Total rows: {len(self.data)}")
        else:
            self.data = pd.DataFrame(columns=self.output_data_columns)
            print("No data processed - empty result.")

    def filter(self, filter_: IedbFilter) -> pd.DataFrame:
        """
        Filters data based on the provided IedbFilter filter and supported fields for querying.
        Supported fields for querying:
            mhc_name
            mhc_class
            chain_1_name
            chain_2_name
            organism
            organism_id
            qualitative_data
            quantitative_data
        Functions performs some additional filtering operations.

        :param filter_: IedbFilter
        :return: pd.DataFrame. Data output has self.output_data_columns column names.
        """
        if self.data is None or self.data.empty:
            return pd.DataFrame(columns=self.output_data_columns)
        
        # Start with all data
        filtered_data = self.data.copy()
        
        # Apply filters based on non-None filter attributes
        if filter_.mhc_name is not None:
            # Filter by MHC name (case-insensitive exact match)
            filtered_data = filtered_data[
                filtered_data['mhc_name'].str.lower() == filter_.mhc_name.lower()
            ]
        
        if filter_.mhc_class is not None:
            # Filter by MHC class (case-insensitive exact match)
            filtered_data = filtered_data[
                filtered_data['mhc_class'].str.lower() == filter_.mhc_class.lower()
            ]
        
        if filter_.organism is not None:
            # Filter by organism/species (case-insensitive exact match)
            filtered_data = filtered_data[
                filtered_data['species'].str.lower() == filter_.organism.lower()
            ]
        
        if filter_.organism_id is not None:
            # Filter by organism ID - note: we don't have this column in output_data_columns
            # This would require additional implementation if organism_id data was available
            pass
        
        # Handle chain filtering - these would require looking up MHC objects
        # Since we don't store chain info directly in the data, we'd need to lookup via names parser
        if filter_.chain_1 is not None or filter_.chain_2 is not None or filter_.chain_any is not None:
            # Filter by checking MHC names against the names parser for chain info
            valid_mhc_names = []
            for mhc_name in filtered_data['mhc_name'].dropna().unique():
                try:
                    # Create a filter for this MHC name and add chain criteria
                    mhc_filter = IedbFilter(
                        mhc_name=mhc_name,
                        chain_1=filter_.chain_1,
                        chain_2=filter_.chain_2,
                        chain_any=filter_.chain_any
                    )
                    # Check if this MHC name matches the chain criteria
                    if self.names.find_one(mhc_filter) is not None:
                        valid_mhc_names.append(mhc_name)
                except (ValueError, Exception):
                    # Skip if multiple matches or other errors
                    continue
            
            # Filter data to only include valid MHC names
            if valid_mhc_names:
                filtered_data = filtered_data[filtered_data['mhc_name'].isin(valid_mhc_names)]
            else:
                # No valid MHC names found, return empty DataFrame
                filtered_data = pd.DataFrame(columns=self.output_data_columns)
        
        # Filter by qualitative data presence
        if filter_.qualitative_data:
            filtered_data = filtered_data[
                filtered_data['qualitative_data'].notna() & 
                (filtered_data['qualitative_data'] != '') &
                (filtered_data['qualitative_data'] != 'null')
            ]
        
        # Filter by quantitative data presence
        if filter_.quantitative_data:
            filtered_data = filtered_data[
                filtered_data['quantitative_data'].notna() & 
                (filtered_data['quantitative_data'] != '') &
                (filtered_data['quantitative_data'] != 'null')
            ]
        
        return filtered_data

    @staticmethod
    def _show_details(data: pd.DataFrame, filter_: IedbFilter, t_elapsed: float) -> None:
        """
        Prints useful information about data that was retrieved by the applied filter.
        Prints information about the filter applied.
        """
        pass

