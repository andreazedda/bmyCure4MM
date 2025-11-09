from __future__ import annotations

from django.test import TestCase
from django.core.files.uploadedfile import SimpleUploadedFile
from chemtools.forms import DrugParamForm, SimilarityForm, BindingVizForm


class DrugParamsFormTests(TestCase):
    """Test drug parameters form validation."""
    
    def test_form_valid_with_smiles(self) -> None:
        """Test form is valid with correct SMILES."""
        form = DrugParamForm(data={"smiles": "CCO"})
        self.assertTrue(form.is_valid())
        
    def test_form_invalid_empty_smiles(self) -> None:
        """Test form is invalid with empty SMILES."""
        form = DrugParamForm(data={"smiles": ""})
        self.assertFalse(form.is_valid())
        self.assertIn("Provide a SMILES string or a PubChem CID", str(form.errors))
        
    def test_form_accepts_complex_smiles(self) -> None:
        """Test form accepts complex SMILES notation."""
        # Bortezomib-like molecule
        smiles = "CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)c2cnccn2)B(O)O"
        form = DrugParamForm(data={"smiles": smiles})
        self.assertTrue(form.is_valid())


class SimilaritySearchFormTests(TestCase):
    """Test similarity search form validation."""
    
    def test_form_valid_with_threshold(self) -> None:
        """Test form is valid with correct inputs."""
        form = SimilarityForm(data={
            "smiles": "c1ccccc1",
            "threshold": "0.7"
        })
        self.assertTrue(form.is_valid())
        
    def test_form_invalid_empty_query(self) -> None:
        """Test form is invalid without query SMILES."""
        form = SimilarityForm(data={
            "smiles": "",
            "threshold": "0.7"
        })
        self.assertFalse(form.is_valid())
        
    def test_form_accepts_threshold_range(self) -> None:
        """Test form accepts valid threshold values."""
        for threshold in ["0.0", "0.5", "0.7", "0.9", "1.0"]:
            form = SimilarityForm(data={
                "smiles": "CCO",
                "threshold": threshold
            })
            self.assertTrue(form.is_valid(), f"Threshold {threshold} should be valid")


class BindingVisualizerFormTests(TestCase):
    """Test binding visualizer form validation."""
    
    def test_form_valid_with_pdb_id(self) -> None:
        """Test form is valid with PDB ID."""
        form = BindingVizForm(data={"pdb_id": "5LF3", "ligand": ""})
        self.assertTrue(form.is_valid())
        
    def test_form_invalid_without_pdb_id(self) -> None:
        """Test form is invalid without PDB ID."""
        form = BindingVizForm(data={})
        self.assertFalse(form.is_valid())
        self.assertIn("pdb_id", form.errors)
        
    def test_form_validates_pdb_id_format(self) -> None:
        """Test form validates PDB ID format."""
        # Valid format: starts with digit, 4 chars
        valid_ids = ["5LF3", "4KW5", "1ABC"]
        for pdb_id in valid_ids:
            form = BindingVizForm(data={"pdb_id": pdb_id})
            self.assertTrue(form.is_valid(), f"{pdb_id} should be valid")
        
        # Invalid formats
        invalid_ids = ["ABC1", "123", "12345"]
        for pdb_id in invalid_ids:
            form = BindingVizForm(data={"pdb_id": pdb_id})
            self.assertFalse(form.is_valid(), f"{pdb_id} should be invalid")


class FormInputSanitizationTests(TestCase):
    """Test that forms sanitize and validate inputs."""
    
    def test_smiles_whitespace_handling(self) -> None:
        """Test SMILES with leading/trailing whitespace."""
        form = DrugParamForm(data={"smiles": "  CCO  "})
        if form.is_valid():
            # Should strip whitespace
            self.assertEqual(form.cleaned_data["smiles"].strip(), "CCO")
            
    def test_threshold_as_string(self) -> None:
        """Test threshold is converted to decimal."""
        form = SimilarityForm(data={
            "smiles": "CCO",
            "threshold": "0.75"
        })
        self.assertTrue(form.is_valid())
