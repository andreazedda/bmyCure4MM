from __future__ import annotations

from django import forms
from django.core.exceptions import ValidationError

from mmportal.forms_mixins import BootstrapValidationMixin
from .validators import validate_smiles_rdkit, sanitize_input, EXAMPLE_SMILES

PDB_REGEX = r"^[0-9][A-Za-z0-9]{3}$"
LIGAND_REGEX = r"^[A-Za-z0-9]{1,3}$"


class DrugParamForm(BootstrapValidationMixin, forms.Form):
    smiles = forms.CharField(
        label="SMILES",
        max_length=4096,
        required=False,
        help_text=(
            f"Provide a canonical SMILES string (≤4096 characters). "
            f"Example: {EXAMPLE_SMILES['ibuprofen']} (Ibuprofen) or "
            f"{EXAMPLE_SMILES['aspirin']} (Aspirin)."
        ),
        widget=forms.TextInput(
            attrs={
                "placeholder": f"e.g., {EXAMPLE_SMILES['ethanol']} (Ethanol)",
                "maxlength": "4096",
                "class": "form-control",
            }
        ),
    )
    cid = forms.IntegerField(
        label="PubChem CID",
        required=False,
        min_value=1,
        help_text=(
            "Optional positive PubChem compound identifier. "
            "Example: 2244 (Aspirin), 3672 (Ibuprofen), 2519 (Caffeine). "
            "Find CIDs at pubchem.ncbi.nlm.nih.gov"
        ),
        widget=forms.NumberInput(
            attrs={
                "placeholder": "e.g., 2244 (Aspirin)",
                "min": "1",
                "step": "1",
                "class": "form-control",
            }
        ),
    )

    def clean(self):
        cleaned = super().clean()
        smiles = (cleaned.get("smiles") or "").strip()
        cid = cleaned.get("cid")
        if not smiles and cid is None:
            raise forms.ValidationError("Provide a SMILES string or a PubChem CID.")
        
        # Validate SMILES if provided
        if smiles:
            smiles = sanitize_input(smiles)
            is_valid, error_msg = validate_smiles_rdkit(smiles)
            if not is_valid:
                raise forms.ValidationError(f"Invalid SMILES: {error_msg}")
        
        cleaned["smiles"] = smiles if smiles else ""
        cleaned["cid"] = cid
        return cleaned


class BindingVizForm(BootstrapValidationMixin, forms.Form):
    pdb_id = forms.RegexField(
        label="PDB ID",
        max_length=4,
        regex=PDB_REGEX,
        error_messages={"invalid": "Enter a valid 4-character PDB ID (e.g., 5LF3)."},
        help_text=(
            "Four-character identifier beginning with a digit (e.g., 5LF3 for Bortezomib complex, "
            "4KW5 for Lenalidomide complex). Find PDB IDs at rcsb.org"
        ),
        widget=forms.TextInput(
            attrs={
                "placeholder": "e.g., 5LF3 (Bortezomib)",
                "pattern": PDB_REGEX,
                "class": "form-control",
                "maxlength": "4",
            }
        ),
    )
    ligand = forms.RegexField(
        label="Ligand code",
        regex=LIGAND_REGEX,
        required=False,
        error_messages={"invalid": "Ligand code must be 1–3 alphanumeric characters."},
        help_text=(
            "Optional, uppercase 1–3 character ligand code. "
            "Examples: BOR (Bortezomib), LEN (Lenalidomide), HEM (Heme). "
            "Find ligand codes in PDB structure pages."
        ),
        widget=forms.TextInput(
            attrs={
                "placeholder": "e.g., BOR (Bortezomib)",
                "pattern": LIGAND_REGEX,
                "maxlength": "3",
                "class": "form-control",
            }
        ),
    )

    def clean_pdb_id(self):
        pdb_id = self.cleaned_data["pdb_id"].upper()
        return sanitize_input(pdb_id)

    def clean_ligand(self):
        ligand = self.cleaned_data.get("ligand", "")
        if ligand:
            ligand = ligand.upper()
            ligand = sanitize_input(ligand)
        return ligand


class SimilarityForm(BootstrapValidationMixin, forms.Form):
    smiles = forms.CharField(
        label="SMILES",
        max_length=4096,
        help_text=(
            f"Canonical SMILES string (≤4096 characters). "
            f"Try: {EXAMPLE_SMILES['bortezomib'][:30]}... (Bortezomib) or "
            f"{EXAMPLE_SMILES['lenalidomide'][:30]}... (Lenalidomide)"
        ),
        widget=forms.TextInput(
            attrs={
                "placeholder": f"e.g., {EXAMPLE_SMILES['caffeine']} (Caffeine)",
                "maxlength": "4096",
                "class": "form-control",
            }
        ),
    )
    
    threshold = forms.DecimalField(
        label="Similarity Threshold",
        min_value=0.0,
        max_value=1.0,
        initial=0.7,
        decimal_places=2,
        required=False,
        help_text=(
            "Minimum similarity score (0.0-1.0). Recommended: 0.7 for similar molecules, "
            "0.9 for very similar, 0.5 for broader search."
        ),
        widget=forms.NumberInput(
            attrs={
                "placeholder": "0.7 (recommended)",
                "min": "0.0",
                "max": "1.0",
                "step": "0.05",
                "class": "form-control",
            }
        ),
    )

    def clean_smiles(self):
        value = (self.cleaned_data.get("smiles") or "").strip()
        if not value:
            raise forms.ValidationError("Enter a SMILES string.")
        
        value = sanitize_input(value)
        is_valid, error_msg = validate_smiles_rdkit(value)
        if not is_valid:
            raise forms.ValidationError(f"Invalid SMILES: {error_msg}")
        
        return value
