from __future__ import annotations

from django import forms

from mmportal.forms_mixins import BootstrapValidationMixin

PDB_REGEX = r"^[0-9][A-Za-z0-9]{3}$"
LIGAND_REGEX = r"^[A-Za-z0-9]{1,3}$"


class DrugParamForm(BootstrapValidationMixin, forms.Form):
    smiles = forms.CharField(
        label="SMILES",
        max_length=4096,
        required=False,
        help_text="Provide a canonical SMILES string (≤4096 characters).",
        widget=forms.TextInput(
            attrs={
                "placeholder": "Canonical SMILES (optional)",
                "maxlength": "4096",
                "class": "form-control",
            }
        ),
    )
    cid = forms.IntegerField(
        label="PubChem CID",
        required=False,
        min_value=1,
        help_text="Optional positive PubChem compound identifier.",
        widget=forms.NumberInput(
            attrs={
                "placeholder": "e.g., 2244",
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
        cleaned["smiles"] = smiles if smiles else ""
        cleaned["cid"] = cid
        return cleaned


class BindingVizForm(BootstrapValidationMixin, forms.Form):
    pdb_id = forms.RegexField(
        label="PDB ID",
        max_length=4,
        regex=PDB_REGEX,
        error_messages={"invalid": "Enter a valid 4-character PDB ID (e.g., 5LF3)."},
        help_text="Four-character identifier beginning with a digit (e.g., 5LF3).",
        widget=forms.TextInput(
            attrs={
                "placeholder": "e.g., 5LF3",
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
        help_text="Optional, uppercase 1–3 character ligand code (e.g., BOR).",
        widget=forms.TextInput(
            attrs={
                "placeholder": "Optional (e.g., BOR)",
                "pattern": LIGAND_REGEX,
                "maxlength": "3",
                "class": "form-control",
            }
        ),
    )

    def clean_pdb_id(self):
        return self.cleaned_data["pdb_id"].upper()

    def clean_ligand(self):
        ligand = self.cleaned_data.get("ligand", "")
        return ligand.upper() if ligand else ""


class SimilarityForm(BootstrapValidationMixin, forms.Form):
    smiles = forms.CharField(
        label="SMILES",
        max_length=4096,
        help_text="Canonical SMILES string (≤4096 characters).",
        widget=forms.TextInput(
            attrs={
                "placeholder": "Canonical SMILES",
                "maxlength": "4096",
                "class": "form-control",
            }
        ),
    )

    def clean_smiles(self):
        value = (self.cleaned_data.get("smiles") or "").strip()
        if not value:
            raise forms.ValidationError("Enter a SMILES string.")
        return value
