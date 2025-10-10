from __future__ import annotations

from django import forms

PDB_REGEX = r"^[0-9][A-Za-z0-9]{3}$"
LIGAND_REGEX = r"^[A-Za-z0-9]{1,3}$"


class DrugParamForm(forms.Form):
    smiles = forms.CharField(
        label="SMILES",
        max_length=4096,
        required=False,
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
        smiles = cleaned.get("smiles", "") or ""
        cid = cleaned.get("cid")
        if smiles:
            smiles = smiles.strip()
            if not smiles:
                smiles = ""
        if not smiles and cid is None:
            raise forms.ValidationError("Provide a SMILES string or a PubChem CID.")
        if not smiles and not cid:
            raise forms.ValidationError("Provide a SMILES string or a PubChem CID.")
        cleaned["smiles"] = smiles
        cleaned["cid"] = cid
        return cleaned


class BindingVizForm(forms.Form):
    pdb_id = forms.RegexField(
        label="PDB ID",
        max_length=4,
        regex=PDB_REGEX,
        error_messages={"invalid": "Enter a valid 4-character PDB ID (e.g., 5LF3)."},
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


class SimilarityForm(forms.Form):
    smiles = forms.CharField(
        label="SMILES",
        max_length=4096,
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
