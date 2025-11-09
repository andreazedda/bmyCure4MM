"""
Validation utilities for chemtools inputs.
Provides SMILES validation, PDB format checking, and security controls.
"""
from __future__ import annotations

from typing import Tuple
from django.core.exceptions import ValidationError


def validate_smiles_basic(smiles: str) -> Tuple[bool, str]:
    """
    Basic SMILES string validation.
    Checks for:
    - Non-empty string
    - Reasonable length
    - No SQL injection patterns
    - Basic SMILES characters
    
    Returns (is_valid, error_message)
    """
    if not smiles or not smiles.strip():
        return False, "SMILES string cannot be empty"
    
    smiles = smiles.strip()
    
    # Check length (very long SMILES can cause performance issues)
    if len(smiles) > 4096:
        return False, "SMILES string too long (max 4096 characters)"
    
    # Basic security check - reject SQL keywords
    dangerous_patterns = ["DROP", "DELETE", "INSERT", "UPDATE", "ALTER", "--", ";"]
    smiles_upper = smiles.upper()
    for pattern in dangerous_patterns:
        if pattern in smiles_upper:
            return False, f"Invalid characters detected: {pattern}"
    
    # Basic SMILES character check (not exhaustive, but catches obvious errors)
    # Valid SMILES uses: C, N, O, S, P, F, Cl, Br, I, @, +, -, =, #, (, ), [, ], %, digits
    import re
    basic_smiles_pattern = r'^[CNOSPFIHBrcbnospfih\[\]()@+=\-#0-9%\.]+$'
    if not re.match(basic_smiles_pattern, smiles):
        return False, "SMILES contains invalid characters. Use standard organic chemistry notation."
    
    return True, ""


def validate_smiles_rdkit(smiles: str) -> Tuple[bool, str]:
    """
    Advanced SMILES validation using RDKit.
    Only runs if RDKit is available.
    
    Returns (is_valid, error_message)
    """
    try:
        from rdkit import Chem
    except ImportError:
        # RDKit not available, fall back to basic validation
        return validate_smiles_basic(smiles)
    
    # First do basic validation
    is_valid, error_msg = validate_smiles_basic(smiles)
    if not is_valid:
        return False, error_msg
    
    # Try to parse with RDKit
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES: RDKit cannot parse this molecule"
        
        # Check for reasonable molecule size (prevent DoS)
        num_atoms = mol.GetNumAtoms()
        if num_atoms > 1000:
            return False, f"Molecule too large ({num_atoms} atoms, max 1000)"
        
        return True, ""
    except Exception as e:
        return False, f"SMILES parsing error: {str(e)}"


def validate_pdb_content(content: bytes, max_size_mb: int = 10) -> Tuple[bool, str]:
    """
    Validate PDB file content.
    Checks for:
    - File size limits
    - Required ATOM records
    - Basic PDB format
    
    Returns (is_valid, error_message)
    """
    # Check size
    size_mb = len(content) / (1024 * 1024)
    if size_mb > max_size_mb:
        return False, f"PDB file too large ({size_mb:.1f}MB, max {max_size_mb}MB)"
    
    # Convert to text
    try:
        text = content.decode('utf-8', errors='ignore')
    except Exception:
        return False, "Cannot decode PDB file as text"
    
    # Check for ATOM records
    lines = text.split('\n')
    has_atom = any(line.startswith('ATOM') or line.startswith('HETATM') for line in lines)
    if not has_atom:
        return False, "PDB file must contain ATOM or HETATM records"
    
    # Basic security check
    dangerous_patterns = ["<script", "DROP TABLE", "DELETE FROM"]
    text_upper = text.upper()
    for pattern in dangerous_patterns:
        if pattern.upper() in text_upper:
            return False, "File contains suspicious content"
    
    return True, ""


def validate_threshold(value: float) -> None:
    """
    Validate similarity threshold.
    Raises ValidationError if invalid.
    """
    if not 0.0 <= value <= 1.0:
        raise ValidationError(
            f"Threshold must be between 0.0 and 1.0 (got {value})"
        )


def sanitize_input(value: str) -> str:
    """
    Sanitize user input to prevent injection attacks.
    Removes dangerous characters while preserving valid chemical notation.
    """
    if not value:
        return ""
    
    # Strip whitespace
    value = value.strip()
    
    # Remove null bytes
    value = value.replace('\x00', '')
    
    # Remove control characters except newline/tab
    value = ''.join(char for char in value if ord(char) >= 32 or char in '\n\t')
    
    return value


# Example SMILES for documentation/testing
EXAMPLE_SMILES = {
    "ethanol": "CCO",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "bortezomib": "CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)c2cnccn2)B(O)O",
    "lenalidomide": "O=C1CCC(=O)N1c2cncc3n2c(nc3N4CCOCC4)N",
}
