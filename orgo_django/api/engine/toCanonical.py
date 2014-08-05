"""
Convert any SMILES string into a canonical SMILES string.

Refer to:
http://openbabel.org/dev-api/canonical_code_algorithm.shtml

Public-facing methods:
    `to_canonical`
"""

import openbabel

def to_canonical(smiles):
    """
    Uses OpenBabel.
    Converts a SMILES string to a canonical SMILES string.
    smiles :: str.
    return :: str.
    """
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "can")
    outMol = openbabel.OBMol()
    obConversion.ReadString(outMol, str(smiles))
    ans = obConversion.WriteString(outMol)
    return ans.strip()
