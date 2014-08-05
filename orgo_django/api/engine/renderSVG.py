"""
renderSVG.py

Contains the function render(smiles), which takes in a string in SMILES format
and outputs a SVG representation (as a string) of that molecule. Uses the
OpenBabel Python library:
http://openbabel.org/wiki/Python
"""

import openbabel
import re

#Set up input and output formats

def render(smiles, hydrogens=False):
    """
    smiles :: str. In SMILES format.
    hydrogens :: bool. Currently, the hydrogens option does nothing.
    return :: str. In SVG format.
    """
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "svg")
    outMol = openbabel.OBMol()
    obConversion.ReadString(outMol, str(smiles))
    ans = obConversion.WriteString(outMol)
    
    ## Make the svg background transparent:
    ## replace fill="rgb(255,255,255)" with fill-opacity="0"
    ans = re.sub("fill=\"white\"", "fill-opacity=\"0\"", ans)

    return ans
