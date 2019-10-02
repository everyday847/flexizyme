# flexizyme
Example workflow to set up flexizyme charging simulations

## What's in this repository

`flexizyme.py` will append a nonstandard residue to 3cun -- it's currently set up to do so for the existing Rosetta residue BZO -- and enumerate torsion angles to determine its most stable local conformation. The resulting structures can then be fed into broader sampling strategies. (We've used `relax`; other methods may be appropriate as well.)

`3cun.pdb` is a copy of the PDB structure 3CUN and is useful input to `flexizyme.py`

`three_p_acid_optimized.pdb` and `three_p_acid_unoptimized.pdb` are the results of running `flexizyme.py`; feel free to check your pyrosetta installation by comparing to these standard outputs. (There will likely be small numeric differences as Rosetta's code improves; there probably won't be large, dramatic changes.)

