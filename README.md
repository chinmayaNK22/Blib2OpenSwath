# Blib2OpenSwath
A Python script to convert the spectral library in BLIB format (.blib) from either Skyline or BiblioSpec to readable library format (OpenSwath or Spectronaut tsv format)

## Installation of requirements

```
>pip install -r requirements.txt
```

## Usage of Blib2OpenSwath

```
>python3 blib2openswath.py -h

usage: blib2openswath.py [-h] [--infile -i [-i ...]] [--fasta -f [-f ...]] [--tol -t [-t ...]] [--mz_type -m [-m ...]]

Blib2OpenSwath: Convert spectral library from BLIB (Skyline) to OpenSwath format

optional arguments:
  -h, --help            show this help message and exit
  --infile -i [-i ...]  The input spectral library file in BLIB format either from Skyline or BiblioSpec output
  --fasta -f [-f ...]   Proteome database in fasta format for mapping peptides sequences
  --tol -t [-t ...]     Library match tolerance in dalton (Da) for fragment m/z annotation (INFO: The tolerance of 0.5 Da and 0.05 Da was set as default)
  --mz_type -m [-m ...] Specify the type of fragment m/z values present in the input spectral library (Ex: "average" or "mono")
```
## Example

```
>python3 blib2openswath.py --infile 
```
## Contact
For more information, post an issue or send an email to chinnu.kemmaai@gmail.com
