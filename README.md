# Blib2OpenSwath
A Python script to convert the spectral library in BLIB format (.blib) from either Skyline or BiblioSpec to readable library format (OpenSwath or Spectronaut tsv format)

## Installation of requirements

```
>pip install -r requirements.txt
```

## Usage

```
>python blib2openswath.py -h

usage: blib2openswath.py [-h] [--infile -i [-i ...]] [--fasta -f [-f ...]] [--tol -t [-t ...]] [--mz_type -m [-m ...]]

Blib2OpenSwath: Convert spectral library from BLIB (Skyline) to OpenSwath format

optional arguments:
  -h, --help            show this help message and exit
  --infile -i [-i ...]  The input spectral library file in BLIB format either from Skyline or BiblioSpec output
  --fasta -f [-f ...]   Proteome database in fasta format for mapping peptides sequences
  --tol -t [-t ...]     Library match tolerance in dalton (Da) for fragment m/z annotation (INFO: The tolerance of 0.5 Da and 0.05 Da was set as default)
  --mz_type -m [-m ...] Specify the type of fragment m/z values present in the input spectral library (Ex: "average" or "mono")
  --lib_fmt -l [-l ...] Specify the spectra library format (Ex: "openswath" or "spectronaut")
  --fmt -o [-o ...]     Specify the output spectra library file format (Ex: "tsv" or "csv")
```
## Example

```
>python blib2openswath.py --infile data\library\Spectral_library.blib --fasta data\database\sequence.fasta --tol 0.02 --mz_type mono 
```
## Contact
For more information, post an issue or send an email to chinnu.kemmaai@gmail.com
