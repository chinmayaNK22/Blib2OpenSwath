import sqlite3
from scripts.read_fasta_file import readfasta
from scripts.protein_digestor import ProteinDigestion
from pyteomics import mass, parser
import numpy as np
from struct import *
from scripts import pyteomics_utils
from scripts import parse_binary
import os
import argparse


parser = argparse.ArgumentParser(description='''Blib2OpenSwath: Convert spectral library from BLIB (Skyline) to OpenSwath format''')

parser.add_argument('--infile', metavar='-i', type=str, nargs='+', help='The input spectral library file in BLIB format either from Skyline or BiblioSpec output')

parser.add_argument('--fasta',metavar='-f', type=str, nargs='+', help='Proteome database in fasta format for mapping peptides sequences')

parser.add_argument('--tol',metavar='-t', type=float, nargs='+', help='Library match tolerance in dalton (Da) for fragment m/z annotation (INFO: The tolerance of 0.5 Da and 0.05 Da was set as default)')

parser.add_argument('--mz_type',metavar='-m', type=str, nargs='+', help='Specify the type of fragment m/z values present in the input spectral library (Ex: "average" or "mono")')

args = parser.parse_args()

def fetch_sqlite3_db(infile, table):

    try:
        con = sqlite3.connect(infile)
        cur = con.cursor()
        if table == 'RefSpectra':
            cur.execute("SELECT * FROM RefSpectra")
            names = list(map(lambda x: x[0], cur.description))
            return names, cur.fetchall()

        elif table == 'RefSpectraPeaks':
            cur.execute("SELECT * FROM RefSpectraPeaks")
            #names = list(map(lambda x: x[0], cur.description))
            return cur.fetchall()

        elif table == 'Modifications':
            cur.execute("SELECT * FROM Modifications")
            return cur.fetchall()

        elif table == 'SpectrumSourceFiles':
            cur.execute("SELECT * FROM SpectrumSourceFiles")
            return cur.fetchall()
    
        con.close()
    except:
        raise Exception(f"The input spectral library {infile} is not a BiblioSpec or Skyline generated spectral library")


def digest_pro(fasta):

    protease = "trypsin"
    missed_cleaves = 2
    min_len = 6
    max_len = 80
    
    print (f"INFO: The input fasta sequences were theoretically digested with {protease} protease enzyme")
    print (f"INFO: Peptides having {min_len} to {max_len} amino acids will be generated with maximum missed cleaves of {missed_cleaves}")
    
    dig_peps = {}
    for rows in readfasta(fasta).read():
        for x in range(0, missed_cleaves+1):
            if protease == 'trypsin':
                peps = ProteinDigestion(rows[1],x,min_len,max_len).trypsin()
                for pep in peps:
                    if pep not in dig_peps:
                        dig_peps[pep] = [rows[0].split(' ')[0]]
                    else:
                        dig_peps[pep].append(rows[0].split(' ')[0])

    print (f'INFO: The theoretical {protease} digestion of {os.path.split(fasta)[1]} has resulted in {len(dig_peps)} tryptic peptides')
    
    return dig_peps

def annotate_lib(infile, fasta):
    name, pepinfo = fetch_sqlite3_db(infile, 'RefSpectra')

    pep_idx = name.index('peptideSeq')
    #print (name, pep_idx)
    output = []
    theo_peps = digest_pro(fasta)
    print (f'The theoretical digestion of {fasta} has resulted in {len(theo_peps)}')
    
    c = 0
    a = 0
    for pinfo in pepinfo:
        peptide = pinfo[pep_idx]
        if peptide in theo_peps:
            c += 1
            output.append(list(map(str,list(pinfo))) + [';'.join(theo_peps[peptide])] + [str(len(list(theo_peps[peptide])))])

        else:
            a += 1
            output.append(list(map(str,list(pinfo))) + ["Not Found","NA"])

    print (f'There were {c} PSMs annotated to the proteome and {a} were not found annotated')   
    outfile = "{0}_protein_annotated.txt".format(infile.rstrip('.blib'))
    with open(outfile, 'w') as outf:
        outf.write('\t'.join(name) + "\tAccession\tNo. of Proteins\n" )
        outf.writelines('\t'.join(i) + '\n' for i in output)


def annotate_pep(peptide, theo_peps):
    if peptide in theo_peps:
        return ';'.join(theo_peps[peptide])

    else:
        return "Not Found"

def fragments(peptide, types=('b', 'y'), maxcharge=1):
    """
    The function generates all possible m/z for fragments of types
    `types` and of charges from 1 to `maxharge`.
    """
    for i in range(1, len(peptide)-1):
        for ion_type in types:
            for charge in range(1, maxcharge+1):
                if ion_type[0] in 'abc':
                    yield mass.calculate_mass(peptide[:i], ion_type=ion_type, charge=charge)
                else:
                    print (peptide[i:])
                    yield mass.calculate_mass(peptide[i:], ion_type=ion_type, charge=charge)


def frags(pep):
    tsg = TheoreticalSpectrumGenerator()
    peptide = AASequence.fromString(pep)
    # standard behavior is adding b- and y-ions of charge 1
    p = Param()
    spec2 = MSSpectrum()
    p.setValue("add_b_ions", "true")
    p.setValue("add_a_ions", "true")
    p.setValue("add_losses", "true")
    p.setValue("add_metainfo", "true")
    tsg.setParameters(p)
    tsg.getSpectrum(spec2, peptide, 1, 1)

    # Iterate over annotated ions and their masses
    #print("Spectrum 2 of", peptide, "has", spec2.size(), "peaks.")
    frag_ions = np.array([], dtype = 'object')
    frag_mz = np.array([], dtype=float)
    for ion, peak in zip(spec2.getStringDataArrays()[0], spec2):
        frag_ions = np.append(frag_ions, ion.decode())
        frag_mz = np.append(frag_mz, peak.getMZ())

    #fragments = np.stack((frag_ions, frag_mz), axis=0)

    return frag_ions, frag_mz

def align_spectra(theo_spectrum, observed_spectrum):
    alignment = []
    spa = SpectrumAlignment()
    p = spa.getParameters()
    # use 0.5 Da tolerance (Note: for high-resolution data we could also use ppm by setting the is_relative_tolerance value to true)
    p.setValue("tolerance", 0.5)
    p.setValue("is_relative_tolerance", "false")
    spa.setParameters(p)
    # align both spectra
    spa.getSpectrumAlignment(alignment, theo_spectrum, observed_spectrum)

    print("Number of matched peaks: " + str(len(alignment)))
    print("ion\ttheo. m/z\tobserved m/z")

def map_frags(infile, fasta, tolerance, mass_type):
    theo_peps = digest_pro(fasta)
    print (f"INFO: A library match tolerance of {tolerance} Da was set for the spectral library generation")

    header, pep = fetch_sqlite3_db(infile, 'RefSpectra')
    numpeak_idx = header.index('numPeaks')
    rt_idx = header.index('retentionTime')
    mz_idx = header.index('precursorMZ') 

    mods = fetch_sqlite3_db(infile, 'Modifications')
    Mods = {}
    all_mods = {}
    for mod in mods:
        all_mods[mod[-1]] = 1
        if mod[1] not in Mods:
            Mods[mod[1]] = [mod[2:]]
        else:
            Mods[mod[1]].append(mod[2:])

    print (f'INFO: There were {len(list(all_mods))} modifications found in the library')

    peaks = fetch_sqlite3_db(infile, 'RefSpectraPeaks')
    rawfiles = fetch_sqlite3_db(infile, 'SpectrumSourceFiles')
    print (f'INFO: The input library {infile} consists of database search results from {len(rawfiles)} raw files')

    sn_threshold = 0
    if mass_type == "average":
        tolerance = 0.5
        sn_threshold = 1
    else:
        tolerance = 0.05
        sn_threshold = 0

    print (f'INFO: Fragments with intensity values above {sn_threshold}% of Signal to Noise threshold will be considered')

    output = []
    frags_count = {}
    all_peps = {}
    for p in pep:
        numpeaks = str(p[numpeak_idx])
        pep = p[1]
        modpep = p[4]
        z = p[3]
        iD = p[0]

        #### Generate theoretical fragment ions with or without modifications accordingly
        if iD in Mods:
            if mass_type == 'average':
                ion_type, theo_mz = pyteomics_utils.calc_theo_mz(pep, z, None, ['-H2O', '-NH3', '-CO'], Mods[iD], True)
            else:
                ion_type, theo_mz = pyteomics_utils.calc_theo_mz(pep, z, None, ['-H2O', '-NH3', '-CO'], Mods[iD], False)

        else:
            if mass_type == 'average':
                ion_type, theo_mz = pyteomics_utils.calc_theo_mz(pep, z, None, ['-H2O', '-NH3', '-CO'], None, True)
            else:
                ion_type, theo_mz = pyteomics_utils.calc_theo_mz(pep, z, None, ['-H2O', '-NH3', '-CO'], None, False)
            
        #theo_ions, theo_mz = frags(modpep)

        spectra = peaks[iD-1] ## fetch peptide corresponding spectra
        ID = spectra[0]
        
        #### Unpack the mz peaks and intensity values stored in binary format
        mz = spectra[1]
        ex_frag_mz = parse_binary.unpack_array(mz, numpeaks, 'peaks')
        
        intensity = spectra[2]
        ex_frag_intensity = parse_binary.unpack_array(intensity, numpeaks, 'intensity')

        max_ints = ex_frag_intensity.max()
        
        matching_values = np.array([], dtype=float)
        matching_ions = np.array([], dtype='object')
        matching_intense = np.array([], dtype='f')
        
        # Loop through the values in thoeretical mz array and identify the matching experimental mz values within given fragment tolerance
        for index, value in np.ndenumerate(theo_mz):
            if np.any(np.abs(ex_frag_mz - value) <= float(tolerance)):
                exp_mz = ex_frag_mz[np.where(np.abs(ex_frag_mz - value) <= float(tolerance))]
                
                if len(exp_mz) > 1:
                    close_exp_idx = np.abs(exp_mz - value).argmin()
                    new_exp_mz = exp_mz[close_exp_idx]
                    if len(np.where(matching_values == new_exp_mz)[0]) == 0: ### Exclude the mz values which already present the the final numpy array
                        matching_values = np.append(matching_values, new_exp_mz)
                        matching_ions = np.append(matching_ions, ion_type[index])
                        if close_exp_idx == len(ex_frag_intensity)-1:
                            matching_intense = np.append(matching_intense, ex_frag_intensity[-1])
                        
                        else:
                            matching_intense = np.append(matching_intense, ex_frag_intensity[close_exp_idx])
                            
                else:
                    exp_mz_idx = np.where(np.abs(ex_frag_mz - value) <= float(tolerance))[0]
                    if len(np.where(matching_values == exp_mz)[0]) == 0: ### Exclude the mz values which already present the the final numpy array
                        matching_values = np.append(matching_values, exp_mz)
                        matching_ions = np.append(matching_ions, ion_type[index])
                        #print (exp_mz_idx, exp_mz, value, ion_type[index], len(ex_frag_intensity), len(ex_frag_mz), ex_frag_mz)
                        if exp_mz_idx == len(ex_frag_intensity)-1 or exp_mz_idx == len(ex_frag_intensity):
                            matching_intense = np.append(matching_intense, ex_frag_intensity[-1])
                        else:
                            try:
                                matching_intense = np.append(matching_intense, ex_frag_intensity[exp_mz_idx])
                            except:
                                matching_intense = np.append(matching_intense, 0)

        if modpep + '@' + str(z) not in all_peps:
            all_peps[modpep + '@' + str(z)] = [str(len(matching_values)) + '@' + str(iD)]
            try:
                if len(matching_ions) >= 3:
                    for index, ions in np.ndenumerate(matching_ions):
                        proteinID = annotate_pep(pep, theo_peps)
                        if (matching_intense[index]/max_ints)*100 > sn_threshold and proteinID != "Not Found": ### Consider fragments with S/N ratio > S/N threshold and matching to only the target fasta sequences
                            transition_group_id = modpep + '+' + str(z)
                            transition_name = transition_group_id + '_' + ions
                            pep = pep
                            modpep = modpep
                            fullpepname = pep
                            rt = str(round(p[rt_idx]*60,4))
                            prec_mz = str(round(p[mz_idx],5))
                            charge = str(z)
                            prod_mz = str(matching_values[index])
                            prod_z = str(ions.split('+')[1][0])
                            lib_intense = str(round(matching_intense[index],4))
                            frag_type = str(ions[0])
                            frag_num = ions.split('+')[0].rstrip('[').lstrip(frag_type)
                            decoy = "0"
                            quant_trans = "1"
                            #print (transition_group_id, transition_name, proteinID, pep, modpep, fullpepname, rt, prec_mz, charge, prod_mz, prod_z, lib_intense, frag_type, frag_num, decoy, quant_trans)
                            output.append([transition_group_id, transition_name, proteinID, pep, modpep, fullpepname, rt, prec_mz, charge, prod_mz, prod_z, lib_intense, frag_type, frag_num, decoy, quant_trans])

            except:
                print (f'Cound not write the peptide transitions of {modpep}')

            
        else:
            all_peps[modpep + '@' + str(z)].append(str(len(matching_values)) + '@' + str(iD))

        
    print (f'INFO: The library contains annotated peaks for {len(all_peps)} peptide precursors')
    print (f'INFO: Peptides with less than 3 transitions were excluded in the output library')
    print (f'INFO: {len(output)} transitions corresponding to {len(all_peps)} peptide precursors were converted to a OpenSwath library format')
    
    outfile = "{0}_OpenSwath.tsv".format(infile.rstrip('blib').rstrip('.'))
    with open(outfile, 'w') as outf:
        outf.write("transition_group_id\ttransition_name\tProteinId\tPeptideSequence\tModifiedPeptideSequence\tFullPeptideName\tRetentionTime\tPrecursorMz\tPrecursorCharge\tProductMz\tProductCharge\tLibraryIntensity\tFragmentIonType\tFragmentSeriesNumber\tIsDecoy\tquantifying_transition\n")
        outf.writelines('\t'.join(i) + '\n' for i in output)

    print (f"INFO: A OpenSwath compatible spectral library {outfile} was generated")

#annotate_lib(infile, fasta)
    
if __name__== "__main__":

    if os.path.isfile(os.path.join(args.infile[0])) and args.infile[0].split('.')[-1] == 'blib':
        map_frags(args.infile[0], args.fasta[0], args.tol[0], args.mz_type[0])
        print (f"FINISHED: Completed the convertion of {args.infile[0]} to OpenSwath format")

    #### Convert multiple BLIB spectral libraries present in a given input folder in loop
    elif os.path.isdir(os.path.join(args.infile[0])):
        infiles = []
        for blibs in os.listdir(os.path.join(args.infile[0])):
            if os.path.isfile(os.path.join(args.infile[0], blibs)) and blibs.split('.')[-1] == 'blib':
                infiles.append(os.path.join(args.infile[0], blibs))

        print (f"INFO: Found {len(infiles)} BLIB spectral library files in the given input path {args.infile[0]}")

        for idx, infile in enumerate(infiles):
            print (f"INFO: Processing file {idx + 1} of {len(infiles)} BLIB spectral libraries in {args.infile[0]}")
            map_frags(infile, args.fasta[0], args.tol[0], args.mz_type[0])
            print (f"FINISHED: Completed the convertion of {os.path.split(infile)[1]} to OpenSwath format")
