from pyteomics import mass
import numpy as np

db = mass.Unimod()
aa_comp = dict(mass.std_aa_comp)

aa_comp['-H2O'] = db.by_title('Dehydrated')['composition']
aa_comp['p'] = db.by_title('Phospho')['composition']
aa_comp['-NH3'] = db.by_title('Ammonia-loss')['composition']
aa_comp['Deamidated'] = db.by_title("Deamidated")['composition']
aa_comp['-CO'] = db.by_title("Asp->Ser")['composition']

##print (mass.calculate_mass(composition=aa_comp["-NH3"]))
##print (mass.calculate_mass(composition=aa_comp["-CO"]))


def add_losses(ion_name, z, frag, frag_mz, losses, exclude, mass_type):
    
    ion_type = []
    theo_frag = []
    
    for loss in losses:
        if loss != exclude:
            loss_mz = mass.calculate_mass(composition=aa_comp[loss], average = bool(mass_type))
            ion_type.append(ion_name + loss)
            theo_frag.append(frag_mz - loss_mz/(z + 1))
        
    return np.asarray(ion_type), np.asarray(theo_frag)
    
def calc_theo_mz(pep, charge, ions, losses, mods, mass_type):

    if losses != None:
        losses == losses #['-H2O', '-NH3', '-CO']
        
    if ions == None:
        ions = ["a","b","y"]
    else:
        ions = list(ions)
        
    theo_frag_mz = np.array([], dtype=float)
    theo_ion_type = np.array([], dtype='object')
    #print (pep, charge, len(pep), mods, bool(mass_type))
    if mods != None and len(mods) != 1: ### For peptides with multiple modifications
        positions = []
        mod_mzs = []
        for m in mods:
            positions.append(int(m[0]))
            mod_mzs.append(float(m[1]))

        for z in range(charge-1):
            b_ion_series = [] ### Generate b-ions
            for idx, aa in enumerate(pep):
                b_ion_series.append(aa)
                b_frag = ''.join(b_ion_series)
                if b_frag != pep:
                    frag_mz = mass.calculate_mass(b_frag, ion_type='b', charge=z+1, average = bool(mass_type))
                    if idx + 1 < min(positions):
                        theo_frag_mz = np.append(theo_frag_mz, frag_mz)
                        theo_ion_type = np.append(theo_ion_type, "b" + str(idx + 1) + "[+" + str(z+1) + "]")
                        if losses != None:
                            frag_ion, frag_mz_loss = add_losses("b" + str(idx + 1) + "[+" + str(z+1) + "]", z, b_frag, frag_mz, losses, None, mass_type)
                            theo_frag_mz = np.append(theo_frag_mz, frag_mz_loss)
                            theo_ion_type = np.append(theo_ion_type, frag_ion)
                            
                    elif idx + 1 >= min(positions):
                        mod_mz = []
                        for mod_idx, p in enumerate(positions):
                            if idx + 1 >= p:
                                mod_mz.append(mod_mzs[mod_idx])

                        #print (b_frag, idx + 1, mod_idx, p, sum(mod_mz))
                        theo_frag_mz = np.append(theo_frag_mz, frag_mz + sum(mod_mz))
                        theo_ion_type = np.append(theo_ion_type, "b" + str(idx + 1) + "[+" + str(z+1) + "]+" + str(sum(mod_mz)))
                        if losses != None:
                            frag_ion, frag_mz_loss = add_losses("b" + str(idx + 1) + "[+" + str(z+1) + "]+" + str(sum(mod_mz)), z, b_frag, frag_mz, losses, None, mass_type)
                            theo_frag_mz = np.append(theo_frag_mz, frag_mz_loss)
                            theo_ion_type = np.append(theo_ion_type, frag_ion)

            y_ion_series = [] ### Generate y-ions
            for idx, aa in enumerate(pep[::-1]):
                y_ion_series.append(aa)
                y_frag = ''.join(y_ion_series)
                if y_frag != pep[::-1]:
                    frag_mz = mass.calculate_mass(y_frag, ion_type='y', charge=z+1, average = bool(mass_type))
                    if idx + 1 > len(pep)-max(positions):
                        mod_mz = []
                        for mod_idx, p in enumerate(positions):
                            if idx + 1 > len(pep[::-1])-p:
                                mod_mz.append(mod_mzs[mod_idx])
                                
                        #print (y_frag, idx + 1, len(pep), positions, len(pep)- max(positions), mod_mz)
                        theo_frag_mz = np.append(theo_frag_mz, frag_mz + sum(mod_mz))
                        theo_ion_type = np.append(theo_ion_type, "y" + str(idx + 1) + "[+" + str(z+1) + "]_" + str(sum(mod_mz)))
                        if losses != None:
                            frag_ion, frag_mz_loss = add_losses("y" + str(idx + 1) + "[+" + str(z+1) + "]_" + str(sum(mod_mz)), z, y_frag, frag_mz, losses, '-CO', mass_type)
                            theo_frag_mz = np.append(theo_frag_mz, frag_mz_loss)
                            theo_ion_type = np.append(theo_ion_type, frag_ion)
                        #print ("y" + str(idx + 1) + "[+" + str(z+1) + "]_" + str(mod_mz), y_frag, pep_mz + mod_mz)
                            
                    if idx + 1 <= len(pep)- max(positions):
                        #print (y_frag, idx + 1, len(pep), positions, len(pep)- max(positions), pep[::-1][idx + 1])
                        theo_frag_mz = np.append(theo_frag_mz, frag_mz)
                        theo_ion_type = np.append(theo_ion_type, "y" + str(idx + 1) + "[+" + str(z+1) + "]")
                        if losses != None:
                            frag_ion, frag_mz_loss = add_losses("y" + str(idx + 1) + "[+" + str(z+1) + "]", z, y_frag, frag_mz, losses, '-CO', mass_type)
                            theo_frag_mz = np.append(theo_frag_mz, frag_mz_loss)
                            theo_ion_type = np.append(theo_ion_type, frag_ion)
                            
    elif mods != None and len(mods) == 1: ### For peptides with single modification
        position = int(mods[0][0])
        mod_mz = float(mods[0][1])
        for z in range(charge-1):
            b_ion_series = []  ### Generate b-ions
            for idx, aa in enumerate(pep):
                b_ion_series.append(aa)
                b_frag = ''.join(b_ion_series)
                if b_frag != pep:
                    frag_mz = mass.calculate_mass(b_frag, ion_type='b', charge=z+1, average = bool(mass_type))
                    if idx + 1 < position:
                        theo_frag_mz = np.append(theo_frag_mz, frag_mz)
                        theo_ion_type = np.append(theo_ion_type, "b" + str(idx + 1) + "[+" + str(z+1) + "]")
                        if losses != None:
                            frag_ion, frag_mz_loss = add_losses("b" + str(idx + 1) + "[+" + str(z+1) + "]", z, b_frag, frag_mz, losses, None, mass_type)
                            theo_frag_mz = np.append(theo_frag_mz, frag_mz_loss)
                            theo_ion_type = np.append(theo_ion_type, frag_ion)
                    if idx + 1 >= position:
                        theo_frag_mz = np.append(theo_frag_mz, frag_mz + mod_mz)
                        theo_ion_type = np.append(theo_ion_type, "b" + str(idx + 1) + "[+" + str(z+1) + "]+" + str(mod_mz))
                        if losses != None:
                            frag_ion, frag_mz_loss = add_losses("b" + str(idx + 1) + "[+" + str(z+1) + "]+" + str(mod_mz), z, b_frag, frag_mz, losses, None, mass_type)
                            theo_frag_mz = np.append(theo_frag_mz, frag_mz_loss)
                            theo_ion_type = np.append(theo_ion_type, frag_ion)
                        #print ("b" + str(idx + 1) + "[+" + str(z+1) + "]_" + str(mod_mz), b_frag, pep_mz + mod_mz)

            y_ion_series = [] ### Generate y-ions
            for idx, aa in enumerate(pep[::-1]):
                y_ion_series.append(aa)
                y_frag = ''.join(y_ion_series)
                if y_frag != pep[::-1]:
                    frag_mz = mass.calculate_mass(y_frag, ion_type='y', charge=z+1, average = bool(mass_type))
                    if idx + 1 > len(pep)-position:
                        theo_frag_mz = np.append(theo_frag_mz, frag_mz + mod_mz)
                        theo_ion_type = np.append(theo_ion_type, "y" + str(idx + 1) + "[+" + str(z+1) + "]_" + str(mod_mz))
                        if losses != None:
                            frag_ion, frag_mz_loss = add_losses("y" + str(idx + 1) + "[+" + str(z+1) + "]_" + str(mod_mz), z, y_frag, frag_mz, losses, '-CO', mass_type)
                            theo_frag_mz = np.append(theo_frag_mz, frag_mz_loss)
                            theo_ion_type = np.append(theo_ion_type, frag_ion)
                        #print ("y" + str(idx + 1) + "[+" + str(z+1) + "]_" + str(mod_mz), y_frag, pep_mz + mod_mz)
                            
                    if idx + 1 <= len(pep)- position:
                        theo_frag_mz = np.append(theo_frag_mz, frag_mz)
                        theo_ion_type = np.append(theo_ion_type, "y" + str(idx + 1) + "[+" + str(z+1) + "]")
                        if losses != None:
                            frag_ion, frag_mz_loss = add_losses("y" + str(idx + 1) + "[+" + str(z+1) + "]", z, y_frag, frag_mz, losses, '-CO', mass_type)
                            theo_frag_mz = np.append(theo_frag_mz, frag_mz_loss)
                            theo_ion_type = np.append(theo_ion_type, frag_ion)
                        #print ("y" + str(idx + 1) + "[+" + str(z+1) + "]", y_frag, pep_mz)
                        
    elif mods == None: ### For peptides with no modifications
        for z in range(charge-1):
            b_ion_series = []  ### Generate b-ions
            for idx, aa in enumerate(pep):
                b_ion_series.append(aa)
                b_frag = ''.join(b_ion_series)
                if b_frag != pep:
                    frag_mz = mass.calculate_mass(b_frag, ion_type='b', charge=z+1, average = bool(mass_type))
                    theo_frag_mz = np.append(theo_frag_mz, frag_mz)
                    theo_ion_type = np.append(theo_ion_type, "b" + str(idx + 1) + "[+" + str(z+1) + "]")
                    if losses != None:
                        frag_ion, frag_mz_loss = add_losses("b" + str(idx + 1) + "[+" + str(z+1) + "]", z, b_frag, frag_mz, losses, None, mass_type)
                        theo_frag_mz = np.append(theo_frag_mz, frag_mz_loss)
                        theo_ion_type = np.append(theo_ion_type, frag_ion)
                    #print ("b" + str(idx + 1) + "[+" + str(z+1) + "]", b_frag, pep_mz)

            y_ion_series = [] ### Generate y-ions
            for idx, aa in enumerate(pep[::-1]):
                y_ion_series.append(aa)
                y_frag = ''.join(y_ion_series)
                if y_frag != pep[::-1]:
                    frag_mz = mass.calculate_mass(y_frag, ion_type='y', charge=z+1, average = bool(mass_type))
                    theo_frag_mz = np.append(theo_frag_mz, frag_mz)
                    theo_ion_type = np.append(theo_ion_type, "y" + str(idx + 1) + "[+" + str(z+1) + "]")
                    if losses != None:
                        frag_ion, frag_mz_loss = add_losses("y" + str(idx + 1) + "[+" + str(z+1) + "]", z, y_frag, frag_mz, losses, '-CO', mass_type)
                        theo_frag_mz = np.append(theo_frag_mz, frag_mz_loss)
                        theo_ion_type = np.append(theo_ion_type, frag_ion)
                    #print ("y" + str(idx + 1) + "[+" + str(z+1) + "]", y_frag, pep_mz)
                    
    return theo_ion_type, theo_frag_mz
