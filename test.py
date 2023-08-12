from get_SA import get_modified_SA, get_SA
import numpy as np
amino_mass_dic = {'G' : 57.02, 'A' : 71.04, 'S' : 87.03, 'P' : 97.05, 'V' : 99.07, 'T' : 101.05, 'C' : 103.01, 'I' : 113.08, 'L' : 113.08, 'N' : 114.04, 'D' : 115.03, 'Q'  : 128.06, 'K' : 128.09, 'E' : 129.04, 'M' : 131.04, 'H' : 137.06, 'F' : 147.07, 'R' : 156.10, 'Y' : 163.06, 'W' : 186.08}

def make_spectrum(sequence):
    neutral_mass = 0
    for aminoacid in sequence:
        neutral_mass += amino_mass_dic[aminoacid]
    #charge = 2라고 가정
    b_ions = []
    y_ions = []
    cur_mass = 0
    for aminoacid in sequence:
        cur_mass += amino_mass_dic[aminoacid]
        b_ions.append(cur_mass/2)
        y_ions.append((neutral_mass - cur_mass)/2)
    b_ions.pop()
    y_ions.pop()
    ion_list = b_ions + y_ions
    # print(ion_list, b_ions, y_ions)
    spectrum = [(i,j*100) for j,i in enumerate(ion_list)]
    # print(spectrum)
    return spectrum


#17 modification이 있다고 가정
def make_mod_spectrum(sequence,mod_idx):
    neutral_mass = 17
    for aminoacid in sequence:
        neutral_mass += amino_mass_dic[aminoacid]
    #charge = 2라고 가정
    # print(neutral_mass)
    b_ions = []
    y_ions = []
    cur_mass = 0
    for idx,aminoacid in enumerate(sequence):
        cur_mass += amino_mass_dic[aminoacid]
        if idx >= mod_idx:
            b_ions.append(((cur_mass+17)/2))
            y_ions.append((neutral_mass - cur_mass-17)/2)
            # print(b_ions[-1], y_ions[-1], (b_ions[-1] + y_ions[-1])*2)
        else:
            b_ions.append((cur_mass)/2)
            y_ions.append((neutral_mass - cur_mass)/2)
            # print(b_ions[-1], y_ions[-1], (b_ions[-1] + y_ions[-1])*2)
    b_ions.pop()
    y_ions.pop()
    ion_list = b_ions + y_ions
    # print(ion_list, b_ions, y_ions)
    spectrum = [(i,j*100) for j,i in enumerate(ion_list)]
    # print(spectrum)
    return spectrum
spec1 = make_spectrum("GAAGSPVVTAGQDNNQ")
spec2 = make_mod_spectrum("GAAGSPVVTAGQDNNQ",3)
spec1.sort()
spec2.sort()
print(spec1)
print(spec2)

print(get_modified_SA(spec2, spec1,17))
print(get_SA(spec2, spec1))
