import math, re

###################################################################################################################################
#Prompts to get the information required to begin the analysis
#Primer and buffer information is input here
#The primer sequences will be analyzed first by calculating their GC content
#Primer melting temperature calculation and adjustment occurs last
###################################################################################################################################
primer1 = input('Input primer 1 sequence here: ')
primer2 = input('Input primer 2 sequence here: ')
oligo1_conc = float(input('Input total single strand oligo concentration in uM here: '))
Mono = float(input('Input total monovalent ion concentration in mM here: '))
Mg = float(input('Input Mg concentration in mM here: '))
dNTPs = float(input('Input dNTP concentration in mM here: '))

gas_constant = 8.31445984848484848484 #J K-1-mol-1 
oligo_c = (oligo1_conc) * 1e-6 #Adjust concentration to mol/L
purine = ['A','G'] #Lists used to convert primer sequence into R/Y purine/pyrimidine nucleotide codes
pyrimidine = ['C','T']
primer1_length, primer2_length = len(primer1), len(primer2)
p1_terminal_base, p1_initial_base = primer1[-1], primer1[0]
p2_terminal_base, p2_initial_base = primer2[-1], primer2[0]

#Calculation of each primer's GC content
def get_primer_gc(primer1, primer2):
    primer1_GC = (sum((1.0 for base in primer1 if base in ['G', 'C'])) / primer1_length) * 100
    primer2_GC = (sum((1.0 for base in primer2 if base in ['G', 'C'])) / primer2_length) * 100
    return primer1_GC, primer2_GC

primer1_gc, primer2_gc = get_primer_gc(primer1, primer2)
###################################################################################################################################
#Determine terminal base compensation parameters for each primer
###################################################################################################################################
'''
For relatively short DNA sequences, such as PCR primers, the initial melting of the helix is influenced by the bases at the ends.
DNA sequences around this length follow appoximately two-state melting. Therefore, the initial 'unzipping' plays a larger factor
in helix melting. Regardless of the underlying mechanism, G or C bases should require more energy to begin unzipping.
'''
def p1_end_comp(p1_initial_base, p1_terminal_base):
    AT_h, AT_s, GC_h, GC_s = 2.53, 5.7, .49, -.7
    dh_comp, ds_comp = 0, 0
    if p1_initial_base in ['A', 'T']:
        dh_comp += AT_h
        ds_comp += AT_s
    if p1_initial_base in ['C', 'G']:
        dh_comp += GC_h
        ds_comp += GC_s
    if p1_terminal_base in ['A', 'T']:
        dh_comp += AT_h
        ds_comp += AT_s
    if p1_terminal_base in ['C', 'G']:
        dh_comp += GC_h
        ds_comp += GC_s
    return dh_comp, ds_comp

def p2_end_comp(p2_initial_base,p2_terminal_base):
    AT_h, AT_s, GC_h, GC_s = 2.53, 5.7, .49, -.7
    dh_comp, ds_comp = 0, 0
    if p2_initial_base in ['A', 'T']:
        dh_comp += AT_h
        ds_comp += AT_s
    if p2_initial_base in ['C', 'G']:
        dh_comp += GC_h
        ds_comp += GC_s
    if p2_terminal_base in ['A', 'T']:
        dh_comp += AT_h
        ds_comp += AT_s
    if p2_terminal_base in ['C', 'G']:
        dh_comp += GC_h
        ds_comp += GC_s
    return dh_comp, ds_comp

p1_term_h, p1_term_s = p1_end_comp(p1_initial_base,p1_terminal_base)
p2_term_h, p2_term_s = p2_end_comp(p2_initial_base,p2_terminal_base)
###################################################################################################################################
#Melting Temperature Calculation
###################################################################################################################################
'''
This section attempts to correct the model's accuracy discrepancy when computing the melting temperatures of oligos with 
different base sequences but the same GC content. Those types of oligos can differ significantly, however, Privalov and 
Crane-Robinson's model doesn't contain any way to account for those differences. I think the discrepancy may be due to intrinsic 
DNA bending, ion binding, and solvent effects. Site specific monovalent and divalent ion binding may break up water binding along 
the minor groove, thus decreasing the entropy cost associated with that water binding. In theory, the longer the alternating purine 
and pyrimidine chain, the more stable the helix is. Conversely, the longer the strings of repeating pyrimidine bases, the less stable 
the helix.
'''
GC_pairs_p1 = re.findall('(?:G|C){1}(?:G|C){1,}', primer1) #Generates a list of repeating guanine and cytosine sequences
GC_pairs_p2 = re.findall('(?:G|C){1}(?:G|C){1,}', primer2) #in each primer. Only includes sequences with two or more repeats.

RY_primer1_str = ''.join(['R' if base in purine else 'Y' for base in primer1]) #Converts each primer into the nucleotide codes
RY_primer2_str = ''.join(['R' if base in purine else 'Y' for base in primer2]) #R and Y for purine and pyrimidine bases

RY_pairs_p1 = re.findall('(?:RY){1,}', RY_primer1_str) #Generates a list of repeating purine and pyrimidine sequences
RY_pairs_p2 = re.findall('(?:RY){1,}', RY_primer2_str) #in each primer. Only includes sequences with two or more repeats.

YY_pairs_p1 = re.findall('(?:Y){2,}', RY_primer1_str) + re.findall('(?:R){2,}', RY_primer1_str) #Generates a list of repeating purine and
YY_pairs_p2 = re.findall('(?:Y){2,}', RY_primer2_str) + re.findall('(?:R){2,}', RY_primer2_str) #pyrimidine sequences in each primer.
'''
The above lists are needed so that their respective influences on helix stability can be measured and quantified. The functions 
below attempts to apply those differences to appropriately modify the melting temperature.
'''
#These variables take those captured sequences from the lists of pairs above and convert them into lengths. The lengths will then
#be used to alter the entropy values during the melting temperature calculation.
p1_GC_len = [len(seq) for seq in GC_pairs_p1]
p1_RY_len = [len(seq) for seq in RY_pairs_p1]
p1_YY_len = [len(seq) for seq in YY_pairs_p1]

p2_GC_len = [len(seq) for seq in GC_pairs_p2]
p2_RY_len = [len(seq) for seq in RY_pairs_p2]
p2_YY_len = [len(seq) for seq in YY_pairs_p2]
#Here the melting temperature of the primers is determined based on the method developed by Privalov and Crane-Robinson
#See Privalov, P. L., & Crane-Robinson, C. (2018). https://doi.org/10.1016/j.pbiomolbio.2018.01.007
heat_capacity = .13 #kJ/K,mol-bp
H_A_25, H_T_25 = 25, 25 #kJ/mol-bp
S_A_25, S_T_25 = 72, 72 #J/K-mol-bp
H_CG_25, S_CG_25 = 18.8, 44.7

def p1_melting_calculation(primer1):

    n_A, n_T = primer1.count('A'), primer1.count('T') #Counts the number of each base in the primer
    n_CG = primer1.count('C') + primer1.count('G') #Needed to calculate total enthalpy and entropy
    
    #Total enthalpy and entropy value calculations for determining the melting temperature. Takes the heat capacity into account.
    delta_H = ((H_A_25 + (heat_capacity * (314.5 - 298.15))) * n_A) + ((H_T_25 + (heat_capacity * (314.5 - 298.15))) * n_T) + ((H_CG_25 + (heat_capacity * (314.5 - 298.15))) * n_CG) + p1_term_h
    delta_S = ((S_A_25 + (heat_capacity * math.log(314.5 / 298.15))) * n_A) + ((S_T_25 + (heat_capacity * math.log(314.5 / 298.15))) * n_T) + ((S_CG_25 + (heat_capacity * math.log((314.5) / 298.15))) * n_CG) + p1_term_s
    '''
The equation below functions as an approximate "best fit" correction for the discrepancy in the model's accuracy at high and low %GC. It 
adjust the constants of a polynomial equation based on chnages in oligo length. The values here were derived from the melting temperatures
from the 92 oligos in Owczarzy's 2004 paper. See Owczarzy, R. et. al (2004). Biochemistry. https://doi.org/10.1021/bi034621r.
'''   
    delta_S -= ((-.001292 * (primer1_length) + .00988) * (primer1_gc) ** 2 + (.0572 * (primer1_length) - .7342) * (primer1_gc) + \
        (-.0366 * (primer1_length) ** 2 + .839 * (primer1_length) + 23.863))
    '''
Additional corrections for the combination of site specific DNA bending, ion binding to specific sequences, plus other assorted solvent
effects which are too complex to fully account for. The following try to account for divalent ion binding to G bases which may displace
part of the water spine therefore reducing entropy. The other YY and RY corrections try to account for the hightened twist-roll coupling
that occurs most often in pyrimidine-purine dimers and least often in purine-pyrimidine steps.
'''    
    if primer1_length < 25: #Corrections are most helpful for oligos shorter than around 25 bp 
        for num in p1_GC_len:
            delta_S -= (.70 * num)
        for num in p1_RY_len:
            delta_S -= (.25 * num)
        for num in p1_YY_len:
            delta_S += (.65 * num)
    tm = (1000 * delta_H) / (delta_S + (gas_constant * (math.log(2 / oligo_c)))) - 298.15

    return tm

def p2_melting_calculation(primer2):

    n_A, n_T, n_CG = primer2.count('A'), primer2.count('T'), primer2.count('C') + primer2.count('G')

    delta_H = ((H_A_25 + (heat_capacity * (314.5 - 298.15))) * n_A) + ((H_T_25 + (heat_capacity * (314.5 - 298.15))) * n_T) + ((H_CG_25 + (heat_capacity * (314.5 - 298.15))) * n_CG) + p2_term_h
    delta_S = ((S_A_25 + (heat_capacity * math.log(314.5 / 298.15))) * n_A) + ((S_T_25 + (heat_capacity * math.log(314.5 / 298.15))) * n_T) + ((S_CG_25 + (heat_capacity * math.log((314.5) / 298.15))) * n_CG) + p2_term_s
    
    delta_S -= ((-.001292 * (primer2_length) + .00988) * (primer2_gc) ** 2 + (.0572 * (primer2_length) - .7342) * (primer2_gc) + \
        (-.0366 * (primer2_length) ** 2 + .839 * (primer2_length) + 23.863))
    if primer2_length < 25:
        for num in p2_GC_len:
            delta_S -= (.70 * num)
        for num in p2_RY_len:
            delta_S -= (.25 * num)
        for num in p2_YY_len:
            delta_S += (.65 * num)
    tm = (1000 * delta_H) / (delta_S + (gas_constant * (math.log(2 / oligo_c)))) - 298.15

    return tm

primer1_melting_temperature = p1_melting_calculation(primer1)
primer2_melting_temperature = p2_melting_calculation(primer2)
###################################################################################################################################
#Salt and Magnesium Adjustments
###################################################################################################################################
#Adjustments and unit conversions for the chosen buffer conditions
Mon = Mono / 2.0 #Divide by two to account for the counterion present, e.g. Cl-, SO4-, etc.
mg_adj = Mg * 1e-3 #Converts to mol/L
mon = Mon * 1e-3
dntps = dNTPs * 1e-3 
ka = 3e4 #Association constant for the Mg2+--dNTP complex. Used to calculate the free magnesium in the buffer
mg = (-(ka * dntps - ka * mg_adj + 1.0) + math.sqrt((ka * dntps - ka * mg_adj + 1.0) ** 2 + 4.0 * ka * mg_adj)) / (2.0 * ka)

#Equations from Owczarzy 2008 that adjusts melting temperatures according to monovalent ion, magnesium, and dNTP concentrations
#See Owczarzy, R., et al. (2008). Biochemistry, https://doi.org/10.1021/bi702363u
def primer1_salt_correction(primer1_melting_temperature, primer1_gc, primer1_length, mg, mon):
    a, b, c, d, e, f, g = 3.92e-5, -9.11e-6, 6.26e-5, 1.42e-5, -4.02e-4, 5.10e-4, 8.31e-5 #Slightly modified constants
    '''
The condition calculating the melting temperature if there is no monovalent salt present must be put first, before calculating R.
This is because the R equation requires dividing by the monovalent ion concentration. And of course, you can't divide by zero.
Additionally, although the condition of R > 6 uses the same equation as when no monovalent ions are present, the two conditions
can't be combined. In situations where no magnesium is present, the equation would try to take the natural log of zero, which
throws an error.
'''
    if mon == 0:
        salt2 = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        return round((1 / salt2) - 273.15, 1)

    R = math.sqrt(mg) / mon #Ratio to determine whether monovalent or divalent ions are dominant in their effects on melting temperature

    if R < 0.22:
        salt1 = (1 / (primer1_melting_temperature + 273.15)) + ((4.29e-5 * (primer1_gc / 100)) - 3.95e-5) * math.log(mon) + 9.40e-6 * (math.log(mon)) ** 2
        return round((1 / salt1) - 273.15, 1)
    elif R < 6.0:
        a = 4.42e-5 * (0.843 - (0.352 * math.sqrt(mon) * math.log(mon)))
        d = 1.02e-5 * ((1.279 - 4.03e-3 * math.log(mon)) - 8.03e-3 * (math.log(mon) ** 2))
        g = 8.71e-5 * ((0.486 - 0.258 * math.log(mon)) + 5.25e-3 * (math.log(mon) ** 3)) 

        salt2 = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        return round((1 / salt2) - 273.15, 1)
    elif R > 6.0:
        salt2 = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
    return round((1 / salt2) - 273.15, 1)

def primer2_salt_correction(primer2_melting_temperature, primer2_gc, primer2_length, mg, mon):
    a, b, c, d, e, f, g = 3.92e-5, -9.11e-6, 6.26e-5, 1.42e-5, -4.02e-4, 5.10e-4, 8.31e-5
    if mon == 0:
        salt2 = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        return round((1 / salt2) - 273.15, 1)

    R = math.sqrt(mg) / mon 
    
    if R < 0.22:
        salt1 = (1 / (primer2_melting_temperature + 273.15)) + ((4.29e-5 * (primer2_gc / 100)) - 3.95e-5) * math.log(mon) + 9.40e-6 * (math.log(mon)) ** 2
        return round((1 / salt1) - 273.15, 1)
    elif R < 6.0:
        a = 4.42e-5 * (0.843 - (0.352 * math.sqrt(mon) * math.log(mon)))
        d = 1.02e-5 * ((1.279 - 4.03e-3 * math.log(mon)) - 8.03e-3 * (math.log(mon) ** 2))
        g = 8.71e-5 * ((0.486 - 0.258 * math.log(mon)) + 5.25e-3 * (math.log(mon) ** 3)) 

        salt2 = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        return round((1 / salt2) - 273.15, 1)
    elif R > 6.0:
        salt2 = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        return round((1 / salt2) - 273.15, 1)

adj_primer1_melting_temperature = primer1_salt_correction(primer1_melting_temperature, primer1_gc, primer1_length, mg, mon)
adj_primer2_melting_temperature = primer2_salt_correction(primer2_melting_temperature, primer2_gc, primer2_length, mg, mon)
###################################################################################################################################
#Printing the results
###################################################################################################################################
print('The GC content of primer 1 is ' + str(round(primer1_gc, 1)) + '%.' + ' The GC content of primer 2 is ' + str(round(primer2_gc, 1)) + '%.')
print('The melting temperatures of your primers are ' + str(adj_primer1_melting_temperature) + 'C for primer 1 and ' + str(adj_primer2_melting_temperature) + 'C for primer 2.')
