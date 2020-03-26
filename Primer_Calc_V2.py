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

#Calculation of each primer's GC content
def get_primer_gc(primer1, primer2):
    primer1_GC = (sum((1.0 for base in primer1 if base in ['G', 'C'])) / primer1_length) * 100
    primer2_GC = (sum((1.0 for base in primer2 if base in ['G', 'C'])) / primer2_length) * 100
    return primer1_GC, primer2_GC

primer1_gc, primer2_gc = get_primer_gc(primer1, primer2)
###################################################################################################################################
#Melting Temperature Calculation
###################################################################################################################################
'''
This section attempts to correct the model's accuracy discrepancy when computing the melting temperatures of oligos with 
different base sequences but the same GC content. Those types of oligos can differ significantly, however, 
Privalov and Crane-Robinson's model doesn't contain any way to account for those differences. I think the discrepancy 
may be due to the presence of alternating purine and pyrimidine bases. The alternating pairs may break up water binding
to adenine and thymine residues, thus decreasing the entropy cost associated with that water binding. In theory, the longer 
the alternating purine and pyrimidine chain, the more stable the helix is. Conversely, the more long strings of repeating
pyrimidine bases, the less stable the helix.
'''
AT_pairs_p1 = re.findall('(?:A|T){1}(?:A|T){1,}', primer1) #Generates a list of repeating adenine and thymine sequences
AT_pairs_p2 = re.findall('(?:A|T){1}(?:A|T){1,}', primer2) #in each primer. Only includes sequences with two or more repeats.

GC_pairs_p1 = re.findall('(?:G|C){1}(?:G|C){1,}', primer1) #Generates a list of repeating guanine and thymine cytosine
GC_pairs_p2 = re.findall('(?:G|C){1}(?:G|C){1,}', primer2) #in each primer. Only includes sequences with two or more repeats.

RY_primer1_str = ''.join(['R' if base in purine else 'Y' for base in primer1])#Converts each primer into the nucleotide codes
RY_primer2_str = ''.join(['R' if base in purine else 'Y' for base in primer2])#R and Y for purine and pyrimidine bases

RY_pairs_p1 = re.findall('(?:RY){2,}', RY_primer1_str) #Generates a list of repeating purine and pyrimidine sequences
RY_pairs_p2 = re.findall('(?:RY){2,}', RY_primer2_str) #in each primer. Only includes sequences with two or more repeats.

YY_pairs_p1 = re.findall('(?:YY){2,}', RY_primer1_str) #Generates a list of repeating pyrimidine sequences
YY_pairs_p2 = re.findall('(?:YY){2,}', RY_primer2_str) ##in each primer. Only includes sequences with two or more repeats.
'''
The above lists are needed so that their respective influences on helix stability can be measured and quantified, assuming
there is some influence. Don't yet know how best to do that. The functions below attempts to apply those differences to appropriately 
modify the melting temperature.
'''
#These variables take those captured sequences from the lists of pairs above and convert them into lengths. The lengths will then
#be used to alter the entropy values during the melting temperature calculation.
p1_AT_len = [len(seq) for seq in AT_pairs_p1]
p1_GC_len = [len(seq) for seq in GC_pairs_p1]
p1_RY_len = [len(seq) for seq in RY_pairs_p1]
p1_YY_len = [len(seq) for seq in YY_pairs_p1]

p2_AT_len = [len(seq) for seq in AT_pairs_p2]
p2_GC_len = [len(seq) for seq in GC_pairs_p2]
p2_RY_len = [len(seq) for seq in RY_pairs_p2]
p2_YY_len = [len(seq) for seq in YY_pairs_p2]
#Here the melting temperature of the primers is determeined based on the method developed by Privalov and Crane-Robinson
#See Privalov, P. L., & Crane-Robinson, C. (2018). https://doi.org/10.1016/j.pbiomolbio.2018.01.007
heat_capacity = .13 #kJ/K,mol-bp
H_A_25, H_T_25 = 25, 25 #kJ/mol-bp
S_A_25, S_T_25 = 72, 72 #J/K-mol-bp
H_CG_25, S_CG_25 = 18.8, 44.7
delta_S_trans = gas_constant * math.log(2 / oligo_c)

def p1_melting_calculation(primer1):

    n_A, n_T = primer1.count('A'), primer1.count('T') #Counts the number of each base in the primer
    n_CG = primer1.count('C') + primer1.count('G') #Needed to calculate total enthalpy and entropy
    
    #Total enthalpy and entropy value calculations for determining the melting temperature. Takes the heat capacity into account.
    delta_H = ((H_A_25 + (heat_capacity * (314.5 - 298.15))) * n_A) + ((H_T_25 + (heat_capacity * (314.5 - 298.15))) * n_T) + ((H_CG_25 + (heat_capacity * (314.5 - 298.15))) * n_CG)
    delta_S = ((S_A_25 + (heat_capacity * math.log(314.5 / 298.15))) * n_A) + ((S_T_25 + (heat_capacity * math.log(314.5 / 298.15))) * n_T) + ((S_CG_25 + (heat_capacity * math.log((314.5) / 298.15))) * n_CG)
    if primer1_length > 15:
        delta_S -= (-.0251 * (primer1_gc) ** 2 + .9798 * (primer1_gc) + 11.852)
    else:
        delta_S -= (-.0064 * (primer1_gc) ** 2 - .0647 * (primer1_gc) + 23.153)
    '''
Corrections for destabilzing synergystic effects of water chaining and the stabilizing effects of disrupting that chaining. While RY/YY
combos are imperfect, as an RY pair may be ATA, which could form water chains and vice versa for YY pairs, adding them both should roughly
cancel out over/underestimations of their stabilization effects.
'''    
    for num in p1_AT_len:
        if num >= 2:
            delta_S += (.45 * num)
    for num in p1_GC_len:
        if num > 2:
            delta_S -= (1.85 * num)
    for num in p1_RY_len:
        if num > 3:
            delta_S -= (.45 * num) * (math.lgamma(primer1_gc) / 100)
    for num in p1_YY_len:
        if num > 3:
            delta_S += (.57 * num)

    tm = (1000 * delta_H) / (delta_S + (gas_constant * (math.log(2 / oligo_c)))) - 298.15
    return tm

def p2_melting_calculation(primer2):

    n_A, n_T, n_CG = primer2.count('A'), primer2.count('T'), primer2.count('C') + primer2.count('G')

    delta_H = ((H_A_25 + (heat_capacity * (314.5 - 298.15))) * n_A) + ((H_T_25 + (heat_capacity * (314.5 - 298.15))) * n_T) + ((H_CG_25 + (heat_capacity * (314.5 - 298.15))) * n_CG)
    delta_S = ((S_A_25 + (heat_capacity * math.log(314.5 / 298.15))) * n_A) + ((S_T_25 + (heat_capacity * math.log(314.5 / 298.15))) * n_T) + ((S_CG_25 + (heat_capacity * math.log((314.5) / 298.15))) * n_CG)
    if primer2_length > 15:
        delta_S -= (-.0251 * (primer2_gc) ** 2 + .9798 * (primer2_gc) + 11.852)
    else:
        delta_S -= (-.0064 * (primer2_gc) ** 2 - .0647 * (primer2_gc) + 23.153)

    for num in p2_AT_len:
        if num >= 2:
            delta_S += (.45 * num)
    for num in p2_GC_len:
        if num > 2:
            delta_S -= (1.85 * num)
    for num in p2_RY_len:
        if num > 3:
            delta_S -= (.45 * num) * (math.lgamma(primer2_gc) / 100)
    for num in p2_YY_len:
        if num > 3:
            delta_S += (.57 * num)

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
    a, b, c, d, e, f, g = 3.92e-5, -9.11e-6, 6.66e-5, 1.42e-5, -4.02e-4, 4.05e-4, 8.31e-5 #Slightly modified constants
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
        a = 3.92e-5 * (0.843 - (0.352 * math.sqrt(mon) * math.log(mon)))
        d = 1.42e-5 * ((1.279 - 4.03e-3 * math.log(mon)) - 8.03e-3 * (math.log(mon) ** 2))
        g = 8.31e-5 * ((0.486 - 0.258 * math.log(mon)) + 5.25e-3 * (math.log(mon) ** 3)) 

        salt2 = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        return round((1 / salt2) - 273.15, 1)
    elif R > 6.0:
        salt2 = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
    return round((1 / salt2) - 273.15,1)

def primer2_salt_correction(primer2_melting_temperature, primer2_gc, primer2_length, mg, mon):
    a, b, c, d, e, f, g = 3.92e-5, -9.11e-6, 6.66e-5, 1.42e-5, -4.02e-4, 4.05e-4, 8.31e-5
    if mon == 0:
        salt2 = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        return round((1 / salt2) - 273.15, 1)

    R = math.sqrt(mg) / mon 
    
    if R < 0.22:
        salt1 = (1 / (primer2_melting_temperature + 273.15)) + ((4.29e-5 * (primer2_gc / 100)) - 3.95e-5) * math.log(mon) + 9.40e-6 * (math.log(mon)) ** 2
        return round((1 / salt1) - 273.15, 1)
    elif R < 6.0:
        a = 3.92e-5 * (0.843 - (0.352 * math.sqrt(mon) * math.log(mon)))
        d = 1.42e-5 * ((1.279 - 4.03e-3 * math.log(mon)) - 8.03e-3 * (math.log(mon) ** 2))
        g = 8.31e-5 * ((0.486 - 0.258 * math.log(mon)) + 5.25e-3 * (math.log(mon) ** 3))

        salt2 = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        return round((1 / salt2) - 273.15, 1)
    elif R > 6.0:
        salt2 = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        return round((1 / salt2) - 273.15,1)

adj_primer1_melting_temperature = primer1_salt_correction(primer1_melting_temperature, primer1_gc, primer1_length, mg, mon)
adj_primer2_melting_temperature = primer2_salt_correction(primer2_melting_temperature, primer2_gc, primer2_length, mg, mon)
###################################################################################################################################
#Printing the results
###################################################################################################################################
print('The GC content of primer 1 is ' + str(round(primer1_gc, 1)) + '%.' + 
' The GC content of primer 2 is ' + str(round(primer2_gc, 1)) + '%.')
print('The melting temperatures of your primers are ' + str(adj_primer1_melting_temperature) + 
'C for primer 1, and ' + str(adj_primer2_melting_temperature) + 'C for primer 2.')