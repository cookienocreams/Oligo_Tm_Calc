import math,re

###################################################################################################################################
#Prompts to get the information required to begin the analysis
#Primer and buffer information is input here
#The primer sequences will be analyzed first by calculating their GC content
#Primer melting temperature calculation and adjustment occurs last
###################################################################################################################################

primer1 = input('Input primer 1 sequence here: ') #Right now there are two inputs for non-complentary primers
primer2 = input('Input primer 2 sequence here: ') #Will either change to a single input or figure out how to expand to batch primers
oligo1_conc = float(input('Input single strand oligo concentration in uM here: '))
Mono = float(input('Input total monovalent ion concentration in mM here: '))
Mg = float(input('Input Mg concentration in mM here: '))
dNTPs = float(input('Input dNTP concentration in mM here: '))

gas_constant = 8.31445984848484848484 #J K-1-mol-1 
oligo_c = (oligo1_conc) * 1e-6 #Adjust concentration to mol/L
purine = ['A','G'] #Lists used to convert primer sequence into R/Y purine/pyrimidine nucleotide codes
pyrimidine = ['C','T']

primer1_length = len(primer1)
primer2_length = len(primer2)

#Captures the first and last base in each primer so that the energies associated with opening the helix can be calculated
p1_initial_base = primer1[0]
p1_terminal_base = primer1[-1]

p2_initial_base = primer2[0]
p2_terminal_base = primer2[-1]

#Calculation of each primer's GC content
def get_primer_gc(primer1, primer2):

    primer1_gc = (sum(1.0 for base in primer1 if base in ['G', 'C']) / primer1_length) * 100
    primer2_gc = (sum(1.0 for base in primer2 if base in ['G', 'C']) / primer2_length) * 100

    return primer1_gc, primer2_gc

primer1_gc, primer2_gc = get_primer_gc(primer1, primer2)

###################################################################################################################################
#Determines terminal base compensation parameters for each primer
###################################################################################################################################
'''
For relatively short DNA sequences, such as PCR primers, the initial melting of the helix is influenced by the bases at the ends.
DNA sequences around this length follow appoximately two-state melting. Therefore, the initial 'unzipping' plays a larger factor
in helix melting. Regardless of the underlying mechanism, G or C bases should require more energy to begin unzipping.
'''
def p1_end_comp(p1_initial_base, p1_terminal_base):

    AT_h, AT_s, GC_h, GC_s = 2.3, 4.1, .1, -2.8 #Compensation values for enthalpy and entropy in kcal/mol and eu
    dh_compensation, ds_compensation = 0, 0

    dh_init = (dh_compensation + AT_h if p1_initial_base in ['A', 'T'] else dh_compensation + GC_h)
    ds_init = (ds_compensation + AT_s if p1_initial_base in ['A', 'T'] else ds_compensation + GC_s)

    dh_term = (dh_compensation + AT_h if p1_terminal_base in ['A', 'T'] else dh_compensation + GC_h)
    ds_term = (ds_compensation + AT_s if p1_terminal_base in ['A', 'T'] else ds_compensation + GC_s)

    dh_compensation = dh_init + dh_term
    ds_compensation = ds_init + ds_term

    return dh_compensation, ds_compensation

def p2_end_comp(p2_initial_base,p2_terminal_base):

    AT_h, AT_s, GC_h, GC_s = 2.3, 4.1, .1, -2.8 #dh in kcal/mol, ds in eu
    dh_compensation, ds_compensation = 0, 0

    dh_init = (dh_compensation + AT_h if p2_initial_base in ['A', 'T'] else dh_compensation + GC_h)
    ds_init = (ds_compensation + AT_s if p2_initial_base in ['A', 'T'] else ds_compensation + GC_s)

    dh_term = (dh_compensation + AT_h if p2_terminal_base in ['A', 'T'] else dh_compensation + GC_h)
    ds_term = (ds_compensation + AT_s if p2_terminal_base in ['A', 'T'] else ds_compensation + GC_s)

    dh_compensation = dh_init + dh_term
    ds_compensation = ds_init + ds_term

    return dh_compensation, ds_compensation

p1_term_h, p1_term_s = p1_end_comp(p1_initial_base,p1_terminal_base)
p2_term_h, p2_term_s = p2_end_comp(p2_initial_base,p2_terminal_base)

###################################################################################################################################
#Melting Temperature Calculation and Salt Adjustments
###################################################################################################################################
'''
This section attempts to correct the model's accuracy discrepancy when computing the melting temperatures of oligos with 
different base sequences but the same GC content. Those types of oligos can differ by as much as 4C or more, however, 
Privalov and Crane-Robinson's model doesn't contain any way to account for those differences. I think the discrepancy 
may be due to the presence of alternating purine and pyrimidine bases. The alternating pairs may break up water binding
to adenine and thymine residues, thus decreasing the entropy cost associated with that water binding. In theory, the longer 
the alternating purine and pyrimidine chain, the more stable the helix is. Conversely, the more long strings of repeating
pyrimidine bases, the less stable the helix.
'''
RY_primer1_str = ''.join(['R' if base in purine else 'Y' for base in primer1])#Converts each primer into the nucleotide codes
RY_primer2_str = ''.join(['R' if base in purine else 'Y' for base in primer2])#R and Y for purine and pyrimidine bases

RY_pairs_p1 = re.findall('(?:RY){2,}', RY_primer1_str) #Generates a list of repeating purine and pyrimidine sequences
RY_pairs_p2 = re.findall('(?:RY){2,}', RY_primer2_str) #in each primer. Only includes sequences with two or more repeats.

AT_pairs_p1 = re.findall('(?:A|T){1}(?:A|T){1,}', primer1) #Generates a list of repeating adenine and thymine sequences
AT_pairs_p2 = re.findall('(?:A|T){1}(?:A|T){1,}', primer2) #in each primer. Only includes sequences with two or more repeats.

GC_pairs_p1 = re.findall('(?:G|C){1}(?:G|C){1,}', primer1) #Generates a list of repeating guanine and thymine cytosine
GC_pairs_p2 = re.findall('(?:G|C){1}(?:G|C){1,}', primer2) #in each primer. Only includes sequences with two or more repeats.

''' The above lists are needed so that their respective influences on helix stability can be measured and quantified, assuming
there is some influence. Don't yet know how best to do that. The functions below attempts to apply those differences to appropriately 
modify the melting temperature.
'''
#These pair length functions take the lists of pairs above and convert those captured sequences into lengths.
def p1_pair_lengths(AT_pairs_p1, GC_pairs_p1, RY_pairs_p1):

    AT_len, GC_len, RY_len = [], [], []

    for seq in AT_pairs_p1:
        AT_len.append(len(seq))
    for seq in GC_pairs_p1:
        GC_len.append(len(seq))
    for seq in RY_pairs_p1:
        RY_len.append(len(seq))
            
    return AT_len, GC_len, RY_len

def p2_pair_lengths(AT_pairs_p2, GC_pairs_p2, RY_pairs_p2):

    AT_len, GC_len, RY_len = [], [], []

    for seq in AT_pairs_p2:
        AT_len.append(len(seq))
    for seq in GC_pairs_p2:
        GC_len.append(len(seq))
    for seq in RY_pairs_p2:
        RY_len.append(len(seq))
            
    return AT_len, GC_len, RY_len

p1_AT_lengths, p1_GC_lengths, p1_RY_lengths = p1_pair_lengths(AT_pairs_p1, GC_pairs_p1, RY_pairs_p1)
p2_AT_lengths, p2_GC_lengths, p2_RY_lengths = p2_pair_lengths(AT_pairs_p2, GC_pairs_p2, RY_pairs_p2)

#Determines the primers' melting temperature based on the method developed by Privalov and Crane-Robinson
#See Privalov, P. L., & Crane-Robinson, C. (2018). https://doi.org/10.1016/j.pbiomolbio.2018.01.007
heat_capacity = .13 #kJ/K,mol-bp
H_AT_25, H_CG_25 = 25, 18.8 #kJ/mol-bp
S_AT_25, S_CG_25 = 72, 44.7 #J/K-mol-bp
delta_S_trans = gas_constant * math.log(2 / oligo_c)

def p1_melting_calculation(primer1, p1_term_h, p1_term_s):

    n_AT = primer1.count('A') + primer1.count('T') #Counts the number of each base in the primer
    n_CG = primer1.count('C') + primer1.count('G') #Used to calculate total enthalpy and entropy
    
    #Total enthalpy and entropy value calculations for determining the melting temperature. Takes the heat capacity into account.
    delta_H = ((H_AT_25 + (heat_capacity * (311.5 - 298.15))) * n_AT) + ((H_CG_25 + (heat_capacity * (311.5 - 298.15))) * n_CG) + p1_term_h
    delta_S = ((S_AT_25 + (heat_capacity * math.log(311.5 / 298.15))) * n_AT) + ((S_CG_25 + (heat_capacity * math.log((311.5) / 298.15))) * n_CG) + p1_term_s
    '''
The constraints here and below are subject to change as they are first guesses based on the limited data availale. As more
data is gathered, the length cut-offs and the correction values will likely change.
'''
    if len(primer1) > 15:
        for num in p1_AT_lengths:
            if num > 5:
                delta_S += 2.0
        for num in p1_GC_lengths:
            if num >= 5:
                delta_S -= .75
        for num in p1_AT_lengths:
            if num >= 4:
                delta_S += .5
        for num in p1_GC_lengths:
            if num >= 4:
                delta_S -= .7
        for num in p1_AT_lengths:
            if num >= 3:
                delta_S += .6
        for num in p1_GC_lengths:
            if num >= 2:
                delta_S -= .55
    else: #Loosens the length restrictions for short sequences
        for num in p1_AT_lengths:
            if num > 3:
                delta_S += 2.0
        for num in p1_GC_lengths:
            if num > 3:
                delta_S -= .75
        for num in p1_AT_lengths:
            if num > 2:
                delta_S += .5
        for num in p1_GC_lengths:
            if num >= 2:
                delta_S -= .6
        for num in p1_AT_lengths:
            if num >= 2:
                delta_S += .55

    tm = (1000 * delta_H) / (delta_S + (gas_constant * (math.log(2 / oligo_c)))) - 298.15

    return tm

def p2_melting_calculation(primer2, p2_term_h, p2_term_s):

    n_AT, n_CG = primer2.count('A') + primer2.count('T'), primer2.count('C') + primer2.count('G')

    delta_H = ((H_AT_25 + (heat_capacity * (311.5 - 298.15))) * n_AT) + ((H_CG_25 + (heat_capacity * (311.5 - 298.15))) * n_CG) + p2_term_h
    delta_S = ((S_AT_25 + (heat_capacity * math.log(311.5 / 298.15))) * n_AT) + ((S_CG_25 + (heat_capacity * math.log((311.5) / 298.15))) * n_CG) + p2_term_s

    if len(primer2) > 15:
        for num in p2_AT_lengths:
            if num > 5:
                delta_S += 2.0
        for num in p2_GC_lengths:
            if num >= 5:
                delta_S -= .75
        for num in p2_AT_lengths:
            if num >= 4:
                delta_S += .5
        for num in p2_GC_lengths:
            if num >= 4:
                delta_S -= .6
        for num in p2_AT_lengths:
            if num >= 3:
                delta_S += .6
        for num in p2_GC_lengths:
            if num >= 2:
                delta_S -= .55
    else:
        for num in p2_AT_lengths:
            if num > 3:
                delta_S += 2.0
        for num in p2_GC_lengths:
            if num > 3:
                delta_S -= .75
        for num in p2_AT_lengths:
            if num > 2:
                delta_S += .5
        for num in p2_GC_lengths:
            if num >= 2:
                delta_S -= .6
        for num in p2_AT_lengths:
            if num >= 2:
                delta_S += .55

    tm = (1000 * delta_H) / (delta_S + (gas_constant * (math.log(2 / oligo_c)))) - 298.15

    return tm

primer1_melting_temperature = p1_melting_calculation(primer1, p1_term_h, p1_term_s)
primer2_melting_temperature = p2_melting_calculation(primer2, p2_term_h, p2_term_s)

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

    a, b, c, d, e, f, g = 3.92e-5, -9.11e-6, 6.06e-5, 1.42e-5, -4.92e-4, 5.35e-4, 8.31e-5
    '''
The condition calculating the melting temperature if there is no monovalent salt present must be put first, before calculating R.
This is because the R equation requires dividing by the monovalent ion concentration. And of course, you can't divide by zero.
Additionally, although the condition of R > 6 uses the same equation as when no monovalent ions are present, the two conditions
can't be combined. In situations where no magnesium is present, the equation would try to take the natural log of zero, which
throws an error.
'''
    if mon == 0:

        salt2 = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        
        return (1 / salt2) - 273.15

    R = math.sqrt(mg) / mon

    if R < 0.22:
        
        salt1 = (1 / (primer1_melting_temperature + 273.15)) + ((4.59e-5 * (primer1_gc / 100)) - 2.90e-5) * math.log(mon) + 8.81e-6 * (math.log(mon)) ** 2
            
        return (1 / salt1) - 273.15

    elif R < 6.0:

        a = 3.92e-5 * (0.843 - (0.352 * math.sqrt(mon) * math.log(mon)))
        d = 1.42e-5 * ((1.279 - 4.03e-3 * math.log(mon)) - 8.03e-3 * (math.log(mon) ** 2))
        g = 8.31e-5 * ((0.486 - 0.258 * math.log(mon)) + 5.25e-3 * (math.log(mon) ** 3)) 

        salt2 = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)

        return (1 / salt2) - 273.15
    
    elif R > 6.0:

        salt2 = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)

        return (1 / salt2) - 273.15

def primer2_salt_correction(primer2_melting_temperature, primer2_gc, primer2_length, mg, mon):

    a, b, c, d, e, f, g = 3.92e-5, -9.11e-6, 6.06e-5, 1.42e-5, -4.92e-4, 5.35e-4, 8.31e-5
    
    if mon == 0:

        salt2 = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        
        return (1 / salt2) - 273.15

    R = math.sqrt(mg) / mon

    if R < 0.22:
        
        salt1 = (1 / (primer2_melting_temperature + 273.15)) + ((4.59e-5 * (primer2_gc / 100)) - 2.90e-5) * math.log(mon) + 8.81e-6 * (math.log(mon)) ** 2
            
        return (1 / salt1) - 273.15

    elif R < 6.0:

        a = 3.92e-5 * (0.843 - (0.352 * math.sqrt(mon) * math.log(mon)))
        d = 1.42e-5 * ((1.279 - 4.03e-3 * math.log(mon)) - 8.03e-3 * (math.log(mon) ** 2))
        g = 8.31e-5 * ((0.486 - 0.258 * math.log(mon)) + 5.25e-3 * (math.log(mon) ** 3)) 

        salt2 = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)

        return (1 / salt2) - 273.15
    
    elif R > 6.0:

        salt2 = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)

        return (1 / salt2) - 273.15

adj_primer1_melting_temperature = primer1_salt_correction(primer1_melting_temperature, primer1_gc, primer1_length, mg, mon)
adj_primer2_melting_temperature = primer2_salt_correction(primer2_melting_temperature, primer2_gc, primer2_length, mg, mon)

#Additional melting temperature adjustment to account for the fact that Privalov and Crane-Robinson's model accuracy can
#vary significantly depending on the GC content of the oligo. It tends to underestimate the Tm of low GC oligos and overestimates
#the Tm of high GC oligos. The highest accuracy is at a GC% around 50. Oligo length is also taken into account here.

def p1_gc_tm_adjustment(adj_primer1_melting_temperature, primer1_gc):

    '''
Correction factors for three GC% tiers. The correction values were calculated using the difference between the predicted melting
temperature using the standard model and the measured melting temperature from the 92 oligos in Owczarzy's 2004 paper. See
Owczarzy, R. et. al (2004). Biochemistry. https://doi.org/10.1021/bi034621r The differences were from around -6C for 20% GC 
to 14C for 80% GC. The melting temperature difference was then divided by the difference between the GC% where the model was most 
accurate, ~50%, and the GC% of the oligos at 20% and 80% GC. Since the melting temperature difference increased with oligo length, 
higher correction values are used to account for that length dependent effect.
'''
    low_gc, mid_gc, high_gc = .166667, .366667, .466667

    if primer1_gc < 45:
        #Multiplies the correction factor by the GC% difference between the primer and the approximate 'highest accuracy' GC%.
        p1_tm = float(adj_primer1_melting_temperature + (low_gc * (52 - primer1_gc)))

    elif primer1_gc >= 45 and primer1_gc < 60:
        if primer1_length >=25 and round(primer1_gc) in [47,48,49,50,51,52]:#The values here are arbitrary and are subject to change

            p1_tm = float(.959 * (adj_primer1_melting_temperature - (mid_gc * (primer1_gc - 50))))
        else:
            p1_tm = float((adj_primer1_melting_temperature - (mid_gc * (primer1_gc - 50))))

    else:
        if primer1_length <= 15: #Length cut-offs were based on the 92 oligos but are subject to change

            high_gc = .266667
            p1_tm = float(adj_primer1_melting_temperature - (high_gc * (primer1_gc - 48))) 
        elif primer1_length <= 25:
            
            p1_tm = float(adj_primer1_melting_temperature - (high_gc * (primer1_gc - 48)))
        else:

            high_gc = .586667
            p1_tm = float(adj_primer1_melting_temperature - (high_gc * (primer1_gc - 48))) 
    return p1_tm

def p2_gc_tm_adjustment(adj_primer2_melting_temperature, primer2_gc):

    low_gc, mid_gc, high_gc = .166667, .366667, .466667

    if primer2_gc < 45:

        p2_tm = float(adj_primer2_melting_temperature + (low_gc * (52 - primer2_gc)))

    elif primer2_gc >= 45 and primer2_gc < 60:
        if primer2_length >=25 and round(primer2_gc) in [47,48,49,50,51,52]:

            p2_tm = float(.959 * (adj_primer2_melting_temperature - (mid_gc * (primer2_gc - 50))))
        else:
            p2_tm = float((adj_primer2_melting_temperature - (mid_gc * (primer2_gc - 50))))

    else:
        if primer2_length <= 15:

            high_gc = .266667
            p2_tm = float(adj_primer2_melting_temperature - (high_gc * (primer2_gc - 48))) 

        elif primer2_length <= 25:
            
            p2_tm = float(adj_primer2_melting_temperature - (high_gc * (primer2_gc - 48)))

        else:

            high_gc = .586667
            p2_tm = float(adj_primer2_melting_temperature - (high_gc * (primer2_gc - 48))) 

    return p2_tm

p1_tm = p1_gc_tm_adjustment(adj_primer1_melting_temperature, primer1_gc)
p2_tm = p2_gc_tm_adjustment(adj_primer2_melting_temperature, primer2_gc)
###################################################################################################################################
#Printing the results
###################################################################################################################################
print('The GC content of primer 1 is ' + str(round(primer1_gc, 1)) + '%.' + 
' The GC content of primer 2 is ' + str(round(primer2_gc, 1)) + '%.')
print('The melting temperatures of your primers are ' + str(round(p1_tm, 1)) + 
'C for primer 1, and ' + str(round(p2_tm, 1)) + 'C for primer 2.')
