import math

###################################################################################################################################
#Prompts to get the information required to begin the analysis
#Primer information is input here
#The primer sequences will be analyzed first by calculating their GC content
#Analysis of primer melting temperatures occurs last
###################################################################################################################################

primer1 = input('Input primer 1 sequence here: ')
primer2 = input('Input primer 2 sequence here: ')
oligo1_conc = float(input('Input total single strand oligo concentration in uM here: '))
Mono = float(input('Input total monovalent ion concentration in mM here: '))
Mg = float(input('Input Mg concentration in mM here: '))
dNTPs = float(input('Input dNTP concentration in mM here: '))

gas_constant = 8.31445984848484848484 #J K-1-mol-1 

primer1_length = len(primer1)
primer2_length = len(primer2)

p1_terminal_base = primer1[-1]
p1_initial_base = primer1[0]

p2_terminal_base = primer2[-1]
p2_initial_base = primer2[0]

oligo_c = (oligo1_conc) * 1e-6

#Calculation of each primer's GC content
def get_primer_gc(primer1, primer2):

    primer1_GC = (sum((1.0 for base in primer1 if base in ['G', 'C'])) / primer1_length) * 100
    primer2_GC = (sum((1.0 for base in primer2 if base in ['G', 'C'])) / primer2_length) * 100

    return primer1_GC, primer2_GC

primer1_gc, primer2_gc = get_primer_gc(primer1, primer2)

###################################################################################################################################
#Determine terminal base compensation parameters for each primer
###################################################################################################################################

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
#Melting Temperature Adjustments
###################################################################################################################################

#Determines the melting temperature of the primers based on the method developed by Privalov and Crane-Robinson
#See Privalov, P. L., & Crane-Robinson, C. (2018). https://doi.org/10.1016/j.pbiomolbio.2018.01.007
def p1_melting_calculation(primer1, p1_term_h, p1_term_s):
    
    heat_capacity = .13 #kJ/K,mol-bp
    H_A_25, H_T_25 = 25.18, 23 #kJ/mol-bp
    S_A_25, S_T_25 = 71.81, 65 #J/K-mol-bp
    H_CG_25 = 18.5
    S_CG_25 = 44.9
    delta_S_trans = gas_constant * math.log(2 / oligo_c)
    n_A, n_T = primer1.count('A'), primer1.count('T')
    n_CG = primer1.count('C') + primer1.count('G')

    delta_H = ((H_A_25 + (heat_capacity * (311.5 - 298.15))) * n_A) + ((H_T_25 + (heat_capacity * (311.5 - 298.15))) * n_T) + ((H_CG_25 + (heat_capacity * (311.5 - 298.15))) * n_CG) + p1_term_h
    delta_S = ((S_A_25 + (heat_capacity * math.log(311.5 / 298.15))) * n_A) + ((S_T_25 + (heat_capacity * math.log(311.5 / 298.15))) * n_T) + ((S_CG_25 + (heat_capacity * math.log((311.5) / 298.15))) * n_CG) + p1_term_s

    tm = (1000 * delta_H) / (delta_S + (gas_constant * (math.log(2 / oligo_c)))) - 298.15

    return tm

def p2_melting_calculation(primer2, p2_term_h, p2_term_s):
    
    heat_capacity = .13 #kJ/K,mol-bp
    H_A_25, H_T_25 = 25.18, 23 #kJ/mol-bp
    S_A_25, S_T_25 = 71.81, 65 #J/K-mol-bp
    H_CG_25 = 18.5
    S_CG_25 = 44.9
    delta_S_trans = gas_constant * math.log(2 / oligo_c)
    n_A, n_T = primer2.count('A'), primer2.count('T')
    n_CG = primer2.count('C') + primer2.count('G')

    delta_H = ((H_A_25 + (heat_capacity * (311.5 - 298.15))) * n_A) + ((H_T_25 + (heat_capacity * (311.5 - 298.15))) * n_T) + ((H_CG_25 + (heat_capacity * (311.5 - 298.15))) * n_CG) + p2_term_h
    delta_S = ((S_A_25 + (heat_capacity * math.log(311.5 / 298.15))) * n_A) + ((S_T_25 + (heat_capacity * math.log(311.5 / 298.15))) * n_T) + ((S_CG_25 + (heat_capacity * math.log((311.5) / 298.15))) * n_CG) + p2_term_s

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

    a, b, c, d, e, f, g = 3.92e-5, -9.11e-6, 6.06e-5, 1.42e-5, -4.82e-4, 5.65e-4, 8.31e-5
    
    if mon == 0:

        salt2 = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)

        return round((1 / salt2) - 273.15, 1)

    R = math.sqrt(mg) / mon 

    if R < 0.22:

        salt1 = (1 / (primer1_melting_temperature + 273.15)) + ((4.59e-5 * (primer1_gc / 100)) - 2.90e-5) * math.log(mon) + 8.81e-6 * (math.log(mon)) ** 2

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

    a, b, c, d, e, f, g = 3.92e-5, -9.11e-6, 6.06e-5, 1.42e-5, -4.82e-4, 5.65e-4, 8.31e-5
    
    if mon == 0:

        salt2 = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)

        return round((1 / salt2) - 273.15, 1)

    R = math.sqrt(mg) / mon 

    if R < 0.22:

        salt1 = (1 / (primer2_melting_temperature + 273.15)) + ((4.59e-5 * (primer2_gc / 100)) - 2.90e-5) * math.log(mon) + 8.81e-6 * (math.log(mon)) ** 2

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