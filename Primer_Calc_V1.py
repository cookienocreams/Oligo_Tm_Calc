import re, math

#Lists with all non-overlapping NN pair possiblities
AA = ['AA', 'TT']
AC = ['AC', 'TG']
AG = ['AG', 'TC']
AT = ['AT']
CA = ['CA', 'GT']
CC = ['CC', 'GG']
CG = ['CG']
GA = ['GA', 'CT']
GC = ['GC']
TA = ['TA']

gas_constant = 1.98720425864083 #cal⋅K−1⋅mol−1
purine = ['A','G'] #Lists used to convert primer sequence into R/Y purine/pyrimidine nucleotide codes
pyrimidine = ['C','T']
#End Compensation parameters
YY_h,RY_h = 1.180248081487909, 0.9577358235882366 #h in kcal/mol, s in eu
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

primer1_length = len(primer1)
primer2_length = len(primer2)

#Captures the first and last base in each primer so that the energies associated with opening the helix can be calculated
p1_RY_seq_str = ''.join(['R' if base in purine else 'Y' for base in primer1])#Converts each primer into purine/pyrimidine strings
p2_RY_seq_str = ''.join(['R' if base in purine else 'Y' for base in primer2])

p1_initial_bases, p1_terminal_bases = p1_RY_seq_str[0:2], p1_RY_seq_str[-2:primer1_length]
p2_initial_bases, p2_terminal_bases = p2_RY_seq_str[0:2], p2_RY_seq_str[-2:primer2_length]

#Splits input primers into doublets for easy NN matching below
duo_primer1 = re.findall('.{1,2}', primer1) + re.findall('.{1,2}', primer1[1:])
duo_primer2 = re.findall('.{1,2}', primer2) + re.findall('.{1,2}', primer2[1:])

oligo_c = (oligo1_conc) * 1e-6 #Adjust concentration to mol/L

#Calculation of each primer's GC content
def get_primer_gc(primer1, primer2):

    primer1_GC = (sum((1.0 for base in primer1 if base in ['G', 'C'])) / primer1_length) * 100
    primer2_GC = (sum((1.0 for base in primer2 if base in ['G', 'C'])) / primer2_length) * 100

    return primer1_GC, primer2_GC

primer1_gc, primer2_gc = get_primer_gc(primer1, primer2)
###################################################################################################################################
#Chi-square Optimized Nearest Neighbor Parameters
###################################################################################################################################
AA_delta_s = -22.15461879548201
AA_delta_h = -7.909076567833574

AC_delta_s = -20.20944930135789
AC_delta_h = -7.564931774075489

AG_delta_s = -20.85324697548337
AG_delta_h = -7.866141135756189

AT_delta_s = -19.931105027077287
AT_delta_h = -7.114678926718001

CA_delta_s = -21.383943453258333
CA_delta_h = -7.814334448491746

CC_delta_s = -19.677205419525077
CC_delta_h = -7.780021121897234

CG_delta_s = -24.335240882791073
CG_delta_h = -9.344024883725648

GA_delta_s = -21.617116510641825
GA_delta_h = -7.934729350061229

GC_delta_s = -23.896038387305747
GC_delta_h = -9.399037969893644

TA_delta_s = -19.12903239340641
TA_delta_h = -6.580371705207274
###################################################################################################################################
#Determination of NN thermodynamic parameters for each primer
###################################################################################################################################

#Calculation of delta H total for each NN pair in primer 1
AA_h_p1 = AA_delta_h * sum(1.0 for NN in duo_primer1 if NN in AA)
AC_h_p1 = AC_delta_h * sum(1.0 for NN in duo_primer1 if NN in AC)
AG_h_p1 = AG_delta_h * sum(1.0 for NN in duo_primer1 if NN in AG)
AT_h_p1 = AT_delta_h * sum(1.0 for NN in duo_primer1 if NN in AT)
CA_h_p1 = CA_delta_h * sum(1.0 for NN in duo_primer1 if NN in CA)
CC_h_p1 = CC_delta_h * sum(1.0 for NN in duo_primer1 if NN in CC)
CG_h_p1 = CG_delta_h * sum(1.0 for NN in duo_primer1 if NN in CG)
GA_h_p1 = GA_delta_h * sum(1.0 for NN in duo_primer1 if NN in GA)
GC_h_p1 = GC_delta_h * sum(1.0 for NN in duo_primer1 if NN in GC)
TA_h_p1 = TA_delta_h * sum(1.0 for NN in duo_primer1 if NN in TA)

#Calculation of delta S total for each NN pair in primer 1
AA_s_p1 = AA_delta_s * sum(1.0 for NN in duo_primer1 if NN in AA)
AC_s_p1 = AC_delta_s * sum(1.0 for NN in duo_primer1 if NN in AC)
AG_s_p1 = AG_delta_s * sum(1.0 for NN in duo_primer1 if NN in AG)
AT_s_p1 = AT_delta_s * sum(1.0 for NN in duo_primer1 if NN in AT)
CA_s_p1 = CA_delta_s * sum(1.0 for NN in duo_primer1 if NN in CA)
CC_s_p1 = CC_delta_s * sum(1.0 for NN in duo_primer1 if NN in CC)
CG_s_p1 = CG_delta_s * sum(1.0 for NN in duo_primer1 if NN in CG)
GA_s_p1 = GA_delta_s * sum(1.0 for NN in duo_primer1 if NN in GA)
GC_s_p1 = GC_delta_s * sum(1.0 for NN in duo_primer1 if NN in GC)
TA_s_p1 = TA_delta_s * sum(1.0 for NN in duo_primer1 if NN in TA)

#Calculation of delta H total for each NN pair in primer 2
AA_h_p2 = AA_delta_h * sum(1.0 for NN in duo_primer2 if NN in AA)
AC_h_p2 = AC_delta_h * sum(1.0 for NN in duo_primer2 if NN in AC)
AG_h_p2 = AG_delta_h * sum(1.0 for NN in duo_primer2 if NN in AG)
AT_h_p2 = AT_delta_h * sum(1.0 for NN in duo_primer2 if NN in AT)
CA_h_p2 = CA_delta_h * sum(1.0 for NN in duo_primer2 if NN in CA)
CC_h_p2 = CC_delta_h * sum(1.0 for NN in duo_primer2 if NN in CC)
CG_h_p2 = CG_delta_h * sum(1.0 for NN in duo_primer2 if NN in CG)
GA_h_p2 = GA_delta_h * sum(1.0 for NN in duo_primer2 if NN in GA)
GC_h_p2 = GC_delta_h * sum(1.0 for NN in duo_primer2 if NN in GC)
TA_h_p2 = TA_delta_h * sum(1.0 for NN in duo_primer2 if NN in TA)

#Calculation of delta S total for each NN pair in primer 2
AA_s_p2 = AA_delta_s * sum(1.0 for NN in duo_primer2 if NN in AA)
AC_s_p2 = AC_delta_s * sum(1.0 for NN in duo_primer2 if NN in AC)
AG_s_p2 = AG_delta_s * sum(1.0 for NN in duo_primer2 if NN in AG)
AT_s_p2 = AT_delta_s * sum(1.0 for NN in duo_primer2 if NN in AT)
CA_s_p2 = CA_delta_s * sum(1.0 for NN in duo_primer2 if NN in CA)
CC_s_p2 = CC_delta_s * sum(1.0 for NN in duo_primer2 if NN in CC)
CG_s_p2 = CG_delta_s * sum(1.0 for NN in duo_primer2 if NN in CG)
GA_s_p2 = GA_delta_s * sum(1.0 for NN in duo_primer2 if NN in GA)
GC_s_p2 = GC_delta_s * sum(1.0 for NN in duo_primer2 if NN in GC)
TA_s_p2 = TA_delta_s * sum(1.0 for NN in duo_primer2 if NN in TA)

sum_p1_delta_h = AA_h_p1 + AC_h_p1 + AG_h_p1 + AT_h_p1 + CA_h_p1 + CC_h_p1 + CG_h_p1 + GA_h_p1 + GC_h_p1 + TA_h_p1
sum_p1_delta_s = AA_s_p1 + AC_s_p1 + AG_s_p1 + AT_s_p1 + CA_s_p1 + CC_s_p1 + CG_s_p1 + GA_s_p1 + GC_s_p1 + TA_s_p1

sum_p2_delta_h = AA_h_p2 + AC_h_p2 + AG_h_p2 + AT_h_p2 + CA_h_p2 + CC_h_p2 + CG_h_p2 + GA_h_p2 + GC_h_p2 + TA_h_p2
sum_p2_delta_s = AA_s_p2 + AC_s_p2 + AG_s_p2 + AT_s_p2 + CA_s_p2 + CC_s_p2 + CG_s_p2 + GA_s_p2 + GC_s_p2 + TA_s_p2

#Calculates the enthalpic compensation for the primer ends and adds it to the total enthalpy
def p1_end_comp(p1_initial_bases,p1_terminal_bases,sum_p1_delta_h):
    if p1_initial_bases in ['YY', 'RR']:
        sum_p1_delta_h += YY_h
    if p1_initial_bases in ['RY', 'YR']:
        sum_p1_delta_h += RY_h
    if p1_terminal_bases in ['YY', 'RR']:
        sum_p1_delta_h += YY_h
    if p1_terminal_bases in ['RY', 'YR']:
        sum_p1_delta_h += RY_h
    return sum_p1_delta_h

def p2_end_comp(p2_initial_bases,p2_terminal_bases,sum_p2_delta_h):
    if p2_initial_bases in ['YY', 'RR']:
        sum_p2_delta_h += YY_h
    if p2_initial_bases in ['RY', 'YR']:
        sum_p2_delta_h += RY_h
    if p2_terminal_bases in ['YY', 'RR']:
        sum_p2_delta_h += YY_h
    if p2_terminal_bases in ['RY', 'YR']:
        sum_p2_delta_h += RY_h
    return sum_p2_delta_h

p1_total_h = p1_end_comp(p1_initial_bases,p1_terminal_bases,sum_p1_delta_h)
p2_total_h = p2_end_comp(p2_initial_bases,p2_terminal_bases,sum_p2_delta_h)
###################################################################################################################################
#Melting Temperature Calculation and Salt Adjustments
###################################################################################################################################

#Determines the melting temperature of the primers
primer1_melting_temperature = (1000 * p1_total_h) / (sum_p1_delta_s + (gas_constant * (math.log(oligo_c)))) - 273.15
primer2_melting_temperature = (1000 * p2_total_h) / (sum_p2_delta_s + (gas_constant * (math.log(oligo_c)))) - 273.15

#Adjustments and unit conversions for the chosen buffer conditions
Mon = Mono / 2.0 #Divide by two to account for the counterion present, e.g. Cl-, SO4-, etc.
mg_adj = Mg * 1e-3 #Converts to mol/L
mon = Mon * 1e-3
dntps = dNTPs * 1e-3 
ka = 3e4 #Association constant for the Mg2+--dNTP complex. Used to calculate the free magnesium in the buffer
mg = (-(ka * dntps - ka * mg_adj + 1.0) + math.sqrt((ka * dntps - ka * mg_adj + 1.0) ** 2 + 4.0 * ka * mg_adj)) / (2.0 * ka)

#Equations from Owczarzy 2008 that adjusts melting temperatures according to monovalent ion, magnesium, and dNTP concentrations
def primer1_salt_correction(primer1_melting_temperature, primer1_gc, primer1_length, mg, mon):

    a, b, c, d, e, f, g = 3.92e-5, -9.11e-6, 6.06e-5, 1.42e-5, -4.82e-4, 5.65e-4, 8.31e-5
    
    if mon == 0:
        salt2 = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        return round((1 / salt2) - 273.15, 1)

    R = math.sqrt(mg) / mon 

    if R < 0.22: #Constant values here were also optimized used Chi-square minimization
        salt1 = (1 / (primer1_melting_temperature + 273.15)) + ((3.087315696964695e-05 * (primer1_gc / 100)) - 2.2560798983637505e-05) * math.log(mon) + 8.081965258030895e-06 * (math.log(mon)) ** 2
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
        salt1 = (1 / (primer2_melting_temperature + 273.15)) + ((3.087315696964695e-05 * (primer2_gc / 100)) - 2.2560798983637505e-05) * math.log(mon) + 8.081965258030895e-06 * (math.log(mon)) ** 2
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
print('The melting temperatures of your primers are ' + str(round(primer1_melting_temperature, 1)) + 
'C for primer 1, and ' + str(round(primer2_melting_temperature, 1)) + 'C for primer 2.')
print('The salt adjusted melting temperatures of your primers are ' + str(adj_primer1_melting_temperature) + 
'C for primer 1, and ' + str(adj_primer2_melting_temperature) + 'C for primer 2.')
