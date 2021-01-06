import re, math

#Lists with all non-overlapping NN pair possiblities
NN_pairs_list = [['AA', 'TT'], ['AC', 'TG'], ['AG', 'TC'], ['AT'], ['CA', 'GT'], ['CC', 'GG'], ['CG'], ['GA', 'CT'], ['GC'], ['TA']]
gas_constant = 1.98720425864083 #cal⋅K−1⋅mol−1
purine, pyrimidine = ['A','G'], ['C','T'] #Lists used to convert primer sequence into R/Y purine/pyrimidine nucleotide codes
###########################################################################################################################################
#Prompts to get the information required to begin the analysis
#Primer information is input here
#The primer sequences will be analyzed first by calculating their GC content
#Analysis of primer melting temperatures occurs last
###########################################################################################################################################
primer1 = input('Input primer 1 sequence here: ')
primer2 = input('Input primer 2 sequence here: ')
oligo1_conc = float(input('Input total single strand oligo concentration in uM here: '))
Mono = float(input('Input total monovalent ion concentration in mM here: '))
Mg = float(input('Input Mg concentration in mM here: '))
dNTPs = float(input('Input dNTP concentration in mM here: '))

primer1_length = len(primer1)
primer2_length = len(primer2)

#Converts primers into their purine/pyrimidine nucleotide codes
primer1_RY_sequence_str = ''.join(['R' if base in purine else 'Y' for base in primer1])
primer2_RY_sequence_str = ''.join(['R' if base in purine else 'Y' for base in primer2])
#Captures the first and last pair of bases in each primer so that the energies associated with opening the helix can be calculated
p1_initial_bases, p1_terminal_bases = primer1_RY_sequence_str[0:2], primer1_RY_sequence_str[-2:primer1_length]
p2_initial_bases, p2_terminal_bases = primer2_RY_sequence_str[0:2], primer2_RY_sequence_str[-2:primer2_length]

#Splits input primers into doublets for easy NN matching below
primer1_NN_list = re.findall('.{1,2}', primer1) + re.findall('.{1,2}', primer1[1:])
primer2_NN_list = re.findall('.{1,2}', primer2) + re.findall('.{1,2}', primer2[1:])

oligo_c = (oligo1_conc) * 1e-6 #Adjusts primer concentration to mol/L

#Calculation of each primer's GC content
primer1_gc = (sum((1.0 for base in primer1 if base in ['G', 'C'])) / primer1_length) * 100
primer2_gc = (sum((1.0 for base in primer2 if base in ['G', 'C'])) / primer2_length) * 100
###########################################################################################################################################
#Modified Nearest Neighbor Parameters based on SantaLucia's 1998 paper
###########################################################################################################################################
AA_delta_s, AA_delta_h = -22.15461879548201, -7.972759800039155
AC_delta_s, AC_delta_h = -20.20944930135789, -7.6063649173383086
AG_delta_s, AG_delta_h = -20.85324697548337, -7.7250305886539525
AT_delta_s, AT_delta_h = -19.931105027077287, -7.069032819297647
CA_delta_s, CA_delta_h = -21.383943453258333, -7.9794539347093165
CC_delta_s, CC_delta_h = -19.677205419525077, -7.612190594454272
CG_delta_s, CG_delta_h = -24.335240882791073, -9.348781072782547
GA_delta_s, GA_delta_h = -21.617116510641825, -8.058362530689527
GC_delta_s, GC_delta_h = -23.896038387305747, -9.457805766372232
TA_delta_s, TA_delta_h = -19.12903239340641, -6.6574837353237415
#Lists of each nearest neighbor's enthalpy and entropy values that will be used to calculate the thermodynamic values of each primer
Delta_S = [AA_delta_s,AC_delta_s,AG_delta_s,AT_delta_s,CA_delta_s,CC_delta_s,CG_delta_s,GA_delta_s,GC_delta_s,TA_delta_s]
Delta_H = [AA_delta_h,AC_delta_h,AG_delta_h,AT_delta_h,CA_delta_h,CC_delta_h,CG_delta_h,GA_delta_h,GC_delta_h,TA_delta_h]
#End Compensation parameters adjusted for primer length and changes in monovalent ion concentration
YY_constant_one = (-.0235 * primer1_length) + .1273
YY_constant_two = (.1639 * primer1_length) - .895
YR_constant_one = (-.296 * math.log(primer1_length)) + .5058
YR_constant_two = (2.0303 * math.log(primer1_length)) - 3.4594
YY_length_adjust_eq = YY_constant_one * math.log(Mono) + YY_constant_two
YR_length_adjust_eq = YR_constant_one * math.log(Mono) + YR_constant_two
YY_h,RY_h = 0.959547884222969 - YY_length_adjust_eq, 0.849386567650066 - YR_length_adjust_eq
###########################################################################################################################################
#Determination of the NN thermodynamic parameters for each primer
###########################################################################################################################################
#Calculation of delta H and delta S total for both primers
def p1_thermodynamic_sum(Delta_H,Delta_S,NN_pairs_list):
    count, p1_total_dh, p1_total_ds, p2_total_dh, p2_total_ds = 0, 0, 0, 0, 0 #Needed to track values of each NN pair
    NN_pair = NN_pairs_list[count] #Indexed to calculate the values of each set of NN in the nearest-neighbor list
    for NN_pair in NN_pairs_list:
        p1_total_dh += Delta_H[count] * sum((1.0 for NN in primer1_NN_list if NN in NN_pair))#Calculates the number of occurances of each NN in 
        p1_total_ds += Delta_S[count] * sum((1.0 for NN in primer1_NN_list if NN in NN_pair))#the primer sequence and multiplies the number of
        p2_total_dh += Delta_H[count] * sum((1.0 for NN in primer2_NN_list if NN in NN_pair))#occurances by the associated delta h/s value
        p2_total_ds += Delta_S[count] * sum((1.0 for NN in primer2_NN_list if NN in NN_pair))
        count += 1 #Adds one to the count to keep track of which NNs have been calculated already
    return p1_total_dh, p1_total_ds, p2_total_dh, p2_total_ds
p1_dh_calc, p1_ds_calc, p2_dh_calc, p2_ds_calc = p1_thermodynamic_sum(Delta_H,Delta_S,NN_pairs_list)

#Adds in the enthalpic compensation for the ends of each primer
p1_dh_init = (p1_dh_calc + YY_h if p1_initial_bases in ['YY','RR'] else p1_dh_calc + RY_h)
p1_dh_total = (p1_dh_init + YY_h if p1_terminal_bases in ['YY','RR'] else p1_dh_init + RY_h)

p2_dh_init = (p2_dh_calc + YY_h if p2_initial_bases in ['YY','RR'] else p2_dh_calc + RY_h)
p2_dh_total = (p2_dh_init + YY_h if p2_terminal_bases in ['YY','RR'] else p2_dh_init + RY_h)
###########################################################################################################################################
#Melting Temperature Calculation and Salt Adjustments
###########################################################################################################################################
#Determines the melting temperature of the primers
primer1_melting_temperature = (1000 * p1_dh_total) / (p1_ds_calc + (gas_constant * (math.log(oligo_c)))) - 273.15
primer2_melting_temperature = (1000 * p2_dh_total) / (p2_ds_calc + (gas_constant * (math.log(oligo_c)))) - 273.15

#Adjustments and unit conversions for the chosen buffer conditions
Mon = Mono / 2.0 #Divide by two to account for the counterion present, e.g. Cl-, SO4-, etc.
mg_adj = Mg * 1e-3 #Converts to mol/L
mon = Mon * 1e-3
dntps = dNTPs * 1e-3 
ka = 3e4 #Association constant for the Mg2+--dNTP complex. Used to calculate the free magnesium in the buffer
mg = (-(ka * dntps - ka * mg_adj + 1.0) + math.sqrt((ka * dntps - ka * mg_adj + 1.0) ** 2 + 4.0 * ka * mg_adj)) / (2.0 * ka)

#Equations from Owczarzy 2008 that adjusts melting temperatures according to monovalent ion, magnesium, and dNTP concentrations
def primer1_salt_correction(primer1_melting_temperature, primer1_gc, primer1_length, mg, mon):
    #General constants and constants to regulate balance between Mg2+ and Na+ ions
    a, b, c, d, e, f, g = 3.919e-5, -2.88e-5, 3.603e-5, 2.322e-5, -3.507e-4, 4.711e-4, 6.52e-5 
    const_a,const_b,const_c,const_d,const_e,const_f,const_g,const_h = -0.1156,-2.0725,-0.1445,6.247e-3,6.131e-3,0.0314,0.5308,4.563e-3
    
    if mon == 0:
        salt2 = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        return round((1 / salt2) - 273.15, 1) #This equation converts the above calculation back into a useable melting temperature
    R = math.sqrt(mg) / mon #R is used to determine whether salt or magnesium is the primary factor affecting melting temperature
    if R < 0.22: #This equation is used if salt is the primary factor
        if Mono > 900:
            salt1 = (1 / (primer1_melting_temperature + 273.15)) + ((4.29e-5 * (primer1_gc / 100)) - 3.95e-5) * math.log(mon) + 9.40e-6 * (math.log(mon)) ** 2
            return round((primer1_melting_temperature + ((1 / salt1) - 273.15)) / 2, 1)
        salt1 = (1 / (primer1_melting_temperature + 273.15)) + ((4.29e-5 * (primer1_gc / 100)) - 3.95e-5) * math.log(mon) + 9.40e-6 * (math.log(mon)) ** 2
        return round((1 / salt1) - 273.15, 1)
    #This equation is used if there is a complex balance between salt and magnesium
    elif R < 6.0:
        a = 3.92e-5 * (const_a * math.log(mg) - (const_b * math.sqrt(mon) * math.log(mon)))
        d = 2.32e-5 * ((const_c * math.log(mg) - const_d * math.log(mon)) - const_e * (math.log(mon) ** 2))
        g = 6.52e-5 * ((const_f * math.log(mg) - const_g * math.log(mon)) + const_h * (math.log(mon) ** 3))
        salt2 = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        return round((1 / salt2) - 273.15, 1)
    #This equation is used if magnesium is the primary factor
    elif R > 6.0:
        salt2 = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        return round((1 / salt2) - 273.15,1)

def primer2_salt_correction(primer2_melting_temperature, primer2_gc, primer2_length, mg, mon):
    a, b, c, d, e, f, g = 3.919e-5, -2.88e-5, 3.603e-5, 2.322e-5, -3.507e-4, 4.711e-4, 6.52e-5
    const_a,const_b,const_c,const_d,const_e,const_f,const_g,const_h = -0.1156,-2.0725,-0.1445,6.247e-3,6.131e-3,0.0314,0.5308,4.563e-3
    
    if mon == 0:
        salt2 = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        return round((1 / salt2) - 273.15, 1)
    R = math.sqrt(mg) / mon 
    if R < 0.22:
        if Mono > 900:
            salt1 = (1 / (primer2_melting_temperature + 273.15)) + ((4.29e-5 * (primer2_gc / 100)) - 3.95e-5) * math.log(mon) + 9.40e-6 * (math.log(mon)) ** 2
            return round((primer2_melting_temperature + ((1 / salt1) - 273.15)) / 2, 1)
        salt1 = (1 / (primer2_melting_temperature + 273.15)) + ((4.29e-5 * (primer2_gc / 100)) - 3.95e-5) * math.log(mon) + 9.40e-6 * (math.log(mon)) ** 2
        return round((1 / salt1) - 273.15, 1)
    elif R < 6.0:
        a = 3.92e-5 * (const_a * math.log(mg) - (const_b * math.sqrt(mon) * math.log(mon)))
        d = 2.32e-5 * ((const_c * math.log(mg) - const_d * math.log(mon)) - const_e * (math.log(mon) ** 2))
        g = 6.52e-5 * ((const_f * math.log(mg) - const_g * math.log(mon)) + const_h * (math.log(mon) ** 3)) 
        salt2 = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        return round((1 / salt2) - 273.15, 1)
    elif R > 6.0:
        salt2 = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        return round((1 / salt2) - 273.15,1)

salt_adj_primer1_melting_temperature = primer1_salt_correction(primer1_melting_temperature, primer1_gc, primer1_length, mg, mon)
salt_adj_primer2_melting_temperature = primer2_salt_correction(primer2_melting_temperature, primer2_gc, primer2_length, mg, mon)
###########################################################################################################################################
#Printing the results
###########################################################################################################################################
print('The GC content of primer 1 is ' + str(round(primer1_gc, 1)) + '%.' + ' The GC content of primer 2 is ' + str(round(primer2_gc, 1)) + '%.')
print('The melting temperatures of your primers are ' + str(salt_adj_primer1_melting_temperature) + 'C for primer 1, and ' + str(salt_adj_primer2_melting_temperature) + 'C for primer 2.')
