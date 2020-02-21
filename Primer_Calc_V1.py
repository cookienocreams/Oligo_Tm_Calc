import re, math

###################################################################################################################################
#Modified NN Parameters based on SantaLucia's 1998 paper
###################################################################################################################################

AA_delta_s = -22.15
AA_delta_h = -7.95
AA_delta_g = AA_delta_h - (273.15*(AA_delta_s/1000))

AC_delta_s = -20.425
AC_delta_h = -7.915
AC_delta_g = AC_delta_h - (273.15*(AC_delta_s/1000))

AG_delta_s = -21.125
AG_delta_h = -8.11
AG_delta_g = AG_delta_h - (273.15*(AG_delta_s/1000))

AT_delta_s = -20.175
AT_delta_h = -7.95
AT_delta_g = AT_delta_h - (273.15*(AT_delta_s/1000))

CA_delta_s = -22.27
CA_delta_h = -8.5
CA_delta_g = CA_delta_h - (273.15*(CA_delta_s/1000))

CC_delta_s = -20.26
CC_delta_h = -8.12
CC_delta_g = CC_delta_h - (273.15*(CC_delta_s/1000))

CG_delta_s = -24.95
CG_delta_h = -10.2
CG_delta_g = CG_delta_h - (273.15*(CG_delta_s/1000))

GA_delta_s = -22.45
GA_delta_h = -8.1
GA_delta_g = GA_delta_h - (273.15*(GA_delta_s/1000))

GC_delta_s = -24.7
GC_delta_h = -9.78
GC_delta_g = GC_delta_h - (273.15*(GC_delta_s/1000))

TA_delta_s = -20.175
TA_delta_h = -7.34
TA_delta_g = TA_delta_h - (273.15*(TA_delta_s/1000))

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

#Makes each input primer a list so that each nucleotide can be counted individually
lst_primer1 = list(primer1)
lst_primer2 = list(primer2)

#Splits input primers into doublets for easy NN matching below
duo_primer1 = re.findall('.{1,2}', primer1) + re.findall('.{1,2}', primer1[1:])
duo_primer2 = re.findall('.{1,2}', primer2) + re.findall('.{1,2}', primer2[1:])

primer1_length = len(primer1)
primer2_length = len(primer2)

p1_terminal_base = lst_primer1[-1]
p1_initial_base = lst_primer1[0]

p2_terminal_base = lst_primer2[-1]
p2_initial_base = lst_primer2[0]

oligo_c = (oligo1_conc / 2.0) * 1e-9

#Calculation of each primer's GC content
def get_primer1_gc(primer1):

    primer1_GC = (sum([1.0 for i in lst_primer1 if i in ['G', 'C']]) / primer1_length) * 100

    return primer1_GC

def get_primer2_gc(primer2):

    primer2_GC = (sum([1.0 for i in lst_primer2 if i in ['G', 'C']]) / primer2_length) * 100

    return primer2_GC

###################################################################################################################################
#Determination of NN thermodynamic parameters for each primer
###################################################################################################################################

#Calculation of delta H total for each NN pair in primer 1
AA_h_p1 = AA_delta_h * sum([1.0 for i in duo_primer1 if i in AA])
AC_h_p1 = AC_delta_h * sum([1.0 for i in duo_primer1 if i in AC])
AG_h_p1 = AG_delta_h * sum([1.0 for i in duo_primer1 if i in AG])
AT_h_p1 = AT_delta_h * sum([1.0 for i in duo_primer1 if i in AT])
CA_h_p1 = CA_delta_h * sum([1.0 for i in duo_primer1 if i in CA])
CC_h_p1 = CC_delta_h * sum([1.0 for i in duo_primer1 if i in CC])
CG_h_p1 = CG_delta_h * sum([1.0 for i in duo_primer1 if i in CG])
GA_h_p1 = GA_delta_h * sum([1.0 for i in duo_primer1 if i in GA])
GC_h_p1 = GC_delta_h * sum([1.0 for i in duo_primer1 if i in GC])
TA_h_p1 = TA_delta_h * sum([1.0 for i in duo_primer1 if i in TA])

#Calculation of delta S total for each NN pair in primer 1
AA_s_p1 = AA_delta_s * sum([1.0 for i in duo_primer1 if i in AA])
AC_s_p1 = AC_delta_s * sum([1.0 for i in duo_primer1 if i in AC])
AG_s_p1 = AG_delta_s * sum([1.0 for i in duo_primer1 if i in AG])
AT_s_p1 = AT_delta_s * sum([1.0 for i in duo_primer1 if i in AT])
CA_s_p1 = CA_delta_s * sum([1.0 for i in duo_primer1 if i in CA])
CC_s_p1 = CC_delta_s * sum([1.0 for i in duo_primer1 if i in CC])
CG_s_p1 = CG_delta_s * sum([1.0 for i in duo_primer1 if i in CG])
GA_s_p1 = GA_delta_s * sum([1.0 for i in duo_primer1 if i in GA])
GC_s_p1 = GC_delta_s * sum([1.0 for i in duo_primer1 if i in GC])
TA_s_p1 = TA_delta_s * sum([1.0 for i in duo_primer1 if i in TA])

#Calculation of delta H total for each NN pair in primer 2
AA_h_p2 = AA_delta_h * sum([1.0 for i in duo_primer2 if i in AA])
AC_h_p2 = AC_delta_h * sum([1.0 for i in duo_primer2 if i in AC])
AG_h_p2 = AG_delta_h * sum([1.0 for i in duo_primer2 if i in AG])
AT_h_p2 = AT_delta_h * sum([1.0 for i in duo_primer2 if i in AT])
CA_h_p2 = CA_delta_h * sum([1.0 for i in duo_primer2 if i in CA])
CC_h_p2 = CC_delta_h * sum([1.0 for i in duo_primer2 if i in CC])
CG_h_p2 = CG_delta_h * sum([1.0 for i in duo_primer2 if i in CG])
GA_h_p2 = GA_delta_h * sum([1.0 for i in duo_primer2 if i in GA])
GC_h_p2 = GC_delta_h * sum([1.0 for i in duo_primer2 if i in GC])
TA_h_p2 = TA_delta_h * sum([1.0 for i in duo_primer2 if i in TA])

#Calculation of delta S total for each NN pair in primer 2
AA_s_p2 = AA_delta_s * sum([1.0 for i in duo_primer2 if i in AA])
AC_s_p2 = AC_delta_s * sum([1.0 for i in duo_primer2 if i in AC])
AG_s_p2 = AG_delta_s * sum([1.0 for i in duo_primer2 if i in AG])
AT_s_p2 = AT_delta_s * sum([1.0 for i in duo_primer2 if i in AT])
CA_s_p2 = CA_delta_s * sum([1.0 for i in duo_primer2 if i in CA])
CC_s_p2 = CC_delta_s * sum([1.0 for i in duo_primer2 if i in CC])
CG_s_p2 = CG_delta_s * sum([1.0 for i in duo_primer2 if i in CG])
GA_s_p2 = GA_delta_s * sum([1.0 for i in duo_primer2 if i in GA])
GC_s_p2 = GC_delta_s * sum([1.0 for i in duo_primer2 if i in GC])
TA_s_p2 = TA_delta_s * sum([1.0 for i in duo_primer2 if i in TA])

sum_p1_delta_h = AA_h_p1 + AC_h_p1 + AG_h_p1 + AT_h_p1 + CA_h_p1 + CC_h_p1 + CG_h_p1 + GA_h_p1 + GC_h_p1 + TA_h_p1
sum_p1_delta_s = AA_s_p1 + AC_s_p1 + AG_s_p1 + AT_s_p1 + CA_s_p1 + CC_s_p1 + CG_s_p1 + GA_s_p1 + GC_s_p1 + TA_s_p1

sum_p2_delta_h = AA_h_p2 + AC_h_p2 + AG_h_p2 + AT_h_p2 + CA_h_p2 + CC_h_p2 + CG_h_p2 + GA_h_p2 + GC_h_p2 + TA_h_p2
sum_p2_delta_s = AA_s_p2 + AC_s_p2 + AG_s_p2 + AT_s_p2 + CA_s_p2 + CC_s_p2 + CG_s_p2 + GA_s_p2 + GC_s_p2 + TA_s_p2

def p1_terminal_comp(p1_initial_base, p1_terminal_base, sum_p1_delta_h, sum_p1_delta_s):

    AT_h, AT_s, GC_h, GC_s = 2.3, 5.9, .25, -.7

    for base in p1_initial_base:   
        if base == 'A' or 'T':

            sum_p1_delta_h += AT_h
            sum_p1_delta_s += AT_s

            return sum_p1_delta_h, sum_p1_delta_s

    for base in p1_initial_base:
        if base == 'G' or 'C':

            sum_p1_delta_h += GC_h
            sum_p1_delta_s += GC_s

            return sum_p1_delta_h, sum_p1_delta_s

    for base in p1_terminal_base:
        if base == 'A' or 'T':

            sum_p1_delta_h += AT_h
            sum_p1_delta_s += AT_s

            return sum_p1_delta_h, sum_p1_delta_s

    for base in p1_terminal_base:
        if base == 'G' or 'C':

            sum_p1_delta_h += GC_h
            sum_p1_delta_s += GC_s

            return sum_p1_delta_h, sum_p1_delta_s

    return sum_p1_delta_h, sum_p1_delta_s

def p2_terminal_comp(p2_initial_base,p2_terminal_base, sum_p2_delta_h, sum_p2_delta_s):

    AT_h, AT_s, GC_h, GC_s = 2.3, 5.9, .25, -.7

    for base in p2_initial_base:
        if base == 'A' or 'T':

            sum_p2_delta_h += AT_h
            sum_p2_delta_s += AT_s

            return sum_p2_delta_h, sum_p2_delta_s

    for base in p2_initial_base:
        if base == 'G' or 'C':

            sum_p2_delta_h += GC_h
            sum_p2_delta_s += GC_s

            return sum_p2_delta_h, sum_p2_delta_s

    for base in p2_terminal_base:
        if base == 'A' or 'T':

            sum_p2_delta_h += AT_h
            sum_p2_delta_s += AT_s

            return sum_p2_delta_h, sum_p2_delta_s

    for base in p2_terminal_base:
        if base == 'G' or 'C':

            sum_p2_delta_h += GC_h
            sum_p2_delta_s += GC_s

            return sum_p2_delta_h, sum_p2_delta_s
            
    return sum_p2_delta_h, sum_p2_delta_s

p1_term_h, p1_term_s = p1_terminal_comp(p1_initial_base,p1_terminal_base, sum_p1_delta_h, sum_p1_delta_s)
p2_term_h, p2_term_s = p2_terminal_comp(p2_initial_base,p2_terminal_base, sum_p2_delta_h, sum_p2_delta_s)

###################################################################################################################################
#Melting Temperature Adjustments
###################################################################################################################################

#Determines the melting temperature of the primers
primer1_gc = get_primer1_gc(primer1)
primer2_gc = get_primer2_gc(primer2)

primer1_melting_temperature = (1000 * p1_term_h) / (p1_term_s + (gas_constant * (math.log(oligo_c)))) - 273.15
primer2_melting_temperature = (1000 * p2_term_h) / (p2_term_s + (gas_constant * (math.log(oligo_c)))) - 273.15

def primer1_salt_correction(Mono, Mg, dNTPs):
    
    Mon = Mono / 2.0
    mg_adj = Mg * 1e-3
    mon = Mon * 1e-3
    dntps = dNTPs * 1e-3 
    ka = 3e4
    mg = (-(ka * dntps - ka * mg_adj + 1.0) + math.sqrt((ka * dntps - ka * mg_adj + 1.0) ** 2 + 4.0 * ka * mg_adj)) / (2.0 * ka)
    R = math.sqrt(mg) / mon 
    b = -9.11e-6 
    c = 6.16e-5
    e = -4.48e-4 
    f = 5.71e-4

    if R < 0.22:

        inverse_primer1_corr = (1 / (primer1_melting_temperature + 273.15)) + ((5.15e-5 * (primer1_gc / 100)) - 2.90e-5) * math.log(mon) + 8.81e-6 * (math.log(mon)) ** 2

        corr = (1 / inverse_primer1_corr) - 273.15

        return round(corr,1)

    elif R < 6.0:

        a = 3.92e-5 * (0.843 - (0.352 * math.sqrt(mon) * math.log(mon)))

        d = 1.42e-5 * ((1.279 - 4.03e-3 * math.log(mon)) - 8.03e-3 * (math.log(mon) ** 2))

        g = 8.31e-5 * ((0.486 - 0.258 * math.log(mon)) + 5.25e-3 * (math.log(mon) ** 3)) 

        inverse_primer1_corr = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        
        corr = (1 / inverse_primer1_corr) - 273.15

        return round(corr,1)
    
    elif R > 6.0:

        a = 3.92e-5

        d = 1.42e-5

        g = 8.31e-5

        inverse_primer1_corr = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)

        corr = (1 / inverse_primer1_corr) - 273.15

        return round(corr,1)

def primer2_salt_correction(Mono, Mg, dNTPs):
    
    Mon = Mono / 2.0
    mg_adj = Mg * 1e-3
    mon = Mon * 1e-3
    dntps = dNTPs * 1e-3 
    ka = 3e4
    mg = (-(ka * dntps - ka * mg_adj + 1.0) + math.sqrt((ka * dntps - ka * mg_adj + 1.0) ** 2 + 4.0 * ka * mg_adj)) / (2.0 * ka)
    R = math.sqrt(mg) / mon 
    b = -9.11e-6 
    c = 6.16e-5
    e = -4.48e-4 
    f = 5.71e-4

    if R < 0.22:
        
        inverse_primer2_corr = (1 / (primer2_melting_temperature + 273.15)) + ((5.15e-5 * (primer2_gc / 100)) - 2.90e-5) * math.log(mon) + 8.81e-6 * (math.log(mon)) ** 2

        corr = (1 / inverse_primer2_corr) - 273.15

        return round(corr,1)

    elif R < 6.0:

        a = 3.92e-5 * (0.843 - (0.352 * math.sqrt(mon) * math.log(mon)))

        d = 1.42e-5 * ((1.279 - 4.03e-3 * math.log(mon)) - 8.03e-3 * (math.log(mon) ** 2))

        g = 8.31e-5 * ((0.486 - 0.258 * math.log(mon)) + 5.25e-3 * (math.log(mon) ** 3)) 

        inverse_primer2_corr = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        
        corr = (1 / inverse_primer2_corr) - 273.15

        return round(corr,1)

    elif R > 6.0:

        a = 3.92e-5

        d = 1.42e-5

        g = 8.31e-5 

        inverse_primer2_corr = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)

        corr = (1 / inverse_primer2_corr) - 273.15

        return round(corr,1)

adj_primer1_melting_temperature = primer1_salt_correction(Mono, Mg, dNTPs)
adj_primer2_melting_temperature = primer2_salt_correction(Mono, Mg, dNTPs)

###################################################################################################################################
#Printing the results
###################################################################################################################################

print('The GC content of primer 1 is ' + str(round(primer1_gc, 1)) + '%.' + 
' The GC content of primer 2 is ' + str(round(primer2_gc, 1)) + '%.')
print('The melting temperatures of your primers are ' + str(round(primer1_melting_temperature, 1)) + 
'C for primer 1, and ' + str(round(primer2_melting_temperature, 1)) + 'C for primer 2.')
print('The salt adjusted melting temperatures of your primers are ' + str(adj_primer1_melting_temperature) + 
'C for primer 1, and ' + str(adj_primer2_melting_temperature) + 'C for primer 2.')
#print(sorted(duo_primer1))
#print(sorted(duo_primer2))