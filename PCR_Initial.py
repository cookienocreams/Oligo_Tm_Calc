import Bio
from Bio import Entrez, SeqIO
import re, math
Entrez.email = "cookienocreams@outlook.com"

###################################################################################################################################
#Modified NN Parameters based on SantaLucia's 1998 paper
###################################################################################################################################

AA_delta_s = -22.15
AA_delta_h = -7.845
AA_delta_g = AA_delta_h - (273.15*(AA_delta_s/1000))

AC_delta_s = -20.125
AC_delta_h = -7.97
AC_delta_g = AC_delta_h - (273.15*(AC_delta_s/1000))

AG_delta_s = -21.125
AG_delta_h = -8.335
AG_delta_g = AG_delta_h - (273.15*(AG_delta_s/1000))

AT_delta_s = -20.315
AT_delta_h = -7.12
AT_delta_g = AT_delta_h - (273.15*(AT_delta_s/1000))

CA_delta_s = -22.27
CA_delta_h = -8.63
CA_delta_g = CA_delta_h - (273.15*(CA_delta_s/1000))

CC_delta_s = -20.65
CC_delta_h = -7.825
CC_delta_g = CC_delta_h - (273.15*(CC_delta_s/1000))

CG_delta_s = -24.95
CG_delta_h = -10.02
CG_delta_g = CG_delta_h - (273.15*(CG_delta_s/1000))

GA_delta_s = -22.47
GA_delta_h = -8.43
GA_delta_g = GA_delta_h - (273.15*(GA_delta_s/1000))

GC_delta_s = -24.7
GC_delta_h = -9.65
GC_delta_g = GC_delta_h - (273.15*(GC_delta_s/1000))

TA_delta_s = -20.365
TA_delta_h = -7.2
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

#These are the ideal variables that will form the standards for the recommendations
ideal_duplex_stabilty = 101 #In Celcius
max_GC_content = 65
min_GC_content = 30
gas_constant = 1.98720425864083 #cal⋅K−1⋅mol−1

###################################################################################################################################
#Prompts to get the information required to begin the analysis
#Target DNA sequence and primer information are input here
#The primer sequences will be analyzed first by calculating their GC content
#Analysis of primer and input sequence melting temperatures occur last
###################################################################################################################################

primer1 = input('Input primer 1 sequence here: ')
primer2 = input('Input primer 2 sequence here: ')
oligo1_conc = float(input('Input total single strand oligo concentration in uM here: '))
oligo2_conc = float(input('Input template concentration in uM here: '))
chrID = input('Input gene ID here: ')
start = input('Put start sequence number here: ')
stop = input('Put stop sequence number here: ')

#Captures the first and last base in each primer so that the energies associated with opening the helix can be calculated
p1_initial_base = primer1[0]
p1_terminal_base = primer1[-1]

p2_initial_base = primer2[0]
p2_terminal_base = primer2[-1]

#Splits input primers into doublets for easy NN matching below
duo_primer1 = re.findall('.{1,2}', primer1) + re.findall('.{1,2}', primer1[1:])
duo_primer2 = re.findall('.{1,2}', primer2) + re.findall('.{1,2}', primer2[1:])

primer1_length = len(primer1)
primer2_length = len(primer2)

#Determines which primer calculation should be used based on the ratio of primer to template
def oligo_calculation(oligo1_conc, oligo2_conc):

    if oligo1_conc > 6 * oligo2_conc:

        oligo_c = (oligo1_conc) * 1e-9 #Converts to mol/L

        return oligo_c

    elif oligo1_conc < 6 * oligo2_conc and oligo1_conc > oligo2_conc:

        oligo_c = (oligo1_conc - (oligo2_conc / 2.0)) * 1e-9

        return oligo_c

    elif oligo1_conc == oligo2_conc:

        oligo_c = (oligo1_conc / 2.0) * 1e-9

        return oligo_c

oligo_c = oligo_calculation(oligo1_conc, oligo2_conc)
template_c = (oligo2_conc / 2.0) * 1e-9

#Calculation of the each primer's GC content
def get_primer_gc(primer1, primer2):

    primer1_GC = (sum(1.0 for base in primer1 if base in ['G', 'C']) / primer1_length) * 100
    primer2_GC = (sum(1.0 for base in primer2 if base in ['G', 'C']) / primer2_length) * 100

    return primer1_GC, primer2_GC

primer1_gc, primer2_gc = get_primer_gc(primer1, primer2)

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

#Calculates the entropic and enthalpic compensation for the ends of each primer
def p1_terminal_comp(p1_initial_base, p1_terminal_base, sum_p1_delta_h, sum_p1_delta_s):

    AT_h, AT_s, GC_h, GC_s = 2.4, 5.7, .45, -1.7 #h in kcal/mol, s in eu
    dh, ds = 0, 0

    dh_init = (dh + AT_h if p1_initial_base in ['A', 'T'] else dh + GC_h)
    ds_init = (ds + AT_s if p1_initial_base in ['A', 'T'] else ds + GC_s)

    dh_term = (dh + AT_h if p1_terminal_base in ['A', 'T'] else dh + GC_h)
    ds_term = (ds + AT_s if p1_terminal_base in ['A', 'T'] else ds + GC_s)

    dh_total = dh_init + dh_term + sum_p1_delta_h
    ds_total = ds_init + ds_term + sum_p1_delta_s

    return dh_total, ds_total

def p2_terminal_comp(p2_initial_base, p2_terminal_base, sum_p2_delta_h, sum_p2_delta_s):
    
    AT_h, AT_s, GC_h, GC_s = 2.4, 5.7, .45, -1.7 #h in kcal/mol, s in eu
    dh, ds = 0, 0

    dh_init = (dh + AT_h if p2_initial_base in ['A', 'T'] else dh + GC_h)
    ds_init = (ds + AT_s if p2_initial_base in ['A', 'T'] else ds + GC_s)

    dh_term = (dh + AT_h if p2_terminal_base in ['A', 'T'] else dh + GC_h)
    ds_term = (ds + AT_s if p2_terminal_base in ['A', 'T'] else ds + GC_s)

    dh_total = dh_init + dh_term + sum_p2_delta_h
    ds_total = ds_init + ds_term + sum_p2_delta_s

    return dh_total, ds_total

p1_total_h, p1_total_s = p1_terminal_comp(p1_initial_base, p1_terminal_base, sum_p1_delta_h, sum_p1_delta_s)
p2_total_h, p2_total_s = p2_terminal_comp(p2_initial_base, p2_terminal_base, sum_p2_delta_h, sum_p2_delta_s)

###################################################################################################################################
#Input Sequence Gathering and Characteriztion
###################################################################################################################################

#This function will be used to aquire the fasta file of the input sequence from the Entrez database
def obtainfasta(chrID, start, stop):

    sequence = str(SeqIO.read(Entrez.efetch(db='nucleotide', id=chrID, rettype='fasta', 
    strand='1', seq_start=start, seq_stop=stop), 'fasta').seq)

    return sequence

#Splits input sequence into doublets for easy NN matching below
sequence = obtainfasta(chrID, start, stop)
duo_sequence = re.findall('.{1,2}', sequence) + re.findall('.{1,2}', sequence[1:])

initial_base = sequence[0]
terminal_base = sequence[-1]

#Calculates the GC content of the input DNA sequence and determines whether that percentage is within the ideal range
GC_calculation = (sum(1.0 for base in sequence if base in ['G', 'C']) / len(sequence)) * 100

def input_sequence_gc(GC_calculation):

    if GC_calculation > max_GC_content or GC_calculation < min_GC_content:

        print('This sequence contains a GC percent of ' +  str(round(GC_calculation, 1)) + '%' + ' and outside the ideal range.')

    else:

        print('The GC content of this sequence is ' + str(round(GC_calculation, 1)) + '%' + ' and within the ideal range.') 

###################################################################################################################################
#Tabulation of NN thermodynamic parameters for the input sequence
###################################################################################################################################

#Calculation of delta G total for each NN in the input sequence
AA_NN_sum = AA_delta_g * sum(1.0 for NN in duo_sequence if NN in AA)
AC_NN_sum = AC_delta_g * sum(1.0 for NN in duo_sequence if NN in AC)
AG_NN_sum = AG_delta_g * sum(1.0 for NN in duo_sequence if NN in AG)
AT_NN_sum = AT_delta_g * sum(1.0 for NN in duo_sequence if NN in AT)
CA_NN_sum = CA_delta_g * sum(1.0 for NN in duo_sequence if NN in CA)
CC_NN_sum = CC_delta_g * sum(1.0 for NN in duo_sequence if NN in CC)
CG_NN_sum = CG_delta_g * sum(1.0 for NN in duo_sequence if NN in CG)
GA_NN_sum = GA_delta_g * sum(1.0 for NN in duo_sequence if NN in GA) 
GC_NN_sum = GC_delta_g * sum(1.0 for NN in duo_sequence if NN in GC)
TA_NN_sum = TA_delta_g * sum(1.0 for NN in duo_sequence if NN in TA)

#Calculation of delta H total for each NN in the input sequence
AA_h_sum = AA_delta_h * sum(1.0 for NN in duo_sequence if NN in AA)
AC_h_sum = AC_delta_h * sum(1.0 for NN in duo_sequence if NN in AC)
AG_h_sum = AG_delta_h * sum(1.0 for NN in duo_sequence if NN in AG)
AT_h_sum = AT_delta_h * sum(1.0 for NN in duo_sequence if NN in AT)
CA_h_sum = CA_delta_h * sum(1.0 for NN in duo_sequence if NN in CA)
CC_h_sum = CC_delta_h * sum(1.0 for NN in duo_sequence if NN in CC)
CG_h_sum = CG_delta_h * sum(1.0 for NN in duo_sequence if NN in CG)
GA_h_sum = GA_delta_h * sum(1.0 for NN in duo_sequence if NN in GA) 
GC_h_sum = GC_delta_h * sum(1.0 for NN in duo_sequence if NN in GC)
TA_h_sum = TA_delta_h * sum(1.0 for NN in duo_sequence if NN in TA)

#Calculation of delta S total for each NN in the input sequence
AA_s_sum = AA_delta_s * sum(1.0 for NN in duo_sequence if NN in AA)
AC_s_sum = AC_delta_s * sum(1.0 for NN in duo_sequence if NN in AC)
AG_s_sum = AG_delta_s * sum(1.0 for NN in duo_sequence if NN in AG)
AT_s_sum = AT_delta_s * sum(1.0 for NN in duo_sequence if NN in AT)
CA_s_sum = CA_delta_s * sum(1.0 for NN in duo_sequence if NN in CA)
CC_s_sum = CC_delta_s * sum(1.0 for NN in duo_sequence if NN in CC)
CG_s_sum = CG_delta_s * sum(1.0 for NN in duo_sequence if NN in CG)
GA_s_sum = GA_delta_s * sum(1.0 for NN in duo_sequence if NN in GA) 
GC_s_sum = GC_delta_s * sum(1.0 for NN in duo_sequence if NN in GC)
TA_s_sum = TA_delta_s * sum(1.0 for NN in duo_sequence if NN in TA)

#Sum of all NN enthalpy values for the input sequence
sum_NN_delta_h = AA_h_sum + AC_h_sum + AG_h_sum + AT_h_sum + CA_h_sum + CC_h_sum + CG_h_sum + GA_h_sum + GC_h_sum + TA_h_sum
avg_delta_H = round((sum_NN_delta_h) / (len(sequence) / 2), 2)

#Sum of all NN entropy values for the input sequence
sum_NN_delta_s = AA_s_sum + AC_s_sum + AG_s_sum + AT_s_sum + CA_s_sum + CC_s_sum + CG_s_sum + GA_s_sum + GC_s_sum + TA_s_sum
avg_delta_S = round((sum_NN_delta_s) / (len(sequence) / 2), 2)

#Sum of all NN Gibb's free energy values for the input sequence
sum_NN_delta_g = AA_NN_sum + AC_NN_sum + AG_NN_sum + AT_NN_sum + CA_NN_sum + CC_NN_sum + CG_NN_sum + GA_NN_sum + GC_NN_sum + TA_NN_sum
avg_delta_G = round((sum_NN_delta_g) / (len(sequence) / 2), 2)

#Calculates the entropic and enthalpic compensation for the ends of the input sequence
def seq_terminal_comp(initial_base, terminal_base, sum_NN_delta_h, sum_NN_delta_s):

    AT_h, AT_s, GC_h, GC_s = 2.4, 5.7, .45, -1.7 #h in kcal/mol, s in eu
    dh, ds = 0, 0

    dh_init = (dh + AT_h if initial_base in ['A', 'T'] else dh + GC_h)
    ds_init = (ds + AT_s if initial_base in ['A', 'T'] else ds + GC_s)

    dh_term = (dh + AT_h if terminal_base in ['A', 'T'] else dh + GC_h)
    ds_term = (ds + AT_s if terminal_base in ['A', 'T'] else ds + GC_s)

    dh_total = dh_init + dh_term + sum_NN_delta_h
    ds_total = ds_init + ds_term + sum_NN_delta_s

    return dh_total, ds_total

total_h, total_s = seq_terminal_comp(initial_base, terminal_base, sum_NN_delta_h, sum_NN_delta_s)

###################################################################################################################################
#Melting Temperature Adjustments
###################################################################################################################################

#Input sequence melting temperature calculation
melting_temperature = (1000 * total_h) / (total_s + (gas_constant * (math.log(template_c)))) - 273.15

#Determines the melting temperature of the primers
primer1_melting_temperature = (1000 * p1_total_h) / (p1_total_s + (gas_constant * (math.log(oligo_c)))) - 273.15
primer2_melting_temperature = (1000 * p2_total_h) / (p2_total_s + (gas_constant * (math.log(oligo_c)))) - 273.15

#PCR buffer conditions
Na = 0
K = 50
Tris = 20
Mg = 1.5
dNTPs = .2

def salt_correction(Na, K, Tris, Mg, dNTPs):
    
    Mon = (Na + K + Tris) / 2.0
    mg_adj = Mg * 1e-3 #Converts to mol/L
    mon = Mon * 1e-3
    dntps = dNTPs * 1e-3 
    ka = 3e4
    mg = (-(ka * dntps - ka * mg_adj + 1.0) + math.sqrt((ka * dntps - ka * mg_adj + 1.0) ** 2 + 4.0 * ka * mg_adj)) / (2.0 * ka)
    R = math.sqrt(mg) / mon 
    b, c, e, f = -9.11e-6, 6.06e-5, -4.82e-4, 5.65e-4 #in K-1

    if R < 0.22:
        
        inverse_mg_corr = (1 / (melting_temperature + 273.15)) + ((4.59e-5 * (GC_calculation / 100)) - 2.90e-5) * math.log(mon) + 8.81e-6 * (math.log(mon)) ** 2

        corr = (1 / inverse_mg_corr) - 273.15

        return corr

    elif R < 6.0:

        a = 3.92e-5 * (0.843 - (0.352 * math.sqrt(mon) * math.log(mon)))
        d = 1.42e-5 * ((1.279 - 4.03e-3 * math.log(mon)) - 8.03e-3 * (math.log(mon) ** 2))
        g = 8.31e-5 * ((0.486 - 0.258 * math.log(mon)) + 5.25e-3 * (math.log(mon) ** 3)) 

        inverse_mg_corr = (1 / (melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((GC_calculation / 100) * (c + d * math.log(mg))) + (1 / (2.0 * ((len(sequence) - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2))

        corr = (1 / inverse_mg_corr) - 273.15

        return corr

    elif R > 6.0:

        a = 3.92e-5
        d = 1.42e-5
        g = 8.31e-5 

        inverse_mg_corr = (1 / (melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((GC_calculation / 100) * (c + d * math.log(mg))) + (1 / (2.0 * ((len(sequence) - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2))

        corr = (1 / inverse_mg_corr) - 273.15

        return corr

def primer1_salt_correction(Na, K, Tris, Mg, dNTPs):
    
    Mon = (Na + K + Tris) / 2.0
    mg_adj = Mg * 1e-3 #Converts to mol/L
    mon = Mon * 1e-3
    dntps = dNTPs * 1e-3 
    ka = 3e4
    mg = (-(ka * dntps - ka * mg_adj + 1.0) + math.sqrt((ka * dntps - ka * mg_adj + 1.0) ** 2 + 4.0 * ka * mg_adj)) / (2.0 * ka)
    R = math.sqrt(mg) / mon 
    b, c, e, f = -9.11e-6, 6.06e-5, -4.82e-4, 5.65e-4 #in K-1

    if R < 0.22:

        inverse_primer1_corr = (1 / (primer1_melting_temperature + 273.15)) + ((4.59e-5 * (primer1_gc / 100)) - 2.90e-5) * math.log(mon) + 8.81e-6 * (math.log(mon)) ** 2

        corr = (1 / inverse_primer1_corr) - 273.15

        return corr

    elif R < 6.0:

        a = 3.92e-5 * (0.843 - (0.352 * math.sqrt(mon) * math.log(mon)))
        d = 1.42e-5 * ((1.279 - 4.03e-3 * math.log(mon)) - 8.03e-3 * (math.log(mon) ** 2))
        g = 8.31e-5 * ((0.486 - 0.258 * math.log(mon)) + 5.25e-3 * (math.log(mon) ** 3)) 

        inverse_primer1_corr = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        
        corr = (1 / inverse_primer1_corr) - 273.15

        return corr
    
    elif R > 6.0:

        a = 3.92e-5
        d = 1.42e-5
        g = 8.31e-5

        inverse_primer1_corr = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)

        corr = (1 / inverse_primer1_corr) - 273.15

        return corr

def primer2_salt_correction(Na, K, Tris, Mg, dNTPs):
    
    Mon = (Na + K + Tris) / 2.0
    mg_adj = Mg * 1e-3 #Converts to mol/L
    mon = Mon * 1e-3
    dntps = dNTPs * 1e-3 
    ka = 3e4
    mg = (-(ka * dntps - ka * mg_adj + 1.0) + math.sqrt((ka * dntps - ka * mg_adj + 1.0) ** 2 + 4.0 * ka * mg_adj)) / (2.0 * ka)
    R = math.sqrt(mg) / mon 
    b, c, e, f = -9.11e-6, 6.06e-5, -4.82e-4, 5.65e-4 #in K-1

    if R < 0.22:
        
        inverse_primer2_corr = (1 / (primer2_melting_temperature + 273.15)) + ((4.59e-5 * (primer2_gc / 100)) - 2.90e-5) * math.log(mon) + 8.81e-6 * (math.log(mon)) ** 2

        corr = (1 / inverse_primer2_corr) - 273.15

        return corr

    elif R < 6.0:

        a = 3.92e-5 * (0.843 - (0.352 * math.sqrt(mon) * math.log(mon)))
        d = 1.42e-5 * ((1.279 - 4.03e-3 * math.log(mon)) - 8.03e-3 * (math.log(mon) ** 2))
        g = 8.31e-5 * ((0.486 - 0.258 * math.log(mon)) + 5.25e-3 * (math.log(mon) ** 3)) 

        inverse_primer2_corr = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        
        corr = (1 / inverse_primer2_corr) - 273.15

        return corr

    elif R > 6.0:

        a = 3.92e-5
        d = 1.42e-5
        g = 8.31e-5 

        inverse_primer2_corr = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)

        corr = (1 / inverse_primer2_corr) - 273.15

        return corr

salty_melting_temperature = salt_correction(Na, K, Tris, Mg, dNTPs)
adj_primer1_melting_temperature = primer1_salt_correction(Na, K, Tris, Mg, dNTPs)
adj_primer2_melting_temperature = primer2_salt_correction(Na, K, Tris, Mg, dNTPs)

def DMSO_correction(salty_melting_temperature, GC_calculation):

    DMSO = 0
    DMSO_factor = .65
    DMSO_melting_temperature = salty_melting_temperature - DMSO_factor * DMSO

    while DMSO_melting_temperature > ideal_duplex_stabilty:

        DMSO += 1
        DMSO_melting_temperature = salty_melting_temperature - DMSO_factor * DMSO

        if DMSO_melting_temperature <= ideal_duplex_stabilty:

            print('The recommended DMSO concentration that should be used is ' + str(DMSO) + '%.')

            break

    return round(DMSO_melting_temperature, 1)

final_melting_temperature = DMSO_correction(salty_melting_temperature, GC_calculation)

###################################################################################################################################
#Printing the results
###################################################################################################################################

print('The GC content of primer 1 is ' + str(round(primer1_gc,1)) + '%.' + 
' The GC content of primer 2 is ' + str(round(primer2_gc,1)) + '%.')
print('The melting temperatures of your primers are ' + str(round(primer1_melting_temperature, 1)) + 
'C for primer 1, and ' + str(round(primer2_melting_temperature, 1)) + 'C for primer 2.')
print('The salt adjusted melting temperatures of your primers are ' + str(round(adj_primer1_melting_temperature, 1)) + 
'C for primer 1, and ' + str(round(adj_primer2_melting_temperature, 1)) + 'C for primer 2.')
input_sequence_gc(GC_calculation)
#print('The calculated delta G total is ' + str(g_total(sum_NN_delta_g)))
#print('The calculated delta H total is ' + str(sum_NN_delta_h))
#print('The calculated delta S total is ' + str(sum_NN_delta_s))
print('The calculated meliting temperature for the input sequence is ' + str(round(melting_temperature, 1)) + 'C.')
print('The salt adjusted meliting temperature for the input sequence is ' + str(round(salty_melting_temperature, 1)) + 'C.')
print('The final adjusted meliting temperature for the input sequence is ' + str(final_melting_temperature) + 'C.')
#print('The calculated average delta G is ' + str(avg_delta_G))
#print('The calculated average delta H is ' + str(avg_delta_H))
#print('The calculated average delta S is ' + str(avg_delta_S))