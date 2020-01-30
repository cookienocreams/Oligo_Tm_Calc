import Bio
from Bio import Entrez, SeqIO
import re, math
Entrez.email = "cookienocreams@outlook.com"

###################################################################################################################################
#Unified NN values from SantaLucia 1998 paper
###################################################################################################################################

AA_delta_s = -22.2
AA_delta_h = -7.9
AA_delta_g = -1.0

AT_delta_s = -20.4
AT_delta_h = -7.2
AT_delta_g = -0.88

TA_delta_s = -21.3
TA_delta_h = -7.2
TA_delta_g = -0.58

CA_delta_s = -22.7
CA_delta_h = -8.5
CA_delta_g = -1.45

GT_delta_s = -22.4
GT_delta_h = -8.4
GT_delta_g = -1.44

CT_delta_s = -21.0
CT_delta_h = -7.8
CT_delta_g = -1.28

GA_delta_s = -22.2
GA_delta_h = -8.2
GA_delta_g = -1.3

CG_delta_s = -27.2
CG_delta_h = -10.6
CG_delta_g = -2.17

GC_delta_s = -24.4
GC_delta_h = -9.8
GC_delta_g = -2.24

GG_delta_s = -19.9
GG_delta_h = -8.0
GG_delta_g = -1.84

#Lists with all non-overlapping NN pair possiblities
AA = ['AA', 'TT']
AT = ['AT']
CA = ['CA', 'AC']
CG = ['CG']
CT = ['CT', 'TC']
GA = ['GA', 'AG']
GC = ['GC']
GT = ['GT', 'TG']
GG = ['GG', 'CC']
TA = ['TA']

#These are the ideal variables that will form the standards for the recommendations
ideal_duplex_stabilty = 101
ideal_primer_stability = 60
max_GC_content = 65
min_GC_content = 30
gas_constant_cal = 1.98720425864083 #cal⋅K−1⋅mol−1
gas_constant_kcal = 0.00198720425864083 #kcal⋅K−1⋅mol−1

###################################################################################################################################
#Prompts to get the information required to begin the analysis
#Target sequence and primers information are input here
#The primer sequences will be analyzed first by calculating their GC content
#Analysis of primer and input sequence melting temperatures occur last
###################################################################################################################################

primer1 = input('Input primer 1 sequence here: ')
primer2 = input('Input primer 2 sequence here: ')
oligo1_conc = int(input('Input oligo 1 concentration here: '))
oligo2_conc = int(input('Input oligo 2 concentration here: '))
chrID = input('Input gene ID here: ')
start = input('Put start sequence number here: ')
stop = input('Put stop sequnce number here: ')

#Makes each input primer a list so that each nucleotide can be counted individually
lst_primer1 = list(primer1)
lst_primer2 = list(primer2)
primer_length = int(len(primer1))
oligo_c = (oligo1_conc) * 1e-9 #Oligo concentration calculation

#Calculation of the each primer's GC content
def get_primer1_gc(primer1):

    primer1_GC = int(round((sum([1.0 for nucl in lst_primer1 if nucl in ['G', 'C']]) / len(lst_primer1)) * 100))

    return primer1_GC

def get_primer2_gc(primer2):

    primer2_GC = int(round((sum([1.0 for nucl in lst_primer2 if nucl in ['G', 'C']]) / len(lst_primer2)) * 100))

    return primer2_GC

#Splits input primers into doublets for easy NN matching below
duo_primer1 = re.findall('.{1,2}', primer1) + re.findall('.{1,2}', primer1[1:])
duo_primer2 = re.findall('.{1,2}', primer2) + re.findall('.{1,2}', primer2[1:])

###################################################################################################################################
#Determination of NN thermodynamic parameters for each primer
###################################################################################################################################

#Calculation of delta H total for each NN pair in primer 1
AA_h_p1 = AA_delta_h * sum([1.0 for i in duo_primer1 if i in AA])
AT_h_p1 = AT_delta_h * sum([1.0 for i in duo_primer1 if i in AT])
CA_h_p1 = CA_delta_h * sum([1.0 for i in duo_primer1 if i in CA])
CG_h_p1 = CG_delta_h * sum([1.0 for i in duo_primer1 if i in CG])
CT_h_p1 = CT_delta_h * sum([1.0 for i in duo_primer1 if i in CT])
GA_h_p1 = GA_delta_h * sum([1.0 for i in duo_primer1 if i in GA])
GC_h_p1 = GC_delta_h * sum([1.0 for i in duo_primer1 if i in GC])
GT_h_p1 = GT_delta_h * sum([1.0 for i in duo_primer1 if i in GT]) 
GG_h_p1 = GG_delta_h * sum([1.0 for i in duo_primer1 if i in GG])
TA_h_p1 = TA_delta_h * sum([1.0 for i in duo_primer1 if i in TA])

sum_p1_delta_h = round(AA_h_p1 + AT_h_p1 + CA_h_p1 + CG_h_p1 + CT_h_p1 + GA_h_p1 + GC_h_p1 + 
GT_h_p1 + GG_h_p1 + TA_h_p1, 2)

#Calculation of delta S total for each NN pair in primer 1
AA_s_p1 = AA_delta_s * sum([1.0 for i in duo_primer1 if i in AA])
AT_s_p1 = AT_delta_s * sum([1.0 for i in duo_primer1 if i in AT])
CA_s_p1 = CA_delta_s * sum([1.0 for i in duo_primer1 if i in CA])
CG_s_p1 = CG_delta_s * sum([1.0 for i in duo_primer1 if i in CG])
CT_s_p1 = CT_delta_s * sum([1.0 for i in duo_primer1 if i in CT])
GA_s_p1 = GA_delta_s * sum([1.0 for i in duo_primer1 if i in GA])
GC_s_p1 = GC_delta_s * sum([1.0 for i in duo_primer1 if i in GC])
GT_s_p1 = GT_delta_s * sum([1.0 for i in duo_primer1 if i in GT]) 
GG_s_p1 = GG_delta_s * sum([1.0 for i in duo_primer1 if i in GG])
TA_s_p1 = TA_delta_s * sum([1.0 for i in duo_primer1 if i in TA])

sum_p1_delta_s = round(AA_s_p1 + AT_s_p1 + CA_s_p1 + CG_s_p1 + CT_s_p1 + GA_s_p1 + GC_s_p1 + 
GT_s_p1 + GG_s_p1 + TA_s_p1, 2)

#Calculation of delta H total for each NN pair in primer 2
AA_h_p2 = AA_delta_h * sum([1.0 for i in duo_primer2 if i in AA])
AT_h_p2 = AT_delta_h * sum([1.0 for i in duo_primer2 if i in AT])
CA_h_p2 = CA_delta_h * sum([1.0 for i in duo_primer2 if i in CA])
CG_h_p2 = CG_delta_h * sum([1.0 for i in duo_primer2 if i in CG])
CT_h_p2 = CT_delta_h * sum([1.0 for i in duo_primer2 if i in CT])
GA_h_p2 = GA_delta_h * sum([1.0 for i in duo_primer2 if i in GA])
GC_h_p2 = GC_delta_h * sum([1.0 for i in duo_primer2 if i in GC])
GT_h_p2 = GT_delta_h * sum([1.0 for i in duo_primer2 if i in GT]) 
GG_h_p2 = GG_delta_h * sum([1.0 for i in duo_primer2 if i in GG])
TA_h_p2 = TA_delta_h * sum([1.0 for i in duo_primer2 if i in TA])

sum_p2_delta_h = round(AA_h_p2 + AT_h_p2 + CA_h_p2 + CG_h_p2 + CT_h_p2 + GA_h_p2 + GC_h_p2 + 
GT_h_p2 + GG_h_p2 + TA_h_p2, 2)

#Calculation of delta S total for each NN pair in primer 2
AA_s_p2 = AA_delta_s * sum([1.0 for i in duo_primer2 if i in AA])
AT_s_p2 = AT_delta_s * sum([1.0 for i in duo_primer2 if i in AT])
CA_s_p2 = CA_delta_s * sum([1.0 for i in duo_primer2 if i in CA])
CG_s_p2 = CG_delta_s * sum([1.0 for i in duo_primer2 if i in CG])
CT_s_p2 = CT_delta_s * sum([1.0 for i in duo_primer2 if i in CT])
GA_s_p2 = GA_delta_s * sum([1.0 for i in duo_primer2 if i in GA])
GC_s_p2 = GC_delta_s * sum([1.0 for i in duo_primer2 if i in GC])
GT_s_p2 = GT_delta_s * sum([1.0 for i in duo_primer2 if i in GT]) 
GG_s_p2 = GG_delta_s * sum([1.0 for i in duo_primer2 if i in GG])
TA_s_p2 = TA_delta_s * sum([1.0 for i in duo_primer2 if i in TA])

sum_p2_delta_s = round(AA_s_p2 + AT_s_p2 + CA_s_p2 + CG_s_p2 + CT_s_p2 + GA_s_p2 + GC_s_p2 + 
GT_s_p2 + GG_s_p2 + TA_s_p2, 2)

#Determines the melting temperature of the primers
primer1_melting_temperature = round((1000 * sum_p1_delta_h) / (sum_p1_delta_s + (gas_constant_cal * (math.log(oligo_c)))) - 273.15)
primer2_melting_temperature = round((1000 * sum_p2_delta_h) / (sum_p2_delta_s + (gas_constant_cal * (math.log(oligo_c)))) - 273.15)

###################################################################################################################################
#Input Sequence Gathering and Characteriztion
###################################################################################################################################

#This function will be used to aquire the fasta file of the input sequence from the Entrez database
def obtainfasta(chrID, start, stop):

    record = str(SeqIO.read(Entrez.efetch(db='nucleotide', id=chrID, rettype='fasta', 
    strand='1', seq_start=start, seq_stop=stop), 'fasta').seq)

    return record

#Calculates the GC content of a given DNA sequnce input and determines whether that percentage is within our ideal range
GC_calculation = round((sum([1.0 for nucl in obtainfasta(chrID,start,stop) if 
nucl in ['G', 'C']]) / len(obtainfasta(chrID,start,stop))) * 100)

def user_input_gc(obtainfasta):

    if GC_calculation > max_GC_content or GC_calculation < min_GC_content:

        print('This sequence contains a GC percent of ' +  str(GC_calculation) + '%' + ' and outside the ideal range.')

    else:

        print('The GC content of this sequence is ' + str(GC_calculation) + '%' + ' and within the ideal range.') 

#Splits input sequence into doublets for easy NN matching below
sequence = obtainfasta(chrID, start, stop)
duo_sequence = re.findall('.{1,2}', sequence) + re.findall('.{1,2}', sequence[1:])

terminal_base = sequence[-1]
initial_base = sequence[0]

###################################################################################################################################
#Determination of NN thermodynamic parameters for the input sequence
###################################################################################################################################

#Calculation of delta G total for each NN using data from van't Hoff script
AA_NN_sum = AA_delta_g * sum([1.0 for i in duo_sequence if i in AA])
AT_NN_sum = AT_delta_g * sum([1.0 for i in duo_sequence if i in AT])
CA_NN_sum = CA_delta_g * sum([1.0 for i in duo_sequence if i in CA])
CG_NN_sum = CG_delta_g * sum([1.0 for i in duo_sequence if i in CG])
CT_NN_sum = CT_delta_g * sum([1.0 for i in duo_sequence if i in CT])
GA_NN_sum = GA_delta_g * sum([1.0 for i in duo_sequence if i in GA])
GC_NN_sum = GC_delta_g * sum([1.0 for i in duo_sequence if i in GC])
GT_NN_sum = GT_delta_g * sum([1.0 for i in duo_sequence if i in GT]) 
GG_NN_sum = GG_delta_g * sum([1.0 for i in duo_sequence if i in GG])
TA_NN_sum = TA_delta_g * sum([1.0 for i in duo_sequence if i in TA])

#Calculation of delta_H_total for each NN
AA_h_sum = AA_delta_h * sum([1.0 for i in duo_sequence if i in AA])
AT_h_sum = AT_delta_h * sum([1.0 for i in duo_sequence if i in AT])
CA_h_sum = CA_delta_h * sum([1.0 for i in duo_sequence if i in CA])
CG_h_sum = CG_delta_h * sum([1.0 for i in duo_sequence if i in CG])
CT_h_sum = CT_delta_h * sum([1.0 for i in duo_sequence if i in CT])
GA_h_sum = GA_delta_h * sum([1.0 for i in duo_sequence if i in GA])
GC_h_sum = GC_delta_h * sum([1.0 for i in duo_sequence if i in GC])
GT_h_sum = GT_delta_h * sum([1.0 for i in duo_sequence if i in GT])
GG_h_sum = GG_delta_h * sum([1.0 for i in duo_sequence if i in GG])
TA_h_sum = TA_delta_h * sum([1.0 for i in duo_sequence if i in TA])

#Calculation of delta_S_total for each NN using data from van't Hoff script
AA_s_sum = AA_delta_s * sum([1.0 for i in duo_sequence if i in AA])
AT_s_sum = AT_delta_s * sum([1.0 for i in duo_sequence if i in AT])
CA_s_sum = CA_delta_s * sum([1.0 for i in duo_sequence if i in CA])
CG_s_sum = CG_delta_s * sum([1.0 for i in duo_sequence if i in CG])
CT_s_sum = CT_delta_s * sum([1.0 for i in duo_sequence if i in CT])
GA_s_sum = GA_delta_s * sum([1.0 for i in duo_sequence if i in GA])
GC_s_sum = GC_delta_s * sum([1.0 for i in duo_sequence if i in GC])
GT_s_sum = GT_delta_s * sum([1.0 for i in duo_sequence if i in GT]) 
GG_s_sum = GG_delta_s * sum([1.0 for i in duo_sequence if i in GG])
TA_s_sum = TA_delta_s * sum([1.0 for i in duo_sequence if i in TA])

#Sum of all NN enthalpy values for the input sequence
sum_NN_delta_h = round(AA_h_sum + AT_h_sum + CA_h_sum + CG_h_sum + CT_h_sum + GA_h_sum + GC_h_sum + 
GT_h_sum + GG_h_sum + TA_h_sum, 2)
avg_delta_H = round((sum_NN_delta_h) / (len(duo_sequence) / 2), 2)

#Sum of all NN entropy values for the input sequence
sum_NN_delta_s = round(AA_s_sum + AT_s_sum + CA_s_sum + CG_s_sum + CT_s_sum + GA_s_sum + GC_s_sum + 
GT_s_sum + GG_s_sum + TA_s_sum, 2)
avg_delta_S = round((sum_NN_delta_s) / (len(duo_sequence) / 2), 2)

#Sum of all NN Gibb's free energy values for the input sequence
sum_NN_delta_g = round(AA_NN_sum + AT_NN_sum + CA_NN_sum + CG_NN_sum + CT_NN_sum + GA_NN_sum + 
GC_NN_sum + GT_NN_sum + GG_NN_sum + TA_NN_sum, 1)
avg_delta_G = round((sum_NN_delta_g) / (len(duo_sequence) / 2), 2)

#Compensation for differences in free energies bewteen the inital bases and terminal bases
def g_total(sum_NN_delta_g):

    if terminal_base == 'G' or terminal_base == 'C' and initial_base == 'G' or initial_base == 'C':

        delta_g_total = sum_NN_delta_g + (.98)*2

        return delta_g_total

    elif initial_base == 'G' or initial_base == 'C' and terminal_base == 'A' or terminal_base == 'T':

        delta_g_total = sum_NN_delta_g + .98 + 1.03

        return delta_g_total

    elif terminal_base == 'A' or terminal_base == 'T' and initial_base == 'A' or initial_base == 'T':

        delta_g_total = sum_NN_delta_g + (1.03)*2

        return delta_g_total

    elif initial_base == 'A' or initial_base == 'T' and terminal_base == 'G' or terminal_base == 'C':

        delta_g_total = sum_NN_delta_g + 1.03 + .98

        return delta_g_total

###################################################################################################################################
#Melting Temperature Adjustments
###################################################################################################################################

#Input sequence melting temperature calculation
oligo_c = (oligo1_conc - (oligo2_conc / 2.0)) * 1e-9
melting_temperature = round((1000 * sum_NN_delta_h) / (sum_NN_delta_s + (gas_constant_cal * (math.log(oligo_c)))) - 273.15)

Na = 0
K = 0
Tris = 25
Mg = 1.5
dNTPs = .2

def salt_correction(Na, K, Tris, Mg, dNTPs):
    
    Mon = Na + K + Tris / 2.0
    mg_adj = Mg * 1e-3
    mon = Mon * 1e-3
    dntps = dNTPs * 1e-3 
    ka = 3e4
    mg = (-(ka * dntps - ka * mg_adj + 1.0) + math.sqrt((ka * dntps - ka * mg_adj + 1.0) ** 2 + 4.0 * ka * mg_adj)) / (2.0 * ka)
    R = math.sqrt(mg) / mon 
    b, c, e, f = -9.11e-6, 6.26e-5, -4.82e-4, 5.25e-4

    if R < 0.22:

        corr = (4.29 * GC_calculation / 100 - 3.95) * 1e-5 * math.log(mon) + 9.40e-6 * math.log(mon) ** 2 

        return round(corr, 1)

    elif R < 6.0:

        a = 3.92e-5 * (0.843 - (0.352 * math.sqrt(mon) * math.log(mon)))

        d = 1.42e-5 * ((1.279 - 4.03e-3 * math.log(mon)) - 8.03e-3 * (math.log(mon) ** 2))

        g = 8.31e-5 * ((0.486 - 0.258 * math.log(mon)) + 5.25e-3 * (math.log(mon) ** 3)) 

        inverse_mg_corr = 1 / melting_temperature + (a + b * math.log(mg)) + int(GC_calculation) / 100 * (c + d * math.log(mg)) + \
        (1 / (2.0 * (len(sequence) - 1))) * (e + f * math.log(mg) + g * (math.log(mg) ** 2))
        
        corr = 1 / inverse_mg_corr

        return round(corr,1)

salty_melting_temperature = salt_correction(Na, K, Tris, Mg, dNTPs)

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

adjusted_melting_temperature = DMSO_correction(salty_melting_temperature, GC_calculation)

###################################################################################################################################
#Printing the results
###################################################################################################################################

print('The GC content of primer 1 is ' + str(get_primer1_gc(primer1)) + '%.' + 
' The GC content of primer 2 is ' + str(get_primer1_gc(primer2)) + '%.')
print('The melting temperature of your primers is ' + str(primer1_melting_temperature) + 
'C for primer 1, and ' + str(primer2_melting_temperature) + 'C for primer 2.')
user_input_gc(obtainfasta)
print('The calculated delta G total is ' + str(g_total(sum_NN_delta_g)))
print('The calculated delta H total is ' + str(sum_NN_delta_h))
print('The calculated delta S total is ' + str(sum_NN_delta_s))
print('The calculated meliting temperature for the input sequnce is ' + str(melting_temperature) + 'C.')
print('The salt adjusted meliting temperature for the input sequnce is ' + str(salty_melting_temperature) + 'C.')
print('The final adjusted meliting temperature for the input sequnce is ' + str(adjusted_melting_temperature) + 'C.')
#print('The calculated average delta G is ' + str(avg_delta_G))
#print('The calculated average delta H is ' + str(avg_delta_H))
#print('The calculated average delta S is ' + str(avg_delta_S))