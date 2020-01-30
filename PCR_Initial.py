import Bio
from Bio import Entrez, SeqIO
Entrez.email = "cookienocreams@outlook.com"

#Importing delta G values calculated from the van't Hoff script
from MST_data_converter import AA_delta_g, AT_delta_g, CA_delta_g, CG_delta_g, CT_delta_g, GA_delta_g, GC_delta_g, GG_delta_g, GT_delta_g, TA_delta_g
from MST_data_converter import AA_delta_h, AT_delta_h, CA_delta_h, CG_delta_h, CT_delta_h, GA_delta_h, GC_delta_h, GG_delta_h, GT_delta_h, TA_delta_h
from MST_data_converter import AA_delta_s, AT_delta_s, CA_delta_s, CG_delta_s, CT_delta_s, GA_delta_s, GC_delta_s, GG_delta_s, GT_delta_s, TA_delta_s
import re, math

#These are the ideal variables that will form the standards for the recommendations
ideal_duplex_stabilty = 55
ideal_primer_stability = 60
max_GC_content = 65
min_GC_content = 25
gas_constant_kcal = 1.98720425864083 #cal⋅K−1⋅mol−1
gas_constant_cal = 0.00198720425864083 #kcal⋅K−1⋅mol−1

#Prompts to get the information required to begin the analysis
#Target sequence and primers are input here
#The primer sequences will be analyzed first by calculating their GC content
primer1 = list(input('Input primer 1 sequence here: '))
primer2 = list(input('Input primer 2 sequence here: '))
oligo1_conc = int(input('Input oligo 1 concentration here: '))
oligo2_conc = int(input('Input oligo 2 concentration here: '))

#Calculation of the each primer's GC content
def get_primer1_gc(primer1):

    primer1_lst = list(primer1)

    primer1_GC = str(int(round((sum([1.0 for nucl in primer1_lst if nucl in ['G', 'C']]) / 
    len(primer1_lst)) * 100)))

    return int(primer1_GC)

def get_primer2_gc(primer2):

    primer2_lst = list(primer2)

    primer2_GC = str(int(round((sum([1.0 for nucl in primer2_lst if nucl in ['G', 'C']]) / 
    len(primer2_lst)) * 100)))

    return int(primer2_GC)

lst_primer1 = list(primer1)
lst_primer2 = list(primer2)
primer_length = int(len(primer1))

#Determines the melting temperature of the primers
def short_oligo_tm_1(lst_primer1):

    w = int(sum([1.0 for base in str(lst_primer1) if base in ['A']]))
        
    x = int(sum([1.0 for base in str(lst_primer1) if base in ['T']]))
        
    y = int(sum([1.0 for base in str(lst_primer1) if base in ['G']]))

    z = int(sum([1.0 for base in str(lst_primer1) if base in ['C']]))

    if primer_length < 14:

        tm_value = 2*(w + x) + 4*(y + z)

    elif primer_length > 13:

        tm_value = 64.9 + ((41 * (y + z - 16.4)) / (w + x + y + z))

    return round(tm_value)

def short_oligo_tm_2(lst_primer2):

    w = int(sum([1.0 for base in str(lst_primer2) if base in ['A']]))
        
    x = int(sum([1.0 for base in str(lst_primer2) if base in ['T']]))
        
    y = int(sum([1.0 for base in str(lst_primer2) if base in ['G']]))

    z = int(sum([1.0 for base in str(lst_primer2) if base in ['C']]))

    if primer_length < 14:

        tm_value = 2*(w + x) + 4*(y + z)

    elif primer_length > 13:

        tm_value = 64.9 + ((41 * (y + z - 16.4)) / (w + x + y + z))

    return round(tm_value)

#Prompts to get the information required to begin the analysis
chrID = input('Input gene ID here: ')
start = input('Put start sequence number here: ')
stop = input('Put stop sequnce number here: ')

#This function will be used to aquire the fasta file from the Entrez database
def obtainfasta(chrID, start, stop):

    record = str(SeqIO.read(Entrez.efetch(db='nucleotide', id=chrID, rettype='fasta', 
    strand='1', seq_start=start, seq_stop=stop), 'fasta').seq)

    return record

#Calculates the GC content of a given DNA sequnce input and determines whether that 
#is within the ideal range
def user_input_gc(obtainfasta):

    GC_calculation = (str(int(round((sum([1.0 for nucl in obtainfasta(chrID,start,stop) if 
    nucl in ['G', 'C']]) / len(obtainfasta(chrID,start,stop))) * 100))))

    if str(GC_calculation) > str(max_GC_content) or str(GC_calculation) < str(min_GC_content):

        print('This sequence contains a GC percent of ' +  str(GC_calculation) + '%' + 
        '  and outside the ideal range.')
    else:

        print('The GC content of this sequence is ' + str(GC_calculation) + '%' + 
        ' and within the ideal range.') 

sequence = obtainfasta(chrID, start, stop)
#Splits input sequence into doublets for easy NN matching below
duo_sequence = re.findall('.{1,2}', sequence) + re.findall('.{1,2}', sequence[1:])

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

#Sum of all NN Gibb's free energy values for the input sequence
sum_NN_delta_g = round(AA_NN_sum + AT_NN_sum + CA_NN_sum + CG_NN_sum + CT_NN_sum + GA_NN_sum + 
GC_NN_sum + GT_NN_sum + GG_NN_sum + TA_NN_sum)
avg_delta_G = round((sum_NN_delta_g) / (len(duo_sequence) / 2), 2)

terminal_base = sequence[-1]
initial_base = sequence[0]

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

#Calculation of delta_H_total for each NN using data from van't Hoff script
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

#Sum of all NN enthalpy values for the input sequence
sum_NN_delta_h = round(AA_h_sum + AT_h_sum + CA_h_sum + CG_h_sum + CT_h_sum + GA_h_sum + GC_h_sum + 
GT_h_sum + GG_h_sum + TA_h_sum, 2)
avg_delta_H = round((sum_NN_delta_h) / (len(duo_sequence) / 2), 2)

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

#Sum of all NN entropy values for the input sequence
sum_NN_delta_s = round(AA_s_sum + AT_s_sum + CA_s_sum + CG_s_sum + CT_s_sum + GA_s_sum + GC_s_sum + 
GT_s_sum + GG_s_sum + TA_s_sum, 2)
avg_delta_S = round((sum_NN_delta_s) / (len(duo_sequence) / 2), 2)

#Input sequence melting temperature calculation
oligo_c = (oligo1_conc - (oligo2_conc / 2.0)) * 1e-9
melting_temperature = round((1000 * sum_NN_delta_h) / (sum_NN_delta_s + (gas_constant_kcal * 
(math.log(oligo_c)))) - 273.15)

#Displaying the results of the above calculations
print('The GC content of primer 1 is ' + str(get_primer1_gc(primer1)) + '%.' + 
' The GC content of primer 2 is ' + str(get_primer2_gc(primer2)) + '%.')
print('The melting temperature of your primers is ' + str(short_oligo_tm_1(lst_primer1)) + 
'C for primer 1, and ' + str(short_oligo_tm_2(lst_primer2)) + 'C for primer 2.')
user_input_gc(obtainfasta)
print('The calculated delta G total is ' + str(g_total(sum_NN_delta_g)))
print('The calculated delta H total is ' + str(sum_NN_delta_h))
print('The calculated delta S total is ' + str(sum_NN_delta_s))
print('The calculated meliting temperature for the input sequnce is ' + str(melting_temperature) + 'C.')
print('The calculated average delta G is ' + str(avg_delta_G))
print('The calculated average delta H is ' + str(avg_delta_H))
print('The calculated average delta S is ' + str(avg_delta_S))