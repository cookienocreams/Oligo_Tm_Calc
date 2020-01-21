import Bio
from Bio import Entrez, SeqIO
Entrez.email = "cookienocreams@outlook.com"

#These are the ideal variables that will form the standards for the recommendations
ideal_duplex_stabilty = 55
ideal_polymerase_stability = 70
ideal_primer_stability = 60
max_GC_content = 65
min_GC_content = 25

#Prompts to get the information required to begin the analysis
chrID = input('Input gene ID here: ')
start = input('Put start sequence number here: ')
stop = input('Put stop sequnce number here: ')

#This function will be used to aquire fasta file from the Entrez database
def obtainfasta(chrID, start, stop):
    record = str(SeqIO.read(Entrez.efetch(db='nucleotide', id=chrID, rettype='fasta', strand='1', seq_start=start, seq_stop=stop), 'fasta').seq)
    return record

#Calculates the GC content of a given DNA sequnce input and determines whether that is within the ideal range
def user_input_gc(obtainfasta):
    GC_calculation = (str(int(round((sum([1.0 for nucl in obtainfasta(chrID,start,stop) if nucl in ['G', 'C']]) / len(obtainfasta(chrID,start,stop))) * 100))))
    if str(GC_calculation) > str(max_GC_content) or str(GC_calculation) < str(min_GC_content):
        print('This sequence contains a GC percent of ' +  str(GC_calculation) + '%' + '  and outside the normal range.')
    else:
        print('The GC content of this sequence is ' + str(GC_calculation) + '%' + ' and within the normal range.') 

#Importing delta G values calculated from the van't Hoff script
from std_NN_parameters import AA_delta_g, AT_delta_g, CA_delta_g, CG_delta_g, CT_delta_g, GA_delta_g, GC_delta_g, GG_delta_g, GT_delta_g, TA_delta_g
from std_NN_parameters import AA_delta_h, AT_delta_h, CA_delta_h, CG_delta_h, CT_delta_h, GA_delta_h, GC_delta_h, GG_delta_h, GT_delta_h, TA_delta_h
from std_NN_parameters import AA_delta_s, AT_delta_s, CA_delta_s, CG_delta_s, CT_delta_s, GA_delta_s, GC_delta_s, GG_delta_s, GT_delta_s, TA_delta_s
import re, math

gas_constant = 1.98720425864083
sequence = obtainfasta(chrID, start, stop)
#Splits sequence into doublets for easy NN matching below
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

#Calculation of delta_G_total for each NN using data from van't Hoff script
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
sum_NN_delta_g = round(AA_NN_sum + AT_NN_sum + CA_NN_sum + CG_NN_sum + CT_NN_sum + GA_NN_sum + GC_NN_sum + GT_NN_sum + GG_NN_sum + TA_NN_sum, 2)
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
sum_NN_delta_h = round(AA_h_sum + AT_h_sum + CA_h_sum + CG_h_sum + CT_h_sum + GA_h_sum + GC_h_sum + GT_h_sum + GG_h_sum + TA_h_sum, 2)

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
sum_NN_delta_s = round(AA_s_sum + AT_s_sum + CA_s_sum + CG_s_sum + CT_s_sum + GA_s_sum + GC_s_sum + GT_s_sum + GG_s_sum + TA_s_sum, 2)

user_input_gc(obtainfasta)
print('The total calculated delta G total is ' + str(g_total(sum_NN_delta_g)))
print('The total calculated delta H total is ' + str(sum_NN_delta_h))
print('The total calculated delta S total is ' + str(sum_NN_delta_s))