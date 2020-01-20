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
duo_sequence = re.findall('.{1,2}', sequence) + re.findall('.{1,2}', sequence[1:])

#Calculation of delta_G_total for each NN using data from van't Hoff script
AA_NN_sum = AA_delta_g * sum([1.0 for i in duo_sequence if i in ['AA']])
TT_NN_sum = AA_delta_g * sum([1.0 for i in duo_sequence if i in ['TT']])
AT_NN_sum = AT_delta_g * sum([1.0 for i in duo_sequence if i in ['AT']])
CA_NN_sum = CA_delta_g * sum([1.0 for i in duo_sequence if i in ['CA']])
TG_NN_sum = CA_delta_g * sum([1.0 for i in duo_sequence if i in ['TG']])
AC_NN_sum = CA_delta_g * sum([1.0 for i in duo_sequence if i in ['AC']])
CG_NN_sum = CG_delta_g * sum([1.0 for i in duo_sequence if i in ['CG']])
CT_NN_sum = CT_delta_g * sum([1.0 for i in duo_sequence if i in ['CT']])
GA_NN_sum = GA_delta_g * sum([1.0 for i in duo_sequence if i in ['GA']])
GC_NN_sum = GC_delta_g * sum([1.0 for i in duo_sequence if i in ['GC']])
GT_NN_sum = GT_delta_g * sum([1.0 for i in duo_sequence if i in ['GT']]) 
GG_NN_sum = GG_delta_g * sum([1.0 for i in duo_sequence if i in ['GG']])
CC_NN_sum = GG_delta_g * sum([1.0 for i in duo_sequence if i in ['CC']])
TA_NN_sum = TA_delta_g * sum([1.0 for i in duo_sequence if i in ['TA']])

sum_NN_delta_g = AA_NN_sum + AT_NN_sum + CA_NN_sum + CG_NN_sum + CT_NN_sum + GA_NN_sum + GC_NN_sum + GT_NN_sum + GG_NN_sum + TA_NN_sum + TT_NN_sum + CC_NN_sum + AC_NN_sum + TG_NN_sum
terminal_base = sequence[-1]
initial_base = sequence[0]

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

user_input_gc(obtainfasta)
print(g_total(sum_NN_delta_g))