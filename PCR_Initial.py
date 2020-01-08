import Bio
from Bio import Entrez, SeqIO
Entrez.email = "cookienocreams@outlook.com"

#These are the ideal variables that will form the standards for the recommendations
ideal_duplex_stabilty = 55
ideal_polymerase_stability = 70
ideal_primer_stability = 60
max_GC_content = 65
min_GC_content = 25
propane_1_2_diol_factor = .8
DMSO_factor = .6
TMAC_factor = 1.2 or .7
betaine_factor = .7
deaza_factor = .85
BSA_factor = 1.05
ammonium_sulfate_factor = 1.2
detergents_factor = 1.1

#Prompts to get 
chrID = input('Input gene ID here: ')
start = input('Put start sequence number here: ')
stop = input('Put stop sequnce number here: ')

#This can be used to aquire a given fasta file from the Entrez database
def obtainfasta(chrID, start, stop):
    record = str(SeqIO.read(Entrez.efetch(db='nucleotide', id=chrID, rettype='fasta', strand='1', seq_start=start, seq_stop=stop), 'fasta').seq)
    return record

#Calculates the GC content of a given DNA sequnce input and determines whether that is within the ideal range
def user_input_gc(obtainfasta):
    GC_calculation = (str(int(round((sum([1.0 for nucl in obtainfasta(chrID,start,stop) if nucl in ['G', 'C']]) / len(obtainfasta(chrID,start,stop))) * 100))))
    #GC_calculation = sum(getFas_list.count(x) for x in ['G', 'C', 'g', 'c']) 
    if str(GC_calculation) > str(max_GC_content) or str(GC_calculation) < str(min_GC_content):
        print('This sequence contains a GC percent of ' +  str(GC_calculation) + '%' + '  and outside the normal range')
    else:
        print('The GC content of this sequence is ' + str(GC_calculation) + '%' + ' and within the normal range') 
    return GC_calculation

#This is where input DNA stabilty will be compared with the "ideal"
def duplex_stability(measurement1):
      if measurement1 >= ideal_duplex_stabilty:
          return 'This DNA is stable'
      elif measurement1 < ideal_duplex_stabilty:
          return 'This DNA is unstable'
#This is where input polymerase will be compared with the "ideal"
def polymerase_stability(measurement2):
      if measurement2 >= ideal_polymerase_stability:
          return 'This polymerase is stable'
      elif measurement2 < ideal_polymerase_stability:
          return 'This polymerase is unstable'
#This is where input primers will be compared with the "ideal"
def polymerase_stability(measurement3):
      if measurement3 >= ideal_primer_stability:
          return 'These primers are stable'
      elif measurement3 < ideal_primer_stability:
          return 'These primers are unstable'

print(str(user_input_gc(obtainfasta)))