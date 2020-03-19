import Bio
from Bio import Entrez, SeqIO
import re, math
Entrez.email = "cookienocreams@outlook.com"

#These are the ideal variables that will form the standards for the recommendations
ideal_duplex_stabilty = 90 #In Celcius
gas_constant = 8.31445984848484848484 #J K-1-mol-1 

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

primer1_length = len(primer1)
primer2_length = len(primer2)

#Determines which primer calculation should be used based on the ratio of primer to template
def oligo_calculation(oligo1_conc, oligo2_conc):

    if oligo1_conc > 6 * oligo2_conc:

        oligo_c = (oligo1_conc) * 1e-6 #Converts to mol/L

        return oligo_c

    elif oligo1_conc < 6 * oligo2_conc and oligo1_conc > oligo2_conc:

        oligo_c = (oligo1_conc - (oligo2_conc / 2.0)) * 1e-6

        return oligo_c

    elif oligo1_conc == oligo2_conc:

        oligo_c = (oligo1_conc / 2.0) * 1e-6

        return oligo_c

oligo_c = oligo_calculation(oligo1_conc, oligo2_conc)
template_c = (oligo2_conc / 2.0) * 1e-6

#Calculation of the each primer's GC content
def get_primer_gc(primer1, primer2):

    primer1_GC = (sum(1.0 for base in primer1 if base in ['G', 'C']) / primer1_length) * 100
    primer2_GC = (sum(1.0 for base in primer2 if base in ['G', 'C']) / primer2_length) * 100

    return primer1_GC, primer2_GC

primer1_gc, primer2_gc = get_primer_gc(primer1, primer2)
###################################################################################################################################
#Input Sequence Gathering and Characteriztion
###################################################################################################################################

#This function will be used to aquire the fasta file of the input sequence from the Entrez database
def obtainfasta(chrID, start, stop):

    sequence = str(SeqIO.read(Entrez.efetch(db='nucleotide', id=chrID, rettype='fasta', 
    strand='1', seq_start=start, seq_stop=stop), 'fasta').seq)

    return sequence

sequence = obtainfasta(chrID, start, stop)
sequence_length = len(sequence)

#Calculates the GC content of the input DNA sequence
GC_calculation = (sum(1.0 for base in sequence if base in ['G', 'C']) / sequence_length) * 100

###################################################################################################################################
#Melting Temperature Calculation
###################################################################################################################################
'''
This section attempts to correct the model's accuracy discrepancy when computing the melting temperatures of oligos with 
different base sequences but the same GC content. Those types of oligos can differ by as much as 4C or more, however, 
Privalov and Crane-Robinson's model doesn't contain any way to account for those differences. I think the discrepancy 
may be due to the presence of alternating purine and pyrimidine bases. The alternating pairs may break up water binding
to adenine and thymine residues, thus decreasing the entropy cost associated with that water binding. In theory, the longer 
the alternating purine and pyrimidine chain, the more stable the helix is. Conversely, the more long strings of repeating
pyrimidine bases, the less stable the helix.
'''
AT_pairs_p1 = re.findall('(?:A|T){1}(?:A|T){1,}', primer1) #Generates a list of repeating adenine and thymine sequences
AT_pairs_p2 = re.findall('(?:A|T){1}(?:A|T){1,}', primer2) #in each primer. Only includes sequences with two or more repeats.

GC_pairs_p1 = re.findall('(?:G|C){1}(?:G|C){1,}', primer1) #Generates a list of repeating guanine and thymine cytosine
GC_pairs_p2 = re.findall('(?:G|C){1}(?:G|C){1,}', primer2) #in each primer. Only includes sequences with two or more repeats.
'''
The above lists are needed so that their respective influences on helix stability can be measured and quantified, assuming
there is some influence. Don't yet know how best to do that. The functions below attempts to apply those differences to appropriately 
modify the melting temperature.
'''
#These variables take those captured sequences from the lists of pairs above and convert into lengths. The lengths will then
#be used to alter the entropy values during the melting temperature calculation.
p1_AT_len = [len(seq) for seq in AT_pairs_p1]
p1_GC_len = [len(seq) for seq in GC_pairs_p1]

p2_AT_len = [len(seq) for seq in AT_pairs_p2]
p2_GC_len = [len(seq) for seq in GC_pairs_p2]

#Determines the primers' melting temperature based on the method developed by Privalov and Crane-Robinson
#See Privalov, P. L., & Crane-Robinson, C. (2018). https://doi.org/10.1016/j.pbiomolbio.2018.01.007
heat_capacity = .13 #kJ/K,mol-bp
H_A_25, H_T_25 = 25, 25 #kJ/mol-bp
S_A_25, S_T_25 = 72, 70.5 #J/K-mol-bp
H_C_25, H_G_25 = 18.8, 18.8
S_C_25, S_G_25 = 44.7, 44.7
delta_S_trans = gas_constant * math.log(2 / oligo_c)

def p1_melting_calculation(primer1):

    n_A, n_T = primer1.count('A'), primer1.count('T') #Counts the number of each base in the primer
    n_C, n_G = primer1.count('C'), primer1.count('G') #Needed to calculate total enthalpy and entropy
    
    #Total enthalpy and entropy value calculations for determining the melting temperature. Takes the heat capacity into account.
    delta_H = ((H_A_25 + (heat_capacity * (311.5 - 298.15))) * n_A) + ((H_T_25 + (heat_capacity * (311.5 - 298.15))) * n_T) + ((H_C_25 + (heat_capacity * (311.5 - 298.15))) * n_C) + ((H_G_25 + (heat_capacity * (311.5 - 298.15))) * n_G)
    delta_S = ((S_A_25 + (heat_capacity * math.log(311.5 / 298.15))) * n_A) + ((S_T_25 + (heat_capacity * math.log(311.5 / 298.15))) * n_T) + ((S_C_25 + (heat_capacity * math.log((311.5) / 298.15))) * n_C) + ((S_G_25 + (heat_capacity * math.log((311.5) / 298.15))) * n_G)
    '''
Additional melting temperature adjustment to account for the fact that Privalov and Crane-Robinson's model accuracy can
vary significantly depending on the length GC content of the oligo. It tends to underestimate the Tm of low GC oligos and overestimates
the Tm of high GC oligos. The highest accuracy is at a GC% around 50. The values here were measured against the melting temperatures 
from the 92 oligos in Owczarzy's 2004 paper. See Owczarzy, R. et. al (2004). Biochemistry. https://doi.org/10.1021/bi034621r.
'''
    for num in p1_AT_len:
        if num >= 2:
            if primer1_gc < 50:
                delta_S += (.77 * (primer1_length/num)) - (math.log(primer1_gc) * math.log(primer1_length / num)) - ((50 / primer1_gc) + (math.log(primer1_length, 10)))
            else:
                delta_S += (.77 * (primer1_length/num)) - (math.log(primer1_gc) * math.log(primer1_length / num)) - ((50 / primer1_gc) - (math.log(primer1_length, 10)))
    for num in p1_GC_len:
        if num > 1:
            if primer1_gc >= 70:
                delta_S -= -20.6 + ((100 - primer1_gc) * .375) #Compensation for non-linearity at high GC%'s
            else:
                delta_S -= (.43 * (primer1_length/num)) - (math.log(primer1_gc) * math.log(primer1_length / num)) + (50 / primer1_gc)

    tm = (1000 * delta_H) / (delta_S + (gas_constant * (math.log(2/oligo_c)))) - 298.15

    return tm

def p2_melting_calculation(primer2):

    n_A, n_T, n_C, n_G = primer2.count('A'), primer2.count('T'), primer2.count('C'), primer2.count('G')

    delta_H = ((H_A_25 + (heat_capacity * (311.5 - 298.15))) * n_A) + ((H_T_25 + (heat_capacity * (311.5 - 298.15))) * n_T) + ((H_C_25 + (heat_capacity * (311.5 - 298.15))) * n_C) + ((H_G_25 + (heat_capacity * (311.5 - 298.15))) * n_G)
    delta_S = ((S_A_25 + (heat_capacity * math.log(311.5 / 298.15))) * n_A) + ((S_T_25 + (heat_capacity * math.log(311.5 / 298.15))) * n_T) + ((S_C_25 + (heat_capacity * math.log((311.5) / 298.15))) * n_C) + ((S_G_25 + (heat_capacity * math.log((311.5) / 298.15))) * n_G)

    for num in p2_AT_len:
        if num >= 2:
            if primer2_gc < 50:
                delta_S += (.77 * (primer2_length/num)) - (math.log(primer2_gc) * math.log(primer2_length / num)) - ((50 / primer2_gc) + (math.log(primer2_length, 10)))
            else:
                delta_S += (.77 * (primer2_length/num)) - (math.log(primer2_gc) * math.log(primer2_length / num)) - ((50 / primer2_gc) - (math.log(primer2_length, 10)))
    for num in p2_GC_len:
        if num >= 1:
            if primer2_gc >= 70:
                delta_S -= -20.6 + ((100 - primer2_gc) * .375)
            else:
                delta_S -= (.43 * (primer2_length/num)) - (math.log(primer2_gc) * math.log(primer2_length / num)) + (50 / primer2_gc)

    tm = (1000 * delta_H) / (delta_S + (gas_constant * (math.log(2/oligo_c)))) - 298.15

    return tm

primer1_melting_temperature = p1_melting_calculation(primer1)
primer2_melting_temperature = p2_melting_calculation(primer2)
###################################################################################################################################
#Buffer Adjustments
###################################################################################################################################

#PCR buffer conditions
Na = 0
K = 50
Tris = 20
Mg = 1.5
dNTPs = .2

#Adjustments and unit conversions for the chosen buffer conditions
Mon = (Na + K + Tris) / 2.0 #Divide by two to account for the counterion present, e.g. Cl-, SO4-, etc.
mg_adj = Mg * 1e-3 #Converts to mol/L
mon = Mon * 1e-3
dntps = dNTPs * 1e-3 
ka = 3e4 #Association constant for the Mg2+--dNTP complex. Used to calculate the free magnesium in the buffer
mg = (-(ka * dntps - ka * mg_adj + 1.0) + math.sqrt((ka * dntps - ka * mg_adj + 1.0) ** 2 + 4.0 * ka * mg_adj)) / (2.0 * ka)

#Input sequence uses different melting temperature caculation since long strands melt differently than short duplexes
#Equations from Owczarzy 2008 that adjusts melting temperatures according to monovalent ion, magnesium, and dNTP concentrations
#See Owczarzy, R., et al. (2008). Biochemistry, https://doi.org/10.1021/bi702363u
def seq_melting_temperature(GC_calculation, mon):

    seq_tm = 449.15 - (2.60 - (GC_calculation / 100)) * (36.0 - (7.04 * math.log10(mon)))
    
    return round(seq_tm - 273.15, 1)

def primer1_salt_correction(primer1_melting_temperature, primer1_gc, primer1_length, mg, mon):

    a, b, c, d, e, f, g = 3.92e-5, -9.11e-6, 6.06e-5, 1.42e-5, -4.92e-4, 5.25e-4, 8.31e-5 #Slightly altered constants
    '''
The condition calculating the melting temperature if there is no monovalent salt present must be put first, before calculating R.
This is because the R equation requires dividing by the monovalent ion concentration. And of course, you can't divide by zero.
Additionally, although the condition of R > 6 uses the same equation as when no monovalent ions are present, the two conditions
can't be combined. In situations where no magnesium is present, the equation would try to take the natural log of zero, which
throws an error.
'''
    if mon == 0:
        salt2 = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        
        return (1 / salt2) - 273.15

    R = math.sqrt(mg) / mon #Ratio to determine whether monovalent or divalent ions are dominant in their effects on melting temperature

    if R < 0.22:
        salt1 = (1 / (primer1_melting_temperature + 273.15)) + ((4.59e-5 * (primer1_gc / 100)) - 2.90e-5) * math.log(mon) + 8.81e-6 * (math.log(mon)) ** 2
            
        return (1 / salt1) - 273.15

    elif R < 6.0:

        a = 3.92e-5 * (0.843 - (0.352 * math.sqrt(mon) * math.log(mon)))
        d = 1.42e-5 * ((1.279 - 4.03e-3 * math.log(mon)) - 8.03e-3 * (math.log(mon) ** 2))
        g = 8.31e-5 * ((0.486 - 0.258 * math.log(mon)) + 5.25e-3 * (math.log(mon) ** 3)) 

        salt2 = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)

        return (1 / salt2) - 273.15
    
    elif R > 6.0:
        salt2 = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)

        return (1 / salt2) - 273.15

def primer2_salt_correction(primer2_melting_temperature, primer2_gc, primer2_length, mg, mon):

    a, b, c, d, e, f, g = 3.92e-5, -9.11e-6, 6.06e-5, 1.42e-5, -4.92e-4, 5.25e-4, 8.31e-5
    
    if mon == 0:
        salt2 = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        
        return (1 / salt2) - 273.15

    R = math.sqrt(mg) / mon

    if R < 0.22:
        salt1 = (1 / (primer2_melting_temperature + 273.15)) + ((4.59e-5 * (primer2_gc / 100)) - 2.90e-5) * math.log(mon) + 8.81e-6 * (math.log(mon)) ** 2
            
        return (1 / salt1) - 273.15

    elif R < 6.0:

        a = 3.92e-5 * (0.843 - (0.352 * math.sqrt(mon) * math.log(mon)))
        d = 1.42e-5 * ((1.279 - 4.03e-3 * math.log(mon)) - 8.03e-3 * (math.log(mon) ** 2))
        g = 8.31e-5 * ((0.486 - 0.258 * math.log(mon)) + 5.25e-3 * (math.log(mon) ** 3)) 

        salt2 = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)

        return (1 / salt2) - 273.15
    
    elif R > 6.0:
        salt2 = (1 / (primer2_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer2_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer2_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)

        return (1 / salt2) - 273.15

salty_melting_temperature = seq_melting_temperature(GC_calculation, mon)

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

adj_primer1_melting_temperature = primer1_salt_correction(primer1_melting_temperature, primer1_gc, primer1_length, mg, mon)
adj_primer2_melting_temperature = primer2_salt_correction(primer2_melting_temperature, primer2_gc, primer2_length, mg, mon)
final_melting_temperature = DMSO_correction(salty_melting_temperature, GC_calculation)
###################################################################################################################################
#Printing the results
###################################################################################################################################

print('The GC content of primer 1 is ' + str(round(primer1_gc,1)) + '%.' + ' The GC content of primer 2 is ' + str(round(primer2_gc,1)) + '%.')
print('The melting temperatures of your primers are ' + str(round(adj_primer1_melting_temperature, 1)) + 'C for primer 1, and ' + str(round(adj_primer2_melting_temperature, 1)) + 'C for primer 2.')
print('The salt adjusted meliting temperature for the input sequence is ' + str(round(salty_melting_temperature, 1)) + 'C.')
if final_melting_temperature < salty_melting_temperature: 
    print('The final adjusted meliting temperature for the input sequence is ' + str(final_melting_temperature) + 'C.')