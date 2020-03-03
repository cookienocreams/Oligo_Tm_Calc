import Bio
from Bio import Entrez, SeqIO
import re, math
Entrez.email = "cookienocreams@outlook.com"

#These are the ideal variables that will form the standards for the recommendations
ideal_duplex_stabilty = 95 #In Celcius
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

#Captures the first and last base in each primer so that the energies associated with opening the helix can be calculated
p1_initial_base = primer1[0]
p1_terminal_base = primer1[-1]

p2_initial_base = primer2[0]
p2_terminal_base = primer2[-1]

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

initial_base = sequence[0]
terminal_base = sequence[-1]

#Calculates the GC content of the input DNA sequence
GC_calculation = (sum(1.0 for base in sequence if base in ['G', 'C']) / sequence_length) * 100

###################################################################################################################################
#Determines terminal base compensation parameters
###################################################################################################################################

#Calculates the entropic and enthalpic compensation for the ends of each primer
def p1_terminal_comp(p1_initial_base, p1_terminal_base):

    AT_h, AT_s, GC_h, GC_s = 2.3, 4.1, .1, -2.8 #h in kcal/mol, s in eu
    dh, ds = 0, 0

    dh_init = (dh + AT_h if p1_initial_base in ['A', 'T'] else dh + GC_h)
    ds_init = (ds + AT_s if p1_initial_base in ['A', 'T'] else ds + GC_s)

    dh_term = (dh + AT_h if p1_terminal_base in ['A', 'T'] else dh + GC_h)
    ds_term = (ds + AT_s if p1_terminal_base in ['A', 'T'] else ds + GC_s)

    dh_total = dh_init + dh_term
    ds_total = ds_init + ds_term

    return dh_total, ds_total

def p2_terminal_comp(p2_initial_base, p2_terminal_base):
    
    AT_h, AT_s, GC_h, GC_s = 2.3, 4.1, .1, -2.8 #h in kcal/mol, s in eu
    dh, ds = 0, 0

    dh_init = (dh + AT_h if p2_initial_base in ['A', 'T'] else dh + GC_h)
    ds_init = (ds + AT_s if p2_initial_base in ['A', 'T'] else ds + GC_s)

    dh_term = (dh + AT_h if p2_terminal_base in ['A', 'T'] else dh + GC_h)
    ds_term = (ds + AT_s if p2_terminal_base in ['A', 'T'] else ds + GC_s)

    dh_total = dh_init + dh_term
    ds_total = ds_init + ds_term

    return dh_total, ds_total

p1_term_h, p1_term_s = p1_terminal_comp(p1_initial_base, p1_terminal_base)
p2_term_h, p2_term_s = p2_terminal_comp(p2_initial_base, p2_terminal_base)

###################################################################################################################################
#Melting Temperature Calculation and Salt Adjustments
###################################################################################################################################

#Determines the primers' melting temperature based on the method developed by Privalov and Crane-Robinson
#See Privalov, P. L., & Crane-Robinson, C. (2018). https://doi.org/10.1016/j.pbiomolbio.2018.01.007
def p1_melting_calculation(primer1, p1_term_h, p1_term_s):
    
    heat_capacity = .13 #kJ/K,mol-bp
    H_A_25, H_T_25 = 25, 24.75 #kJ/mol-bp
    S_A_25, S_T_25 = 72, 71.55 #J/K-mol-bp
    H_CG_25, S_CG_25 = 18.9, 44.7
    delta_S_trans = gas_constant * math.log(2 / oligo_c)
    n_A, n_T = primer1.count('A'), primer1.count('T')
    n_CG = primer1.count('C') + primer1.count('G')

    delta_H = ((H_A_25 + (heat_capacity * (311.5 - 298.15))) * n_A) + ((H_T_25 + (heat_capacity * (311.5 - 298.15))) * n_T) + ((H_CG_25 + (heat_capacity * (311.5 - 298.15))) * n_CG) + p1_term_h
    delta_S = ((S_A_25 + (heat_capacity * math.log(311.5 / 298.15))) * n_A) + ((S_T_25 + (heat_capacity * math.log(311.5 / 298.15))) * n_T) + ((S_CG_25 + (heat_capacity * math.log(311.5 / 298.15))) * n_CG) + p1_term_s

    tm = (1000 * delta_H) / (delta_S + (gas_constant * (math.log(2 / oligo_c)))) - 298.15

    return tm

def p2_melting_calculation(primer2, p2_term_h, p2_term_s):
    
    heat_capacity = .13 #kJ/K,mol-bp
    H_A_25, H_T_25 = 25, 24.75 #kJ/mol-bp
    S_A_25, S_T_25 = 72, 71.55 #J/K-mol-bp
    H_CG_25, S_CG_25 = 18.9, 44.7
    delta_S_trans = gas_constant * math.log(2 / oligo_c)
    n_A, n_T = primer2.count('A'), primer2.count('T')
    n_CG = primer2.count('C') + primer2.count('G')

    delta_H = ((H_A_25 + (heat_capacity * (311.5 - 298.15))) * n_A) + ((H_T_25 + (heat_capacity * (311.5 - 298.15))) * n_T) + ((H_CG_25 + (heat_capacity * (311.5 - 298.15))) * n_CG) + p2_term_h
    delta_S = ((S_A_25 + (heat_capacity * math.log(311.5 / 298.15))) * n_A) + ((S_T_25 + (heat_capacity * math.log(311.5 / 298.15))) * n_T) + ((S_CG_25 + (heat_capacity * math.log((311.5) / 298.15))) * n_CG) + p2_term_s

    tm = (1000 * delta_H) / (delta_S + (gas_constant * (math.log(2 / oligo_c)))) - 298.15

    return tm

primer1_melting_temperature = p1_melting_calculation(primer1, p1_term_h, p1_term_s)
primer2_melting_temperature = p2_melting_calculation(primer2, p2_term_h, p2_term_s)

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

    a, b, c, d, e, f, g = 3.92e-5, -9.11e-6, 6.26e-5, 1.42e-5, -4.82e-4, 5.25e-4, 8.31e-5
    
    if mon == 0:

        salt2 = (1 / (primer1_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer1_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer1_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
        
        return (1 / salt2) - 273.15

    R = math.sqrt(mg) / mon

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

    a, b, c, d, e, f, g = 3.92e-5, -9.11e-6, 6.26e-5, 1.42e-5, -4.82e-4, 5.25e-4, 8.31e-5
    
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
adj_primer1_melting_temperature = primer1_salt_correction(primer1_melting_temperature, primer1_gc, primer1_length, mg, mon)
adj_primer2_melting_temperature = primer2_salt_correction(primer2_melting_temperature, primer2_gc, primer2_length, mg, mon)

#Additional melting temperature adjustment to account for the fact that Privalov and Crane-Robinson's model accuracy can
#vary significantly depending on the GC content of the oligo. It tends to underestimate the Tm of low GC oligos and overestimates
#the Tm of high GC oligos. The highest accuracy is at a GC% around 50. Oligo length is also taken into account here.
def gc_tm_adjustment(adj_primer1_melting_temperature, adj_primer2_melting_temperature, primer1_gc, primer2_gc):
    
    low_gc, mid_gc, high_gc = .166667, .376667, .466667 #Correction factors for three GC% tiers

    if primer1_gc < 45:

        p1_tm = float(adj_primer1_melting_temperature + (low_gc * (52 - primer1_gc)))

    elif primer1_gc >= 45 and primer1_gc <= 65:
        if primer1_length >=25:
            if round(primer1_gc) in [47,48,49,50,51,52]:
                
                p1_tm = float(.959 * (adj_primer1_melting_temperature - (mid_gc * (primer1_gc - 50))))
        else:
            p1_tm = float((adj_primer1_melting_temperature - (mid_gc * (primer1_gc - 50))))

    else:
        if primer1_length <= 15:

            high_gc = .266667
            p1_tm = float(adj_primer1_melting_temperature - (high_gc * (primer1_gc - 48))) 
        elif primer1_length <= 25:
            
            p1_tm = float(adj_primer1_melting_temperature - (high_gc * (primer1_gc - 48)))
        else:

            high_gc = .566667
            p1_tm = float(adj_primer1_melting_temperature - (high_gc * (primer1_gc - 48))) 
        return p1_tm

    if primer2_gc < 45:

        p2_tm = float(adj_primer2_melting_temperature + (low_gc * (52 - primer2_gc)))

    elif primer2_gc >= 45 and primer2_gc <= 65:
        if primer2_length >=25:
            if round(primer2_gc) in [47,48,49,50,51,52]:

                p2_tm = float(.959 * (adj_primer2_melting_temperature - (mid_gc * (primer2_gc - 50))))
        else:
            p2_tm = float((adj_primer2_melting_temperature - (mid_gc * (primer2_gc - 50))))

    else:
        if primer2_length <= 15:

            high_gc = .266667
            p2_tm = float(adj_primer2_melting_temperature - (high_gc * (primer2_gc - 48))) 

        elif primer2_length <= 25:
            
            p2_tm = float(adj_primer2_melting_temperature - (high_gc * (primer2_gc - 48)))

        else:

            high_gc = .566667
            p2_tm = float(adj_primer2_melting_temperature - (high_gc * (primer2_gc - 48))) 

    return p1_tm, p2_tm

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

final_tm_p1, final_tm_p2 = gc_tm_adjustment(adj_primer1_melting_temperature, adj_primer2_melting_temperature, primer1_gc, primer2_gc)
final_melting_temperature = DMSO_correction(salty_melting_temperature, GC_calculation)

###################################################################################################################################
#Printing the results
###################################################################################################################################

print('The GC content of primer 1 is ' + str(round(primer1_gc,1)) + '%.' + 
' The GC content of primer 2 is ' + str(round(primer2_gc,1)) + '%.')
print('The melting temperatures of your primers are ' + str(round(final_tm_p1, 1)) + 
'C for primer 1, and ' + str(round(final_tm_p2, 1)) + 'C for primer 2.')
print('The salt adjusted meliting temperature for the input sequence is ' + str(round(salty_melting_temperature, 1)) + 'C.')
if final_melting_temperature < salty_melting_temperature: 
    print('The final adjusted meliting temperature for the input sequence is ' + str(final_melting_temperature) + 'C.')