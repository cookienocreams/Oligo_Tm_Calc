import re, math

###########################################################################################################################################
#Prompts to get the information required to begin the analysis
#Primer information is input here
###########################################################################################################################################

input_primer_list = input('Input primer sequences separated by commas here: ').replace(" ", "").split(',')
oligo_conc = float(input('Input total single strand oligo concentration in uM here: '))
#This section gathers information from the user about the conditions to primer will be used in so that a more accurate Tm can be calculated
Mono = float(input('Input total monovalent ion concentration in mM here: '))
Mg = float(input('Input Mg concentration in mM here: '))
dNTPs = float(input('Input dNTP concentration in mM here: '))
oligo_c = (oligo_conc) * 1e-6 #Adjusts primer concentration to mol/L
list_of_melting_temperatures = [] #Melting temperatures will be added to this list as they are calculated
list_of_gc_content = [] #The percent GC is added to this list as they are calculated

###########################################################################################################################################
#Primer Melting Temperature Calculator
###########################################################################################################################################

def calculate_melting_temperature(input_primer_list,oligo_c,Mono,Mg,dNTPs):
    for sequence in input_primer_list:#Loops through the input sequences and calculates each one's melting temperature
        #List of all the nearest-neighbor pairs
        NN_pairs_list = [['AA', 'TT'], ['AC', 'TG'], ['AG', 'TC'], ['AT'], ['CA', 'GT'], ['CC', 'GG'], ['CG'], ['GA', 'CT'], ['GC'], ['TA']]
        purine, pyrimidine = ['A','G'], ['C','T'] #Lists used to convert primer sequence into R/Y purine/pyrimidine nucleotide codes
        gas_constant = 1.98720425864083 #cal*K−1*mol−1

        primer_length = len(sequence)
        primer_RY_sequence_str = ''.join(['R' if base in purine else 'Y' for base in sequence]) #Converts primers into their purine/pyrimidine nucleotide codes
        #Captures the first and last pair of bases in each primer so that the energies associated with opening the helix can be calculated
        primer_initial_bases, primer_terminal_bases = primer_RY_sequence_str[0:2], primer_RY_sequence_str[-2:primer_length]
        primer_NN_list = re.findall('.{1,2}', sequence) + re.findall('.{1,2}', sequence[1:]) #Splits input primers into doublets for easy NN matching below
        primer_gc = (sum((1.0 for base in sequence if base in ['G', 'C'])) / primer_length) * 100 #Calculation of each primer's GC content

        #Modified Nearest Neighbor Parameters based on SantaLucia's 1997 paper Thermodynamics and NMR of internal G·T mismatches in DNA. https://doi.org/10.1021/bi962590c
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

        #Lists of each nearest neighbor's entropy and enthalpy values that will be used to calculate the thermodynamic values of each primer
        Delta_S = [AA_delta_s,AC_delta_s,AG_delta_s,AT_delta_s,CA_delta_s,CC_delta_s,CG_delta_s,GA_delta_s,GC_delta_s,TA_delta_s]
        Delta_H = [AA_delta_h,AC_delta_h,AG_delta_h,AT_delta_h,CA_delta_h,CC_delta_h,CG_delta_h,GA_delta_h,GC_delta_h,TA_delta_h]

        #End Compensation parameters adjusted for primer length and changes in monovalent ion concentration
        YY_constant_one = (-.0235 * primer_length) + .1273
        YY_constant_two = (.1639 * primer_length) - .895
        YR_constant_one = (-.296 * math.log(primer_length)) + .5058
        YR_constant_two = (2.0303 * math.log(primer_length)) - 3.4594
        YY_length_adjust_eq = YY_constant_one * math.log(Mono) + YY_constant_two
        YR_length_adjust_eq = YR_constant_one * math.log(Mono) + YR_constant_two
        YY_h,RY_h = 0.959547884222969 - YY_length_adjust_eq, 0.849386567650066 - YR_length_adjust_eq

        def primer_thermodynamic_sum(Delta_H,Delta_S,NN_pairs_list): #Calculates the total enthalpy and entropy of each sequence
            count, primer_total_dh, primer_total_ds = 0, 0, 0 #Needed to track values of each NN pair
            NN_pair = NN_pairs_list[count] #Indexed to calculate the values of each set of NN in the nearest-neighbor list
            for NN_pair in NN_pairs_list:
                primer_total_dh += Delta_H[count] * sum((1.0 for NN in primer_NN_list if NN in NN_pair))#Calculates the number of occurances of each NN in the primer sequence
                primer_total_ds += Delta_S[count] * sum((1.0 for NN in primer_NN_list if NN in NN_pair))#and multiplies the number of occurances by the associated delta h/s value
                count += 1 #Adds one to the count to keep track of which NNs have been calculated already
            return primer_total_dh, primer_total_ds
        
        primer_total_dh, primer_total_ds = primer_thermodynamic_sum(Delta_H,Delta_S,NN_pairs_list)
        primer_dh_init = (primer_total_dh + YY_h if primer_initial_bases in ['YY','RR'] else primer_total_dh + RY_h)
        primer_dh_total = (primer_dh_init + YY_h if primer_terminal_bases in ['YY','RR'] else primer_dh_init + RY_h)

        #Determines the melting temperature of the primer at 1M monovalent salt, still needs adjustment to account for actual salt concentration + Mg2+ and dNTPs
        primer_melting_temperature = (1000 * primer_dh_total) / (primer_total_ds + (gas_constant * (math.log(oligo_c)))) - 273.15

        #Equations from Owczarzy 2008 that adjusts melting temperatures according to monovalent ion, magnesium, and dNTP concentrations
        #See Owczarzy, R., et al. (2008). Biochemistry, https://doi.org/10.1021/bi702363u
        def primer_salt_correction(primer_melting_temperature):
            #Adjustments and unit conversions for the chosen buffer conditions
            Mon = Mono / 2.0 #Divide by two to account for the counterion present, e.g. Cl-, SO4-, etc.
            mg_adj = Mg * 1e-3 #Converts to mol/L
            mon = Mon * 1e-3 #Converts to mol/L
            dntps = dNTPs * 1e-3 #Converts to mol/L
            ka = 3e4 #Association constant for the Mg2+--dNTP complex. Used to calculate the free magnesium in the buffer
            mg = (-(ka * dntps - ka * mg_adj + 1.0) + math.sqrt((ka * dntps - ka * mg_adj + 1.0) ** 2 + 4.0 * ka * mg_adj)) / (2.0 * ka)

            #General constants and constants to regulate balance between Mg2+ and Na+ ions
            a, b, c, d, e, f, g = 3.919e-5, -2.88e-5, 3.603e-5, 2.322e-5, -3.507e-4, 4.711e-4, 6.52e-5 
            const_a,const_b,const_c,const_d,const_e,const_f,const_g,const_h = -0.1156,-2.0725,-0.1445,6.247e-3,6.131e-3,0.0314,0.5308,4.563e-3
        
            if mon == 0:
                salty = (1 / (primer_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
                return round((1 / salty) - 273.15, 1) #This equation converts the above calculation back into a useable melting temperature
            R = math.sqrt(mg) / mon #R is used to determine whether salt or magnesium is the primary factor affecting melting temperature
            if R < 0.22: #This equation is used if salt is the primary factor
                if Mono > 900:
                    salty = (1 / (primer_melting_temperature + 273.15)) + ((4.29e-5 * (primer_gc / 100)) - 3.95e-5) * math.log(mon) + 9.40e-6 * (math.log(mon)) ** 2
                    #return round((primer_melting_temperature + ((1 / salty) - 273.15)) / 2, 1)
                    return round((1 / salty) - 273.15, 1)
                salty = (1 / (primer_melting_temperature + 273.15)) + ((4.29e-5 * (primer_gc / 100)) - 3.95e-5) * math.log(mon) + 9.40e-6 * (math.log(mon)) ** 2
                return round((1 / salty) - 273.15, 1)
            #This equation is used if there is a complex balance between salt and magnesium
            elif R < 6.0:
                a = 3.92e-5 * (const_a * math.log(mg) - (const_b * math.sqrt(mon) * math.log(mon)))
                d = 2.32e-5 * ((const_c * math.log(mg) - const_d * math.log(mon)) - const_e * (math.log(mon) ** 2))
                g = 6.52e-5 * ((const_f * math.log(mg) - const_g * math.log(mon)) + const_h * (math.log(mon) ** 3))
                salty = (1 / (primer_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
                return round((1 / salty) - 273.15, 1)
            #This equation is used if magnesium is the primary factor
            elif R > 6.0:
                salty = (1 / (primer_melting_temperature + 273.15)) + a + (b * math.log(mg)) + ((primer_gc / 100) * (c + d * math.log(mg))) + (1 / (2.0 * (primer_length - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
                return round((1 / salty) - 273.15,1)
            
        salt_adj_primer_melting_temperature = primer_salt_correction(primer_melting_temperature)

        list_of_melting_temperatures.append(salt_adj_primer_melting_temperature) #Adds calculated melting temperature to temperatures list
        list_of_gc_content.append(primer_gc) #Adds calculated percent GC to GC list

#Printing the results
def display_results(list_of_gc_content, list_of_melting_temperatures):
    count = 0
    primer_tracker = 1
    for gc_percent in list_of_gc_content:
        primer_gc_percent = list_of_gc_content[count]
        primer_Tm = list_of_melting_temperatures[count]
        print('The GC content of primer '+ str(primer_tracker) +' is ' + str(round(primer_gc_percent, 1)) + '% ' +'and the melting temperature is ' + str(primer_Tm) + 'C')
        count += 1
        primer_tracker += 1

def main():
    calculate_melting_temperature(input_primer_list,oligo_c,Mono,Mg,dNTPs)
    display_results(list_of_gc_content, list_of_melting_temperatures)

if __name__=='__main__':
    main()
