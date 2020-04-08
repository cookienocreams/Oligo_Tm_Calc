import random

nucleotides = 'ATGC'

def sequence_list(nucleotides):
    lst = [] #List which will hold the generated sequences
    count = 0 #Needed to keep track of when to increase the desired GC percent
    gc = 20 #Starting GC percent
    while len(lst) < 14:
        sequence = ''.join(random.choices(nucleotides, k = 10)) #Randomly chooses a base from the nucleotides string k times
        sequence_gc = (sum(1.0 for base in sequence if base in ['G', 'C']) / len(sequence)) * 100
        if round(sequence_gc) == gc: #Checks if the random sequence's GC% equals the target GC%. If yes, it's added to the list.
            lst.append(sequence)
            count += 1
            if count % 2 == 0: #Increases GC percent after two sequences are added
                gc += 10 #This value changes with different lengths. It's 10 for 10, 20, 30 bp oilgos, 8 for 25 bp, and 6 and 7 for 15 bp
    return lst

sequence = sequence_list(nucleotides)
#Writes the generated sequences into a new file. Adds sequences one line at a time.
with open('Random Sequences.txt', 'w') as filehandle:
    for seq in sequence:
        filehandle.write('%s\n' % seq)