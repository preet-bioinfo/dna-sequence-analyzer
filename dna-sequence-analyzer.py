#DNA sequence analyzer
#this program analyzer basic properties of DNA sequence

#------function 1: count nucleotides-----
#count how many times each bases (A,T,G,C)appear in DNA sequence
def count_nucleotides(dna):
    dna=dna.upper()
    count={
        'A':dna.count('A'),#count Adenine
        'T':dna.count('T'),#count thymine
        'G':dna.count('G'),#count guanine
        'C':dna.count('C')#count cytosine
    }
    return count#return the dictinoary of counts
#test it
sequence="ATGCATGCatg"
result=count_nucleotides(sequence)
print(result)



#--------function 2:GC content-------
#GC content = percentage of G and C in the sequence
#high GC = more stable DNA (G-C has 3 hydrogen bonds and A-T has 2)
def gc_content(dna):
    dna=dna.upper()
    gc=dna.count('G')+dna.count('C') #total G and C bases
    total=len(dna)                   #total length of sequence
    percentage=(gc/total)*100        #calculate percentage
    return round(percentage,2)       #round 2 decimal places

#----test----
sequence="ATGCATGCatg" #sample DNA sequence (mixed case to test upper()

result=count_nucleotides(sequence)
print("Nucleotide count"),result

print("GC content:",gc_content(sequence,),"%")


#-----function 3: complement strand -------
#each bases in DNA pairs with its complement
# A pairs with T , T pairs with A ,G pairse with C, C pairs with G
def complement(dna):
    dna=dna.upper()
    pairs={'A':'T','T':'A','G':'C','C':'G'}  #base pairing rules
    comp=''  #empty string to build complement
    for base in dna:            #go through eah base one by one
        comp=comp+pairs[base]   #add its complement
    return comp


#---test it ---
print("complement:",complement(sequence))


#------function 4 : Reverse complement-----
#DNA is antiparallel - so the acttual complementary strand
#is read in reverse direction
#used in PCR primer design, BLAST , gene finding
def reverse_complement(dna):
    comp=complement(dna)  #first get the complement
    rev_comp=comp[::-1]   #then reverse it using slicing
    return rev_comp

#-test it -
print("reverse complement:",reverse_complement(sequence))

#----------function 5: Transcription-----------
#DNA is tanscribed into mRNA in the nucleus
#only difference :Thymine (T) is replaced by Uracil(U)
def transcribe(dna):
    dna=dna.upper()
    mrna=dna.replace('T','U')    #replace every T with U
    return mrna

#--test it --
print("mRNA:",transcribe(sequence))



#-----function 6: Translation------
#mRNA is read in group of 3 bases called codons
#each codon codes for one amino acid
#UAA , UAG ,UGA = stop codons(end of protein)

def translate(mrna):
    codon_table = {
        'AUG': 'Met', 'UUU': 'Phe', 'UUC': 'Phe',
        'UUA': 'Leu', 'UUG': 'Leu', 'UCU': 'Ser',
        'UCC': 'Ser', 'UAU': 'Tyr', 'UAC': 'Tyr',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGG': 'Trp',
        'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln',
        'CAG': 'Gln', 'AAU': 'Asn', 'AAC': 'Asn',
        'AAA': 'Lys', 'AAG': 'Lys', 'GAU': 'Asp',
        'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala',
        'GCG': 'Ala', 'GGU': 'Gly', 'GGC': 'Gly',
        'GGA': 'Gly', 'GGG': 'Gly', 'CCU': 'Pro',
        'CCC': 'Pro', 'ACU': 'Thr', 'ACC': 'Thr',
        'UAA': 'STOP', 'UAG': 'STOP', 'UGA': 'STOP'
    }
    protein = []
    for i in range(0, len(mrna) - 2, 3):  # read every 3 bases
        codon = mrna[i:i+3]                # extract codon
        amino_acid = codon_table.get(codon, '?')  # lookup amino acid
        if amino_acid == 'STOP':           # stop if stop codon
            break
        protein.append(amino_acid)
    return ' - '.join(protein)

# Test it
mrna_seq = transcribe(sequence)
print("Protein:", translate(mrna_seq))


#updated by preet-bioinfo