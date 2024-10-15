# This script generates all possible repair templates (including deletions of START and STOP codons)
# Previous Acr sequence:
#GeneA = "aataaaagtatcaacaaaaaattgttaatatacctctatactttaacgtcaaggagaaaaaaccccggattctagaactagtggatcccccgggaaaaaaATGAACATCAACGACCTCATCCGCGAGATAAAGAACAAAGACTATACCGTGAAGCTAAGCGGGACTGACAGCAACAGTATTACACAGCTAATTATTCGTGTCAACAACGACGGAAACGAGTATGTAATCAGCGAGAGCGAGAATGAGAGCATTGTGGAGAAATTCATAAGTGCCTTCAAGAATGGATGGAATCAGGAGTATGAGGACGAGGAGGAGTTCTATAACGACATGCAGACCATCACTCTTAAGAGTGAGCTTAACTGAnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn"
# Alternate Acr sequence (first copy)
GeneA = "aataaaagtatcaacaaaaaattgttaatatacctctatactttaacgtcaaggagaaaaaaccccggattctagaactagtggatcccccgggaaaaaaATGAACATCAACGACCTCATCCGCGAGATAAAGAACAAAGACTATACCGTGAAGCTAAGCGGGACTGACAGCAACTCCATAACGCAGCTGATAATACGCGTAAATAACGACGGAAACGAGTATGTAATCAGCGAGAGCGAGAATGAGAGCATTGTGGAGAAATTCATAAGTGCCTTCAAGAATGGATGGAATCAGGAGTATGAGGACGAGGAGGAGTTCTATAACGACATGCAGACCATCACTCTTAAGAGTGAGCTTAACTGAnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn"
# Second Acr copy
GeneB = "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnATGAATATTAATGATTTGATTAGAGAAATTAAAAATAAGGATTACACTGTTAAATTGTCTGGTACAGATTCTAATTCAATCACTCAATTGATCATCAGAGTTAATAATGATGGTAATGAATACGTTATTTCTGAATCAGAAAACGAATCTATCGTTGAAAAGTTTATTTCAGCTTTTAAAAACGGTTGGAACCAAGAATACGAAGATGAAGAAGAATTTTACAATGATATGCAAACTATTACATTGAAATCAGAATTAAATTAAtatataaactcatttacttatgtaggaataaagagtatcatctttcaaacgcccagcggtagtacaattcaaagtagtaggtaccaatggtagtactagt"

file = open("Acr_ALTERNATIVE_dels.txt", "a")

homology_len = 100
gene_length = len(GeneA) - (homology_len*2)
codons = gene_length // 3

RT_list = list()
deletion_length = 3

# Make every possible in-frame deletion starting at the START codon of GeneA
for num_codons_del in range((codons) - 1):
    num_RT = codons - ((deletion_length//3) - 1)
    for index in range(num_RT):
        deletion_pt = (index*3) + homology_len
        geneA_segment = GeneA[deletion_pt - 100 : deletion_pt]
        geneB_segment = GeneB[deletion_pt + deletion_length : deletion_pt + deletion_length + 100]
        #print("For deletion length " + str(deletion_length) + " make a deletion at index " + str(deletion_pt))
        #print(geneA_segment)
        #print(geneB_segment)
        RT = geneA_segment + geneB_segment
        RT_list.append(RT)
        print( str(deletion_pt-100) + " " + str((deletion_pt + deletion_length - 99)) + " " + RT )
        #print(RT)

    deletion_length += 3

for element in RT_list:
    file.write(element)
    file.write('\n')
file.close()
