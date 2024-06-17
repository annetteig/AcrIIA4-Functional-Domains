# This script matches the sequence of each repair template to the genomic reads and counts how often a match is found

import json

template_count = dict()
sequences = list()

# Match repair template to genomic miseq_reads
def match_maker(template, miseq_read, starts_w_hom, ends_w_hom):
    if starts_w_hom and ends_w_hom:
        if template == miseq_read:
            return True
    elif starts_w_hom:
        miseq_read = miseq_read[0:len(template)]
        if template == miseq_read:
            return True
    elif ends_w_hom:
        miseq_read = miseq_read[len(miseq_read) - len(template) : len(miseq_read)]
        if template == miseq_read:
            return True
    elif template.upper() in miseq_read.upper():
        return True

    return False

# Load file containing all repair template sequences
def load_templates():
    temp_file = open("Acr_ALT_ALL_dels.txt", "r")
    for template in temp_file.readlines():
        template_count[template.strip()] = 0

# Load file containing all genomic read sequences
def load_seqs():
    global sequences
    seqs_file = open("Replicate1Glucose_genomicreads.txt")
    for sequence in seqs_file.readlines():
        sequences.append(sequence.strip())

# Match and count how many times each repair template sequence perfectly appears in genomic reads
i = 0
for template in template_count.keys():
    original_template = template
    starts_w_hom = template[0].islower()
    ends_w_hom = template[-1].islower()
    if starts_w_hom or ends_w_hom:
        template = template.replace("a", "")
        template = template.replace("t", "")
        template = template.replace("g", "")
        template = template.replace("c", "")
    print("Checking: " + str(i))
    for sequence in sequences:
        if match_maker(template,sequence,starts_w_hom,ends_w_hom):
            template_count[original_template] += 1
    i+=1


