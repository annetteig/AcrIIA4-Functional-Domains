# This script matches the sequence of each repair template to the plasmid reads and counts every occurrence
import json

template_count = dict()
sequences = list()

# Check if miseq plasmid read contains template sequence
def match_maker(template, miseq_read):
    '''if miseq_read.startswith(template):
        print("Found a match: " + template + "\n- " + miseq_read)
        return True
    return False
    '''
    if miseq_read in template:
        return True
    return False


# Load template file containing all repair template sequences
def load_templates():
    temp_file = open("Acr_ALT_ALL_dels.txt", "r")
    for template in temp_file.readlines():
        template_count[template.strip().upper()] = 0

# Load miseq plasmid reads
def load_seqs():
    global sequences
    seqs_file = open("TESTAGAIN.txt")
    for sequence in seqs_file.readlines():
        sequences.append(sequence.strip().upper())
load_templates()
load_seqs()

# Match template and read sequences and count occurrences
i = 0
for template in template_count.keys():
    if i == 80:
        pass
    original_template = template
    print("Checking: " + str(i))
    for sequence in sequences:
        if match_maker(template,sequence):
            template_count[original_template] += 1
    i+=1

# Save to file
with open("TESTcounts.txt","w") as result_file:
    result_file.write(json.dumps(template_count))

