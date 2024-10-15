# This script intakes the file with RT sequences and read counts, and 'annotates' it with the start and end index of the deletion each repair template represents
sequences = list()

# Annotations = file containing RT sequence and start and end indexes
annotations = open("annotations.txt", "r").readlines()
counts = open("5foaB_counts.txt", "r").readlines()

print(len(counts))
print(len(annotations))

# Match sequences and append deletion annotations
i = 0
for count in counts:
    count = count.upper()
    count_seq = count[0:200].strip()
    print(str(i))
    for annotation in annotations:
        annotation = annotation.upper()
        annotation_seq = annotation.split(" ")[2].strip()
        if (count_seq in annotation_seq):
            sequences.append(count.strip() + " " + " " + annotation)
    i+=1

print(len(sequences))
with open("5foaB_annotatedcounts.txt", "w") as results_file:
    for i in sequences:
        results_file.write(i)
