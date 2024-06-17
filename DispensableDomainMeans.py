# This script randomly selects N amount of values from Data_list and calculates the mean; this repeats 1000 times
import random
File = open("Rep3Values.txt", "r")
Data_read = File.read()
Data_split = Data_read.split('\n')
Data_list = [eval(number) for number in Data_split]

B1B2loop_area = 57
B3a2loop_area = 17
Cterminal_area = 28

Allmeans = list()

for i in range(1000):
    Picks = random.sample(Data_list, B3a2loop_area)
    Mean = sum(Picks) / B3a2loop_area
    Allmeans.append(Mean)
    print(Mean)
    i += 1

print(Allmeans)
