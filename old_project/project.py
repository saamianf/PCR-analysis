print ('PCR plate analysis')

import pandas as pd
import numpy as np

print ("This program is designed to calculate the power",
       "values for your PCR machine's output data for one",
       "96-well plate. Additional plates can then be added",
       "at the end. This data can then be directly inputted",
       "to your platform of choice for data visualization.")

print ('Assumptions: \n\t1. Each primer takes up two full rows.',
       '\n\t2. The first two rows of the initial plate are',
       'the internal control samples for comparison against',
       'all additional primers. \n\t3. For each primer, there',
       'are three consecutive conditions whose power is',
       'calculated together using their averaged dTs.')

plates = int(input('\nHow many plates do you have to compute? (enter integer value) '))
assert type(plates) == int

def avg(primerset):
    avgs = [] 
    for i in range(int(len(primerset)/2)):
        av = (primerset[i] + primerset[i+12])/2
        avgs.append(av)
    return avgs

def calcdT(primer, control):
    dTval = primer - control
    return dTval

def calcdTTs(dTvalue, dTav):
    ddTval = abs(dTvalue - dTav)
    return ddTval

def calcpower(power):
    finalval = 2**(-1*power)
    return finalval

int_con_vals = []
def calc2rows(plate):
    primers = []
    if count == 0:
        int_con_vals.append(avg(plate[0:24]))
        int_con_vals.append(avg(plate[0:24]))
        int_con_vals.append(avg(plate[0:24]))
        int_con_vals.append(avg(plate[0:24]))
        for i in range(24, len(plate), 24):
            primer = plate[i:i+24]
            primers.append(avg(primer))
    if count >= 1:
        for i in range(0, len(plate), 24):
            primer = plate[i:i+24]
            primers.append(avg(primer))            
    dTs = []
    for i in range(len(primers)):
        dTpart = []
        for y in range(len(primers[i])):
            dT = calcdT(primers[i][y], int_con_vals[i][y])
            dTpart.append(dT)
        dTs.append(dTpart)
    ddTs = []
    for i in range(len(dTs)):
        ddTpart = []
        dTavg = (dTs[i][0] + dTs[i][1] + dTs[i][2])/3
        for y in range(len(dTs[i])):
            ddT = calcdTTs(dTs[i][y], dTavg)
            ddTpart.append(ddT)
        ddTs.append(ddTpart)
    finalans = []
    for i in range(len(ddTs)):
        answers = []
        for y in range(len(ddTs[i])):
            answer = calcpower(ddTs[i][y])
            answers.append(answer)
        finalans.append(answers)
    return finalans

count = 0
outputs = []
for i in range(plates):
    pcr_file = input("Enter your raw PCR data file's name:")
    info = (pd.read_excel(pcr_file, sheet_name='0', usecols=[7])).values
    raw = [40 if pd.isna(x) else x for j in info for x in j]
    check = calc2rows(raw)
    outputs.append(check)
    print ('Finished computing for plate #', str(count+1))
    count = count + 1

finaloutputs = [x for j in outputs for x in j]
for i in range(len(finaloutputs)):
    finaloutputs[i].insert(0, str('Primer ' + str(i+2)))

df = pd.DataFrame(finaloutputs, columns = ['Primer', 'Val1', 'Val2', 'Val3', 'Val4', 'Val5', 'Val6',
                  'Val7', 'Val8', 'Val9', 'Val10', 'Val11', 'Val12'])

df.to_csv('outputs.csv')

print ("\nEach primer was calculated against Primer 1 (control). \nCheck 'outputs.csv' for computed Power values!")
