import csv
import sys

ogfilename = sys.argv[1]
newfilename = sys.argv[2]

with open(ogfilename, 'r')as ogfile:
    csvFile = csv.reader(ogfile)
    lastline = list(csvFile)[-1]
    
newfile = open(newfilename, "w")
#the starting index of 2 omits the x0 and dx0 that will be at the beginning of the row in the input
for i in range(2, len(lastline)-1):
    newfile.write(lastline[i])
    newfile.write(", ")
    
newfile.write(lastline[-1])
newfile.close()