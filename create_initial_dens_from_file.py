import sys

ogfilename = sys.argv[1]
newfilename = sys.argv[2]

ogfile = open(ogfilename, "r")

lastline = ogfile.readlines()[-1]

ogfile.close()

newfile = open(newfilename, "w")
newfile.write(lastline)

newfile.close()
