import sys

with open(sys.argv[2]) as qc:
    next(qc)
    for lin in qc:
        lin = lin.strip("\n")
        lin = lin.split("\t")
        clean = int(lin[3])
with open(sys.argv[3],'w') as f:
    with open(sys.argv[1]) as fu:
        title = next(fu)
        f.write(title)
        for line in fu:
            line = line.strip("\n")
            line = line.split("\t")
            support = int(line[6])
            if vars().has_key('support'):
                if clean == 0:
                    norm = support
                else:
                    norm = int(support * 50 / clean)
                if norm >= 200:
                    line[6] = str(norm)
                    line[7] = "Positive"
                    for word in line:
                        f.write(word + '\t')
                    f.write("\n")
                else:
                    line[6] = str(norm)
                    line[7] = "Negetive"
                    #for word in line:
                    #    f.write(word + '\t')
                    #f.write("\n")
