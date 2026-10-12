file = open('data/galex.dat','r')
lines = file.readlines()[1:]
file.close()

new = open('data/galex_new.dat','w')
for line in lines:
    line = line.split()
    new.write(line[0]+'\t'+line[1]+'\t')
    new.write('\t'.join(line[25:29]))
    new.write('\n')
new.close()