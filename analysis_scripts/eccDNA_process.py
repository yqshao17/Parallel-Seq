infile = sys.argv[1] #*.concordant_freq3.2SPLIT-1M.inoneline.txt (circle_finder output)
outfile = sys.argv[2]
with open(infile) as input:
    output=open(outfile+'.tmp', 'w')
    for line in input:
        items=line.split('\t')
        barcode=items[3].split(':')[0]
        umi=items[3].split(':')[1]
            
        if items[0]==items[10] and items[0]==items[20] and items[6]==items[16] and len(items[7])<=12 and len(items[17])<=12 and len(items[27])<=12:
                if (items[6]=='+' and items[26]=='-') or (items[6]=='-' and items[26]=='+'):
                    if items[16]=='+' and items[18]=='second' and int(items[11])<int(items[1]) and  int(items[21])>=int(items[11]) and int(items[22])<=int(items[2]):
                        output.write('\t'.join([items[0],items[11],items[2]])+'\t'+barcode+'\t'+umi+'\t'+items[3]+'\n')
                    elif items[6]=='+' and items[8]=='second' and int(items[1])<int(items[11]) and  int(items[21])>=int(items[1]) and int(items[22])<=int(items[12]):
                        output.write('\t'.join([items[0],items[1],items[12]])+'\t'+barcode+'\t'+umi+'\t'+items[3]+'\n')
                    elif items[16]=='-' and items[18]=='second' and int(items[11])<int(items[1]) and  int(items[21])>=int(items[11]) and int(items[22])<=int(items[2]):
                        output.write('\t'.join([items[0],items[11],items[2]])+'\t'+barcode+'\t'+umi+'\t'+items[3]+'\n')
                    elif items[6]=='-' and items[8]=='second' and int(items[1])<int(items[11]) and  int(items[21])>=int(items[1]) and int(items[22])<=int(items[12]):
                        output.write('\t'.join([items[0],items[1],items[12]])+'\t'+barcode+'\t'+umi+'\t'+items[3]+'\n')
    output.close()
os.system('bash jt.sh %s %s'%(outfile+'.tmp', outfile))