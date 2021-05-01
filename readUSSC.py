# -*- coding: utf-8 -*-
import sys
import re
# annovarda verilen komut bu bundan test ediyorum.
#  ./annotate_variation.pl -out ex1 -build hg19 example/ex1.avinput humandb/


# default genannotation yapacak filter ve region gelince opsiyonel olacak
genomebinsize = 10000

# bu da opsiyonel olacak şu an bu değer çünkü sadece genanno yapıyoruz
neargene = 1000#for upstream/downstream annotation of variants, specify the distance threshold between variants and genes


def checkRecordValid(record,argsDic):
    if argsDic['dbtype'] == 'refGene':
        if len(record) != 16 and len(record) == 15:
            record.insert(0,0)
            return True
        elif len(record) == 16:
            return True

        print("invalid record in " + argsDic['dbloc'] + "expecting 15 or 16 tab-delimited fields in refGene file line \n")
        return False


def createDbDicts(lines,argsDic):
    lineCnt = 0
    dbtype = argsDic['dbtype']
    dbloc = argsDic['dbloc']

    genedb = {} # bunlar veritabanındaki bilgilerle doldurulacak
    geneidmap = {} # bunlar veritabanındaki bilgilerle doldurulacak
    name2count = {}
    cdslen = {} # bunlar veritabanındaki bilgilerle doldurulacak
    mrnalen = {} # bunlar veritabanındaki bilgilerle doldurulacak
    iscoding = {}
    genecount = 0
    ncgenecount = 0

    badgene = {}

    file = open("resultUSCS.txt","a")
    for line in lines:
        lineCnt +=1
        record = re.split("\t",line)
        if checkRecordValid(record,argsDic):
            if dbtype == 'refGene':
                name,chr,dbstrand,txstart,txend,cdsstart,cdsend,exoncount,exonstart,exonend,id,name2,cdsstartstat,cdsendstat,exonframes = record[1:]

                # bu değerlerin int olması gerekiyor.
                txstart=int(txstart)
                txend=int(txend)
                cdsstart=int(cdsstart)
                cdsend=int(cdsend)
                exoncount=int(exoncount)

		        #handle situations where the same transcript is mapped to several chromosomes or regions (for example, NM_019105 is mapped to chr6, chr6_cox_hap1, chr6_qbl_hap2; NM_002538 is mapped to chr5 positive and negative strand and also in chr5_h2_hap1)
                if re.search(r'hap\d+$',chr): # 3115 tane buluyor annovar ile doğru
                    continue

                chr = re.search(r'[^chr]',chr).group() #some genomes like zebrafish does not start with chr in their chromosome names.

                name = "#".join((name, chr, str(txstart))) #for transcript-based gene definition, it is possible for the same transcript to map to multiple locations, therefore the "name" is not unique.

                if not(dbstrand == "+" or dbstrand == "-"):
                    raise Exception("Error: invalid dbstrand information found in %s (dbstrand has to be + or -): %d\n" % (dbloc, lineCnt))

                exonSIntList = exonstart.split(",")
                exonSIntList.pop() # son "," yüzünden boş eleman geliyor onu pop ediyoruz.
                exonSIntList = map(int, exonSIntList)

                exonEIntList = exonend.split(",")
                exonEIntList.pop()  # son , yüzünden boş eleman geliyor onu pop ediyoruz.
                exonEIntList = map(int, exonEIntList)

                if not int(exoncount) == len(exonSIntList):
                    raise Exception("Error: invalid record found in %s (exoncount discordance): %d vs %d\n" % (dbloc, int(exoncount), len(exonSIntList)))

                if not len(exonSIntList) == len(exonEIntList):
                    raise Exception("Error: invalid record found in %s (exonstart and exonend count discordance): %d vs %d\n" % (dbloc, len(exonSIntList), len(exonEIntList)))

                txstart += 1
                cdsstart += 1

                #LOGIC here:
                #first calcluate mRNA length, and if the transcript maps to multiple locations with discordant mRNA length, only consider the leftmost chromosome and leftmost coordinate (because the FASTA file is sorted in this manner)

                cdsLength = 0
                mrnaLength = 0

                for i in range(len(exonSIntList)):
                    mrnaLength += exonEIntList[i] - exonSIntList[i] +1

                for i in range(len(exonSIntList)):
                    if cdsstart >= exonSIntList[i] and cdsstart <= exonEIntList[i]:
                        if(cdsend <= exonEIntList[i]):
                            cdsLength = cdsend - cdsstart +1
                            break
                        else :
                            cdsLength += exonEIntList[i] - cdsstart + 1
                            continue
                    if cdsLength and cdsend < exonSIntList[i]:
                        raise Exception("FATAL ERROR: impossible scenario for %s ,n %s  (cdsend is less than exon start)", name, dbloc)
                    elif cdsLength and cdsend <= exonEIntList[i]:
                        cdsLength += exonEIntList[i] - exonSIntList[i] + 1
                        break
                    elif cdsLength and cdsend - exonSIntList[i]:
                        cdsLength += exonEIntList[i] - exonSIntList[i] + 1

                file.write(name+"\n")

                if cdsstart != cdsend + 1: # coding gene
                    if name in mrnalen and mrnalen[name] != mrnaLength:
                        print("WARNING %s occurs more than once in %s with different mRNA length. The first occurences with identical mRNA length will be uesd in analysis.", name, dbloc)
                        continue
                    if name in cdslen and cdslen[name] != cdsLength:
                        print()
                        continue
                    if not name2 in iscoding:
                        iscoding[name2] = 0
                    iscoding[name2] += 1
                else: # none coding gene
                    1

                cdslen[name] = cdsLength
                mrnalen[name] = mrnaLength


                bin1 = int((txstart - neargene) / genomebinsize)
                bin2 = int((txend + neargene) / genomebinsize)

                for nextbin in range(bin1,bin2):
                    genedbKey = chr + str(nextbin)
                    if(genedbKey in genedb):
                        genedb[genedbKey].append([name,dbstrand,txstart,txend,cdsstart,cdsend,exonSIntList,exonEIntList,name2])
                    else:
                        genedb[genedbKey]= [name,dbstrand,txstart,txend,cdsstart,cdsend,exonSIntList,exonEIntList,name2]
                geneidmap[name] = name2
                genecount +=1
                if( name2 in name2count):
                    name2count[name2] += 1
                else:
                    name2count[name2] = 1
                if cdsstart == cdsend +1 :
                    ncgenecount +=1
        else:
            continue
    file.close()

    #for key in genedb:
    #    newgenedb = {}
    #    geneinfo = genedb[key]
    #    for




def processDBFile(argsDic):
    dbloc = argsDic['dbloc']
    lines = []
    with open(dbloc) as dbFile:
        for line in dbFile:
            lines.append(line.rstrip())
    dbtype = argsDic['dbtype']
    record = re.split("\t",lines[0])
    createDbDicts(lines,argsDic)


def checkArgumetns(args):
    checkedArgs = {}
    argsCnt = len(args)
    if argsCnt>1:
        args = args[1:]
        argCount = argsCnt/2

        if argCount > 0:
            for i in range(argsCnt/2):
                if args[i*2][:2] == "--":
                    checkedArgs[args[i*2][2:]] = args[i*2+1]
        else:
            print("please give arguments.")
    else:
        print("please give argument.")

    return checkedArgs


def main():
    argsDic = checkArgumetns(sys.argv)
    processDBFile(argsDic)


if(__name__== "__main__"):
    main()
