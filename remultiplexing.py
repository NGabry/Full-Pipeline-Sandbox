import sys
import os
import subprocess
import csv


def remultiplex(bamfile):


    convert = str('samtools fastq '+bamfile+' > raw_seqs.fastq')
    subprocess.run(convert, shell=True)

    #INSERT BARCODES HERE IF THEY ARE NOT ION XPRESS 1 - 10
    Fwd_bcs = {'01':'CTAAGGTAACGAT',
               '02':'TAAGGAGAACGAT',
               '03':'AAGAGGATTCGAT',
               '04':'TCGCAATTACGAT',
               '05':'CAGAAGGAACGAT',
               '06':'CTGCAAGTTCGAT',
               '07':'TTCGTGATTCGAT',
               '08':'TTCCGATAACGAT',
               '09':'TGAGCGGAACGAT',
               '10':'CTGACCGAACGAT',}

    f = open('fwd_sort.txt', 'w')
    for k in Fwd_bcs:
        fwd_sort = str('cutadapt -g '+Fwd_bcs[k]+' --discard-untrimmed -e 0 -j 0 --overlap 10 -o F'
                   +k+'.fastq raw_seqs.fastq')
        subprocess.run(fwd_sort, shell=True, stdout=f)

    os.system('find . -type f -size 0 -delete')

    #INSERT BARCODES HERE IF THEY ARE NOT ION XPRESS 11 - 20
    Rev_bcs = {'11':'GATTCGAGGA',
               '12':'GAACCACCTA',
               '13':'GTCCGTTAGA',
               '14':'GACACTCCAA',
               '15':'GACCTCTAGA',
               '16':'GTCATCCAGA',
               '17':'GACGAATAGA',
               '18':'GCAATTGCCT',
               '19':'GTCCGACTAA',
               '20':'GATGGATCTG'}

    fastqs = []
    for file in os.listdir('.'):
        if len(file) == 9 and file.endswith('fastq'):
            fastqs.append(file)

    f = open('rev_sort.txt', 'w')
    for fastq in fastqs:
        for k in Rev_bcs:
            rev_sort = str('cutadapt -a '+Rev_bcs[k]+' --discard-untrimmed -e 0 -j 0 --overlap 10 -o '
                           +fastq[0:3]+'R'+k+'.fastq '+fastq)
            subprocess.run(rev_sort, shell=True, stdout=f)

    os.system('find . -type f -size 0 -delete')


    reader = csv.reader(open('20bp_bcs_1-20.csv', 'r', encoding='utf-8-sig'))
    bcs_20bp = {}
    for row in reader:
        k, v = row
        bcs_20bp[k] = v

        fastqs = []
    for fastq in os.listdir('.'):
        if len(fastq) == 12 and fastq.endswith('fastq'):
            fastqs.append(fastq)

    fastqs.sort()


    for fastq in fastqs:
        for k ,v in bcs_20bp.items():
            if fastq[0:6] == k:
                with open(fastq, 'r') as f, open('all_seqs.fastq', 'a+') as out:
                    count = 0
                    for idx, line in enumerate(f.read().splitlines()):
                        count += 1
                        if count == 2:
                            print(v+line, file=out)
                        elif count == 4:
                            print('IIIIIIIIIIIIIIIIIIII'+line, file=out)
                        else:
                            print(line, file=out)
                        if count == 4:
                            count = 0
            else:
                continue

    os.system('mkdir fastqs')
    os.system('mv F*.fastq fastqs| mv fwd_sort.txt fastqs| mv rev_sort.txt fastqs | rm raw_seqs.fastq')


if __name__ == '__main__':
    main()
