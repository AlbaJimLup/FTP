#!/usr/bin/env python

# AUTHOR: Johan Zicola
# DATE: 2018-07-19

import argparse
import os
import random
import sys

'''
	USAGE: This script contains several functions to generate fasta 
        and fastq files. It can generate a random fasta file, generate 
        random fastq files (paired end or single end) of specified length
        and insert size (for paired end). The script can also generate
        single-end and paired-end fastq files based on a specified fasta
        reference file. Coverage and read size can be specified.
'''



#Functions

def generate_DNA(sequence_size):
    seq = []
    while sequence_size != 0:
        nucleotide = random.choice('ACGT')
        seq.append(nucleotide)
        sequence_size -= 1

    return ''.join(seq)


def generate_fasta(name_seq, sequence_size):
    print('>'+name_seq)
    print(generate_DNA(sequence_size))


def reverse_complement(sequence):
    my_dna = str(sequence.upper())
    new_str1 = str.replace(my_dna, 'A','t')
    new_str2 = str.replace(new_str1, 'T','a')
    new_str3 = str.replace(new_str2, 'C','g')
    new_str4 = str.replace(new_str3, 'G','c')
    new_str5 = new_str4.upper()
    new_str6 = new_str5[::-1]
    return new_str6


def generate_fastq(sequence, header="@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:1297 1:N:0:TGACCAAT", splitFQs=False):
    # Head contains: Instrument name, run ID, flowcell ID, flowcell lane,
    # tile number within flowcell lane, x-coordinate, y-coordinate,
    # member of a pair (1 or 2), filtered or not (Y or N), control bits, index sequence
    #header = "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:1297 1:N:0:TGACCAAT"
    # Assuming best quality Phred score for Illumina 1.8+
    
    quality = len(sequence) * 'I'
    
    if splitFQs:
        return("{}\n{}\n+\n{}\n".format(header, sequence, quality))
    else:
        print(header)
        print(sequence)
        print('+')
        print(quality)


def generate_random_fastq_SE(sequence_size, nb_seq):
    index = 1

    while nb_seq > 0:
        sequence = generate_DNA(sequence_size)
        header = "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:"+str(index)+" 1:N:0:TGACCAAT"
        generate_fastq(sequence, header)
        nb_seq -=1
        index += 1

        
def generate_random_fastq_PE(sequence_size, nb_seq):

    nb_seq_initial = nb_seq
    
    index = 1
    while nb_seq > 0:
        sequence = generate_DNA(sequence_size)
        header = "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:"+str(index)+"#0/1"
        generate_fastq(sequence, header)
        nb_seq -=1
        index += 1

    #Reset nb_seq
    nb_seq = nb_seq_initial
    
    index = 1
    while nb_seq > 0:
        sequence = generate_DNA(sequence_size)
        header = "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:"+str(index)+"#0/2"
        generate_fastq(sequence, header)
        nb_seq -=1
        index += 1


def generate_mapped_fastq_SE(ref_fasta, sequence_size, coverage):
    #Retrieve sequences from fasta file
    fasta =  open(ref_fasta,"r")
    fasta = fasta.read()
    fasta = fasta.split("\n")

    list_seq = []

    for seq in fasta:
        if seq[:1] != ">":
            seq = seq.strip()
            list_seq.append(seq.upper())

    #Remove empty line if there is
    list_seq = [_f for _f in list_seq if _f]
    


    #Create all possible kmers for each sequence of the fasta file
    
    kmers = []
    
    while coverage > 0:
        start = 0 + coverage
        for dna in list_seq:
            while start + sequence_size < len(dna):
                    kmer = dna[start:start+sequence_size]
                    kmers.append(kmer)
                    start += sequence_size
        coverage -= 1

        
    index = 1
    for kmer in kmers:
        header = "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:"+str(index)+" 1:N:0:ATGT"
        generate_fastq(kmer, header)
        index += 1



# Function generates a set of paired-end reads. The first read maps to the
# forward strand and the second on the reverse strand. The output needs to 
# be processed to separate reads 1 and 2 (e.g. with tail and head functions)

def generate_mapped_fastq_PE(ref_fasta, sequence_size, insertion_size, coverage, fqRoute, fqUnit, splitFQs):
    #Retrieve sequences from fasta file
    fasta =  open(ref_fasta,"r")
    fasta = fasta.read()
    fasta = fasta.split("\n")
    # sys.stderr.write(str((fasta)))
    # sys.stderr.write(str(len(fasta))+"\n")

    list_seqID = []
    list_seq = []
    indiv_seq = ""

    # for seq in fasta:
    #     if seq[:1] != ">":
    #         seq = seq.strip()
    #         list_seq.append(seq.upper())
    seqID_line = False
    for i, seq in enumerate(fasta):
        # account for fasta sequences that are split over multiple lines. 
        # check to see line is ID or sequence
        # sys.stderr.write("\n{}, {}\n".format(i, seq))
        if seq[:1] == ">":
            seqID_line = True
            # will use list of seq IDs for headers later
            # these seq IDs appear at the same index as their respective sequences
            #    in list_seq, which is why we use "i" within the while loop below
            list_seqID.append(seq[1:])
        else:
            seqID_line = False
        
        # if sequence, concatenate all lines for that sequence (in case it is split)
        if not seqID_line:
            seq = seq.strip()
            indiv_seq += seq
        
        # if reach an ID line, will have concatenated all lines of previous sequence, can append
        #   also must account for final sequence which does not have a succeeding ID line
        #   finally, do not include the initialized empty string from before the for loop (i > 0)
        if (seqID_line or i == len(fasta)-1) and i > 0:
            list_seq.append(indiv_seq.upper())
            indiv_seq = ""
            
    # sys.stderr.write("\n\n{}\n".format(list_seq))
    
    # sys.stderr.write("{}  {}\n".format(len(list_seqID), len(list_seq)))
    
    #Remove empty line if there is
    list_seq = [_f for _f in list_seq if _f]
    # sys.stderr.write(str(len(list_seq))+"\n")
    # import numpy as np
    # from statistics import mode
    # sys.stderr.write(str(np.mean([len(x) for x in list_seq]))+"\n")
    # sys.stderr.write(str(np.var([len(x) for x in list_seq]))+"\n")
    # sys.stderr.write(str(np.max([len(x) for x in list_seq]))+"\n")
    # sys.stderr.write(str(np.min([len(x) for x in list_seq]))+"\n")
    # sys.stderr.write(str(mode([len(x) for x in list_seq]))+"\n")
    
    # Initialize the 2 lists
    kmers1 = []
    kmers2 = []
    
    kmers1_coords = []
    kmers2_coords = []
    
    # trying to fix erroneous insertion size
    insert_wrong = insertion_size - (sequence_size*2)
    # just change the original value directly so that the input will use the correct
    #    definition of insert size (ie the full length between the adapters, instead
    #    of just the length of bps that occurs between the for and rev reads)
    insertion_size = insertion_size - (sequence_size*2)
    coverage_real = coverage
    
    # Get read sequences starting from index coverage-1 so that the last coverage
    # the index 0 (coverage = 1 - 1). Reads containing read1 (kmer1), an insert, 
    # read2 (kmer2). The coverage will switch the reads of 1 bp along the sequence.
    list_seqID_counts = []
    while coverage > 0:
        # start = -1 + coverage
        for i, dna in enumerate(list_seq):
            # sys.stderr.write("{} {} {}\n".format(i, len(dna), list_seqID[i]))
            start = -1 + coverage
            k1for = 0
            k2for = 0
            while (start + sequence_size*2 + insertion_size) < len(dna):
              
                # add seqID to a list as many times as kmers are generated for that seq
                # later will use this when assigning the headers
                list_seqID_counts.append( list_seqID[i] )
                
                # generate 0 or 1 randomly to decide which will be for or rev
                randOrient = random.randint(0,1)
                if randOrient == 0:
                    # if random is 0, kmer1 is for, kmer2 is rev
                    kmer1 = dna[start:start+sequence_size]
                    kmers1.append(kmer1)
                    
                    kmer2 = dna[start+insertion_size+sequence_size:start+insertion_size+sequence_size*2]
                    kmer2 = reverse_complement(kmer2)
                    kmers2.append(kmer2)
                    
                    kmers1_coords.append((start, start+sequence_size))
                    # if start == 99 and start+sequence_size == 175:
                    #     sys.stderr.write("{} {} {}\n".format(start, start+sequence_size, list_seqID[i]))
                    kmers2_coords.append((start+insertion_size+sequence_size, start+insertion_size+sequence_size*2))
                    k1for += 1
                elif randOrient == 1:
                    # if random is 1 kmer2 is for, kmer1 is rev
                    kmer2 = dna[start:start+sequence_size]
                    kmers2.append(kmer2)
                    
                    kmer1 = dna[start+insertion_size+sequence_size:start+insertion_size+sequence_size*2]
                    kmer1 = reverse_complement(kmer1)
                    kmers1.append(kmer1)
                    
                    kmers2_coords.append((start, start+sequence_size))
                    kmers1_coords.append((start+insertion_size+sequence_size, start+insertion_size+sequence_size*2))
                    # if start+insertion_size+sequence_size == 99 and start+insertion_size+sequence_size*2 == 175:
                    #   sys.stderr.write("{} {} {}\n".format(start+insertion_size+sequence_size, start+insertion_size+sequence_size*2, list_seqID[i]))
                    k2for += 1
                ktot = k1for + k2for
                # start += (sequence_size*2)+insertion_size
                
                #start += ((sequence_size*2)+insertion_size)/2
                start += sequence_size
                # start += 1
            #sys.stderr.write("tot = {}\nk1for = {} ({}%)\nk2for = {} ({}%)\n\n".format(ktot, k1for, float(k1for)/ktot, k2for, float(k2for)/ktot))
            
            # Jesse's modification to fix insert size but now coverage is wrong            
            # while (start + sequence_size*2 + insert_wrong) < len(dna):
            #     kmer1 = dna[start:start+sequence_size]
            #     kmers1.append(kmer1)
            # 
            #     kmer2 = dna[start+insert_wrong+sequence_size:start+insert_wrong+sequence_size*2]
            #     kmer2 = reverse_complement(kmer2)
            #     kmers2.append(kmer2)

            #    start += insertion/sequence_size#1#(sequence_size*2)+insert_wrong
                
        # coverage -= 1
        coverage -= 2
        # coverage -= (sequence_size*2)+insertion_size
    
    # if indicated, write fastq to separate files for for and rev
    # instead of to stdout, as by default
    if splitFQs:
        import gzip
        # r1 = gzip.open(ref_fasta.replace(".fasta", ".R1.fq.gz").replace(".fa", ".R1.fq.gz"), "wt")
        # r2 = gzip.open(ref_fasta.replace(".fasta", ".R2.fq.gz").replace(".fa", ".R2.fq.gz"), "wt")
        r1 = gzip.open("{}/{}.{}PE_{}x.R1.fq.gz".format(fqRoute, fqUnit, sequence_size, coverage_real), "wt")
        r2 = gzip.open("{}/{}.{}PE_{}x.R2.fq.gz".format(fqRoute, fqUnit, sequence_size, coverage_real), "wt")
    
    
    # To link the paired reads, they need to contain in their name the same
    # Y coordinate (based on observation 35517.end1.fq and 35517.end2.fq headers
    index = 1
    for i, kmer1 in enumerate(kmers1):
        header = "@{}:1:R1_{}-{}_R2_{}-{}:1:1:0:{}#0/1".format(list_seqID_counts[index-1], 
          kmers1_coords[i][0], kmers1_coords[i][1], 
          kmers2_coords[i][0], kmers2_coords[i][1], 
          index)
        # if kmers1_coords[i][0] == 99 and kmers1_coords[i][1] == 175:
        #     sys.stderr.write("{}\n".format(header))
        # header = "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:"+str(index)+"#0/1"
        if splitFQs:
            r1.write( generate_fastq(kmer1, header, splitFQs) )
        else:
            generate_fastq(kmer1, header)
        index += 1

    index = 1
    for i, kmer2 in enumerate(kmers2):
        header = "@{}:1:R1_{}-{}_R2_{}-{}:1:1:0:{}#0/2".format(list_seqID_counts[index-1], 
          kmers1_coords[i][0], kmers1_coords[i][1], 
          kmers2_coords[i][0], kmers2_coords[i][1], 
          index)
        # header = "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:"+str(index)+"#0/2"
        if splitFQs:
            r2.write( generate_fastq(kmer2, header, splitFQs) )
        else:
            generate_fastq(kmer2, header)
        index += 1
    
    
    if splitFQs:
        r1.close()
        r2.close()



def main():
    if len(sys.argv) == 1 or  sys.argv[1] == '-h':
        print('''
    fastq_generator
    
    Work with Python2.7
    
    Author: Johan Zicola (johan.zicola@gmail.com)
    
    Date: 2017-09-09
        
    usage: fastq_generator.py [-h] function

    functions available: generate_fasta, generate_random_fastq_SE, generate_random_fastq_PE, generate_mapped_fastq_SE, generate_mapped_fastq_PE. 

    For more information for each function, enter function name followed by "-h" (e.g.  python fastq_generator.py generate_fasta -h')''')
        sys.exit()
    
    if sys.argv[1] == "generate_fasta":
        #Display help if no argument of "-h"
        if len(sys.argv) == 2 or  sys.argv[2] == '-h':
            print('''
    generate_fasta creates a de novo fasta file of chosen
    name and size. The sequences are generated randomly.
    
    usage: generate_fasta name_seq sequence_size
    
    arguments:
    name_seq = name of the sequence (string)
    sequence_size = length of the sequence in bp (integer)''')
            sys.exit()
        else:
            name_seq = str(sys.argv[2])
        
        #Test if sequence_size argument is given and if integer
        try:
            sys.argv[3] = int(sys.argv[3])
        except IndexError:
            sys.exit("Error: sequence_size missing")
        except ValueError:
            sys.exit("Error: sequence_size should be an integer")
        
        #Call the function
        sequence_size = sys.argv[3]
        generate_fasta(name_seq, sequence_size)
    
    elif sys.argv[1] == "generate_random_fastq_SE":
        if len(sys.argv) == 2 or  sys.argv[2] == '-h':
            print('''
    generate_random_fastq_SE generates random sets of single-end reads
    
    usage: generate_random_fastq_SE sequence_size nb_seq
                     
    arguments:
    sequence_size = length of the sequence in bp (integer)       
    nb_seq = number of sequences required (integer)''') 
            sys.exit()
        else:
            sequence_size = sys.argv[2]
            try:
                sequence_size = int(sequence_size)
            except ValueError:
                sys.exit("Error: sequence_size should be an integer")  
        
        try: 
            sys.argv[3] = int(sys.argv[3])
        except IndexError:
            sys.exit("Error: nb_seq missing")
        except ValueError:
            sys.exit("Error: nb_seq should be an integer")        

        #Call the function
        nb_seq = sys.argv[3]
        generate_random_fastq_SE(sequence_size, nb_seq)
 
    elif sys.argv[1] == "generate_random_fastq_PE":
        if len(sys.argv) == 2 or  sys.argv[2] == '-h':
            print('''
    generate_random_fastq_PE generates random sets of paired-end reads
    
    usage: generate_random_fastq_SE sequence_size nb_seq
                     
    arguments:
    sequence_size = length of the sequence in bp (integer)       
    nb_seq = number of sequences required (integer)''') 
            sys.exit()
        else:
            sequence_size = sys.argv[2]
            try:
                sequence_size = int(sequence_size)
            except ValueError:
                sys.exit("Error: sequence_size should be an integer")  
        
        try: 
            sys.argv[3] = int(sys.argv[3])
        except IndexError:
            sys.exit("Error: nb_seq missing")
        except ValueError:
            sys.exit("Error: nb_seq should be an integer")        

        #Call the function
        nb_seq = sys.argv[3]
        generate_random_fastq_PE(sequence_size, nb_seq)
     
    elif sys.argv[1] == "generate_mapped_fastq_SE":
        if len(sys.argv) == 2 or  sys.argv[2] == '-h':
            print('''
    generate_mapped_fastq_SE generates sets of single-end reads from a reference fasta file
    
    usage: generate_mapped_fastq_SE ref_fasta sequence_size coverage
                     
    arguments:
    ref_fasta = reference file in fasta format (fasta file)
    sequence_size = length of the sequence in bp (integer)
    coverage = number of reads covering each bp (integer)''')
            sys.exit()
        else:
            ref_fasta = sys.argv[2]
        
        try: 
            sys.argv[3] = int(sys.argv[3])
        except IndexError:
            sys.exit("Error: sequence_size missing")
        except ValueError:
            sys.exit("Error: sequence_size should be an integer") 
    
        
        try: 
            sys.argv[4] = int(sys.argv[4])
        except IndexError:
            sys.exit("Error: coverage missing")
        except ValueError:
            sys.exit("Error: coverage should be an integer")
        
        #Call the function
        sequence_size = sys.argv[3]
        coverage = sys.argv[4]
        generate_mapped_fastq_SE(ref_fasta, sequence_size, coverage)
        
    elif sys.argv[1] == "generate_mapped_fastq_PE":
        if len(sys.argv) == 2 or  sys.argv[2] == '-h':
            print('''
    generate_mapped_fastq_PE generates sets of paired-end reads from a reference fasta file
    
    usage: generate_mapped_fastq_PE ref_fasta sequence_size  insertion_size coverage [outRoute] [outFQunit]
     
    arguments:
    ref_fasta = reference file in fasta format (fasta file)
    sequence_size = length of the sequence in bp (integer)       
    insertion_size = distance between the paired reads in bp (integer)
    coverage = number of reads covering each bp (integer)''')
            sys.exit()
        else:
            ref_fasta = sys.argv[2]
        
        # test argument sequence_size
        try: 
            sys.argv[3] = int(sys.argv[3])
        except IndexError:
            sys.exit("Error: sequence_size missing")
        except ValueError:
            sys.stderr.write(sys.argv[3])
            sys.exit("Error: sequence_size should be an integer") 
        
        # test argument insertion_size
        try: 
            sys.argv[4] = int(sys.argv[4])
        except IndexError:
            sys.exit("Error: insertion_size missing")
        except ValueError:
            sys.exit("Error: insertion_size should be an integer") 
        
        # test argument coverage
        try: 
            sys.argv[5] = int(sys.argv[5])
        except IndexError:
            sys.exit("Error: coverage missing")
        except ValueError:
            sys.exit("Error: coverage should be an integer")      
        
        # test argument split fastqs
        try:
            # if you want it to split the fq files into for and rev,
            # include 2 additional agruments:
            #   directory for outputs
            #   unit number ($i in the bash script)
            # the other relevant variables come from the other agruments
            fqRoute = sys.argv[6]
            fqUnit = sys.argv[7]
            splitFQs = True
        except IndexError:
            # by default splitFQs is false, writes to stdout
            fqRoute = None
            fqUnit = None
            splitFQs = False
        
        #Call the function
        sequence_size = sys.argv[3]
        insertion_size = sys.argv[4]
        coverage = sys.argv[5]
        generate_mapped_fastq_PE(ref_fasta, sequence_size, insertion_size, coverage, fqRoute, fqUnit, splitFQs)    
        
        
        
        
    
if __name__ == "__main__":
    sys.exit(main())

    
    
