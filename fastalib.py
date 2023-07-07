from Bio import SeqIO
from Bio.Seq import Seq
import re

## writes a fasta dictionary to file ##
def write_file (dict,filename):   
    with open(filename, "w") as outfile:
        for id in dict.keys():
            outfile.write(">"+id+"\n"+dict[id]+"\n")

## converts fasta to 1l and writes into file ##            
def fasta_1l(fasta):
    genes={}
    for seq_record in SeqIO.parse(fasta, "fasta"):
        name=re.sub(".*/| .*","",str(seq_record.id))
        name=re.sub(".p","_p/",name)
        sequence=re.sub("\*|\--","",str(seq_record.seq))
        genes[name]=sequence
        write_file(genes,str(fasta+".1l"))
    return genes
    

## custom fasta parser with base python ##
def fasta_parser ( fasta_file ):
    seqs = {}
    header = ""
    sequence = ""   

    with open ( fasta_file, "r" ) as file:
        for fasta in file:
            fasta=fasta.strip ("\n")
        
            if fasta.startswith(">"):
                header=re.sub(">","",fasta)
                seqs[header]=sequence
      
            elif not re.search ("\--",fasta):
                sequence += fasta
                seqs[header] = sequence
                sequence = ""
    return seqs


## sizer of fasta sequences, takes dictionary as input ##
def fasta_sizer ( fasta_file, mode = True):
    sizer = {}

    for id,seq in fasta_1l(fasta_file):
        sizer[id] = len( seq )
        sort = dict ( sorted ( sizer.items(), key = lambda x:x[1], reverse = mode ) )

    for pr,size in sort.items():
        print ( pr+"\t",size )
        
    return


def fasta_filter ( fasta_file, filter_list ):
    list = {}

    with open ( filter_list, "r" ) as file:
        for l in file:
            l = re.sub ( "\n", "", l )
            if l not in list.keys():
                list[l]=1

    for id,seq in fasta_parser(fasta_file).items():
        if id not in list.keys():
            print (">"+id + "\n" + fasta_parser(fasta_file)[id] )

    return


def extract_fasta ( fasta_file, extract_list ):
    list = {}

    with open ( extract_list, "r" ) as file:
        for l in file:
            l = re.sub ( "\n", "", l )
            if l not in list.keys():
                list[l]=1

    for id,seq in fasta_parser(fasta_file).items():
        if id in list.keys():
            print (">"+id + "\n" + fasta_parser(fasta_file)[id] )

    return


def extract_fasta_headers ( fasta_file ):
    for id,seq in fasta_parser(fasta_file).items():
        if id in list.keys():
            id = re.sub(" .*", "", id )
            print (id )

    return
