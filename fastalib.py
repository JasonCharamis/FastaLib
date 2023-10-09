
from natsort import natsorted
import argparse
import re, os, core


def isfile(input_file, field=0):
    # Function to check if input is a file or a string

    if os.path.isfile(input_file):
        with open(input_file, "r") as file:
            # Check if single or multi-column file, if latter is true select the one with gene list (good for using in associative lists of genes)
            lines = file.readlines()
            glist = []

            for line in lines:
                num_columns = len(line.strip().split('\t'))  # Split the first line into columns based on the delimiter

                if num_columns == 1:
                    glist.append(line.strip())
                    
                elif num_columns > 1:
                    glist.append(line.split('\t')[field].strip())

            sorted_list = natsorted( glist )

            return sorted_list

    elif isinstance(input_file, str):
        glists = input_file

    return glists



class FASTA:

    __slots__ = ['id','seq']
    
    def __init__ (self, header, sequence):
        self.id = str(header)
        self.seq = str(sequence)
        

    def __str__ ( self ) :
        return f">{self.id}\n{self.seq}"

    
    def fasta_parser(fasta_file, sprint = False):
        fasta_instances = []

        with open(fasta_file, 'r') as fasta:
            lines = fasta.readlines()
            seqid = ""
            sequence = ""
            entries = {}

            for line in lines:
                line = line.strip('\n')  # Remove newline character
               
                if line.startswith('>'):
                    seqid = re.sub(".*name=|\-0000\d|\.t\d","", line[1:])
                    sequence = ""
                        
                elif not re.search("--|#", line):  # Ignore lines containing "--" or "#"
                    sequence += re.sub("\*$|\.$","",line)

                if seqid and sequence:
                    entries[seqid]=sequence

        for seqids, sequences in entries.items():
            fasta_instance = FASTA(seqids, sequences)
            fasta_instances.append(fasta_instance)

            #        sorted_list = natsorted(fasta_instances, key=lambda instance: getattr(instance, id))
                    
        return fasta_instances
   
        
    def fasta_sizer ( fasta_file): ## sizer of fasta sequences, takes dictionary as input ##
        lst = []
        fasta_instances = FASTA.fasta_parser(fasta_file)

        # Define a custom sorting key function
        for fasta_instance in fasta_instances:
            fasta_size = '\t'.join([fasta_instance.id, str(len(fasta_instance.seq))])
            lst.append (fasta_size)
           
        sorted_list = natsorted(lst, key=lambda x: x[1], reverse=False) # sort based on sequence length
       
        return sorted_list

       
     def fasta_extract(fasta_file, gene_list):
        subset = []
        fasta_instances = FASTA.fasta_parser(fasta_file)

        if isinstance(isfile(gene_list),list):       
            subset = list(set([ fasta_instance for fasta_instance in fasta_instances for g in isfile(gene_list) if re.search(g, fasta_instance.id) ]))

        else:
            subset = list(set([ fasta_instance for fasta_instance in fasta_instances if re.search(isfile(gene_list), fasta_instance.id) ]))

        if len(subset) > 0:            
            return subset
        
        else:
            print ( "Gene list is empty or not provided." )

            
    def fasta_remove ( fasta_file, gene_list ):
        subset = []
        fasta_instances = FASTA.fasta_parser(fasta_file)

        if isinstance(isfile(gene_list),list):       
            subset = list(set([ fasta_instance for fasta_instance in fasta_instances for g in isfile(gene_list) if not re.search(g, fasta_instance.id) ]))

        else:
            subset = list(set([ fasta_instance for fasta_instance in fasta_instances if not re.search(isfile(gene_list), fasta_instance.id) ]))
                
        if len(subset) > 0:            
            return subset
        
        else:
            print ( "Gene list is empty or not provided." )

    def translate(fasta_file):
        codon2aa = {
           'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
           'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
           'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
           'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
           'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
           'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
           'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
           'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
           'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
           'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
           'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
           'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
           'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
           'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
           'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
           'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W'
       }

       
        fasta_instances = FASTA.fasta_parser(fasta_file)

        protein_fasta = []
       
        for fasta_instance in fasta_instances:
            fasta_instance.seq = fasta_instance.seq.upper() ## Convert CDS to uppercase

            if not re.search ("A|T|G|C", fasta_instance.seq):
                print ( fasta_instance.id + " is not a nucleotide sequence. Found other than A,T,G,C." )

            else:
                aminoacids = []
    
                for i in range(0,len(fasta_instance.seq),3):
                    codon = fasta_instance.seq[i:i+3]
                    
                    if codon in codon2aa:
                        aminoacids.append(codon2aa[codon])

                    else:
                        aminoacids.append("N")

            protein_sequence = FASTA ( fasta_instance.id, ''.join(aminoacids) )   
            protein_fasta.append( protein_sequence )

        if len(protein_fasta ) > 0 :
            return protein_fasta

        else:
            print ( "Provided fasta file is empty.")



## Implementation ##

def kargs():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Library for efficiently manipulating fasta files.')

    ## options to parse and convert multi-line fasta to one-line fasta
    parser.add_argument('-f','--fasta', type=str, help='FASTA file')
    parser.add_argument('-o','--one_line', action="store_true", help='Option to extract FASTA sequences.')

    ## options to extract or remove genes
    parser.add_argument('-e','--extract', action="store_true", help='Option to extract FASTA sequences.')
    parser.add_argument('-r','--remove', action="store_true", help='Option to remove FASTA sequences.')
    parser.add_argument('-g','--gene_list', type=str, help='Gene list to extract or FASTA sequences.')
    parser.add_argument('--field', type=int, help='Field to extract geneids from gene_list.')

    ## option to compute size of fasta sequences
    parser.add_argument('-s','--size', action="store_true", help='Option to print sizes of FASTA sequences.')

    ## option to convert CDS to PEP 
    parser.add_argument('-tr','--cds2pep', action="store_true", help='Option to print sizes of FASTA sequences.')   
        
    args = parser.parse_args()

    if not any(vars(args).values()):
        parser.print_help()
        print("Error: No arguments provided.")

    return parser.parse_args()


def main():
    
    parser = argparse.ArgumentParser(description='Library for efficiently manipulating fasta files.')
    args = kargs()
    
    if args.fasta:

        inp = re.sub (".aa$|.fa$|.faa$|.fna$|.fasta$|.1l$","",args.fasta)

        if args.one_line:
            with open(f"{inp}.fasta.1l", "w") as f:
                for out in FASTA.fasta_parser(args.fasta):
                    print ( out, file = f )

        elif args.size:
            with open(f"{inp}.fasta.sizes", "w") as f:
                for out in FASTA.fasta_sizer(args.fasta):
                    print ( out, file = f )

        elif args.extract:
            if args.gene_list:
                with open(f"{inp}.{args.gene_list}.extracted.fasta", "w") as f:
                    for out in FASTA.fasta_extract(args.fasta, args.gene_list):
                        print ( out, file = f )
            else:
                print ( "Please provide a gene list.")
                
        elif args.remove:
            with open(f"{inp}.{args.gene_list}.removed.fasta", "w") as f:
                for out in FASTA.fasta_remove(args.fasta, args.gene_list, args.field):
                    print ( out, file = f )
               
        elif args.cds2pep:
            with open(f"{inp}.pep.fasta", "w") as f:
                for out in FASTA.translate(args.fasta):
                    print ( out, file = f )              
    
    else:
        print("Please select a potential FASTA operation.")
               

if __name__ == "__main__":
    main()


