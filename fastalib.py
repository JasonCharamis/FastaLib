
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
    else:
        print("Please select a potential FASTA operation.")
               

if __name__ == "__main__":
    main()


