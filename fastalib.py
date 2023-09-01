from natsort import natsorted
import argparse
import re, os


def check_isfile ( input_file ): ## Function to check if input is a file or a string
    if os.path.isfile(input_file):
        glist = []
        
        with open(input_file, "r") as file:
            lines = file.readlines()
            for line in lines:
                glist.append(line.strip())
        return glist
        
    elif isinstance(input_file, str):
        glists = input_file
        return glists


class FASTA:

    __slots__ = ['id','seq']
    
    def __init__ (self, header, sequence):
        self.id = header
        self.seq = sequence


    def __str__ ( self ) :
        return f">{self.id}\n{self.seq}"
    

    def fasta_parser ( fasta_file ): ## opens fasta file and creates new FASTA instances 
        fasta_instances = []

        with open ( fasta_file , 'r' ) as fasta:
            lines = fasta.readlines()

            for line in lines:
                line = line.strip('\n')
                sequence = ""

                if line.startswith('>'):
                    id = re.sub ( ">", "", line)
                        
                elif not re.search ("\--|#",line):
                    sequence += line

                if id and sequence:
                    fasta_instance = FASTA(id, sequence)                   
                    fasta_instances.append(fasta_instance)
                    
            return ( fasta_instances )

        
    def fasta_sizer ( fasta_file): ## sizer of fasta sequences, takes dictionary as input ##
        lst = []
        fasta_instances = FASTA.fasta_parser(fasta_file)

        # Define a custom sorting key function
        for fasta_instance in fasta_instances:
            fasta_size = '\t'.join([fasta_instance.id, str(len(fasta_instance.seq))])
            lst.append (fasta_size)
           
        sorted_list = sorted(lst, key=lambda x: x[1], reverse=True) # sort based on sequence length
       
        return sorted_list

       
    def fasta_extract ( fasta_file, gene_list ):
        subset = []
    
        fasta_instances = FASTA.fasta_parser(fasta_file)

        # Define a custom sorting key function
        for fasta_instance in fasta_instances:
            
            if check_isfile ( gene_list ):
                if fasta_instance.id in check_isfile ( gene_list ) :
                    subset.append ( fasta_instance )
            else:
                if re.search ( gene_list, fasta_instance.id ):
                    subset.append ( fasta_instance )
              
        if len(subset) > 0:            
            return subset
        
        else:
            print ( "Gene list is empty or not provided." )

            
    def fasta_remove ( fasta_file, gene_list ):
        subset = []
        fasta_instances = FASTA.fasta_parser(fasta_file)

        # Define a custom sorting key function
        for fasta_instance in fasta_instances:            
            if check_isfile ( gene_list ):
                if fasta_instance.id not in check_isfile ( gene_list ):
                    subset.append ( fasta_instance )
            else:
                if fasta_instance.id != gene_list:
                    subset.append ( fasta_instance )
              
        if len(subset) > 0:            
            return subset
        
        else:
            print ( "Gene list is empty or not provided." )

            

## Implementation ##

def kargs():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Library for efficiently manipulating fasta files.')
    parser.add_argument('-f','--fasta', type=str, help='FASTA file')
    parser.add_argument('-o','--one_line', action="store_true", help='Option to extract FASTA sequences.')
    parser.add_argument('-e','--extract', type=str, help='Option to extract FASTA sequences.')
    parser.add_argument('-r','--remove', type=str, help='Option to remove FASTA sequences.')
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

        inp = re.sub (".fa$|.faa$|.fna$|.fasta$","",args.fasta)
        
        if args.one_line:
            with open(f"{inp}.fasta.1l", "w") as f:
                for out in FASTA.fasta_parser(args.fasta):
                    print ( out, file = f )

        elif args.size:
            with open(f"{inp}.fasta.sizes", "w") as f:
                for out in FASTA.fasta_sizer(args.fasta):
                    print ( out )
                    #print ( out, file = f )

        elif args.extract:
            with open(f"{inp}.{args.extract}.extracted.fasta", "w") as f:
                for out in FASTA.fasta_extract(args.fasta, args.extract):
                    print ( out )

        elif args.remove:
            with open(f"{inp}.{args.remove}.removed.fasta", "w") as f:
                for out in FASTA.fasta_remove(args.fasta, args.remove):
                    print ( out )
    else:
        print("Please provide either --one_line, --extract, --remove, or --size option.")
               

if __name__ == "__main__":
    main()
