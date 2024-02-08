

#!/usr/bin python3

import os
from natsort import natsorted
from xopen import xopen
import argparse
import re
from collections import defaultdict


def isfile ( gene_list, delimiter = '\t' ):
    try:
        with open ( gene_list, "r" ) as file:
            num_columns = len(file.readline().split(delimiter)) ## Assuming that each line has the same number of columns            
            lines = file.readlines()
            
            instances = []
            
            for line in lines:
                instances.append(line)

    except FileNotFoundError:    
        instances = gene_list

    return instances



class FASTA:

    __slots__ = ['id','seq']
   
    def __init__ (self, header = "", sequence = ""):
        self.id = str(header)
        self.seq = str(sequence)
        

    def __str__ ( self ) :
        return f">{self.id}\n{self.seq}"
    

    def fasta_parser(fasta_file, sprint = False) -> list:

        '''

        Parse a FASTA file and return a list of FASTA instances.

        Parameters:
        - fasta_file (str): Path to the FASTA file.
        - sprint (bool): If True, print the parsed sequences.

        Returns:
        - list: List of FASTA instances.

        '''      
        
        fasta_instances = []

        with xopen(fasta_file, 'r') as fasta:
            lines = fasta.readlines()
            seqid = ""
            sequence = ""
            entries = defaultdict()
            encountered_ids = set()

            for line in lines:
                line = line.strip('\n')  # Remove newline character

                if line.startswith('>'):
                    seqid = re.sub("\s+", "_", line[1:])
                    seqid = re.sub("\r", "\n", line[1:])
                    sequence = ""

                    if seqid in encountered_ids:
                        print(f"Error: Duplicate seqid found - '{seqid}'.")

                    else:
                        encountered_ids.add(seqid)

                    if len(seqid) > 50: # Check if seqid length is greater than 50
                        print(f"Error: Seqid '{seqid}' has length greater than 50.")
                            
                elif not re.search("--|#", line):  # Ignore lines containing "--" or "#"
                    sequence += re.sub("\*$|\.$", "", re.sub("\W", "", line))

                    if seqid and sequence:
                        entries[seqid] = sequence

        for seqids, sequences in entries.items():
            fasta_instance = FASTA(seqids, sequences)
            fasta_instances.append(fasta_instance)
            
        return fasta_instances
    

    def output_one_by_one ( fasta_file ):

        '''
        Option to output each sequence one-by-one in a new FASTA file.

        '''
        
        fasta_instances = FASTA.fasta_parser ( fasta_file )

        for fasta_instance in fasta_instances:

            with open ( fasta_instance.id, "w" ) as out:
                out.write ( str(fasta_instance) )

    
    def fasta_sizer ( fasta_file) -> list:
        
        '''
        Size and sort FASTA sequences based on sequence length.

        Parameters:
        - fasta_file (str): Path to the FASTA file.

        Returns:
        - list: Sorted list of sequence sizes.
        
        '''
        
        lst = []
        fasta_instances = FASTA.fasta_parser(fasta_file)

        # Define a custom sorting key function
        for fasta_instance in fasta_instances:
            fasta_size = '\t'.join([fasta_instance.id, str(len(fasta_instance.seq))])
            lst.append (fasta_size)
            
        sorted_list = natsorted(lst, key=lambda x: x[1], reverse=False) # sort based on sequence length
        
        return sorted_list

       
    def replace_names ( fasta_file, gene_list ) -> list:

        '''
        Replaces original IDs of fasta sequences with new associated IDs.
        
        Parameters:
        - fasta_file (str): Path to the FASTA file.
        - gene_list (str): Path to the list of user-specified gene IDs.

        Returns:
        - fasta: Fasta file with replace IDs for specified sequences.
        
        '''
       
        fasta_instances = FASTA.fasta_parser ( fasta_file )

        new_fasta_instances = []
        
        for fasta_instance in fasta_instances:
            for g in isfile( gene_list ):
                if re.search ( fasta_instance.id , g ):
                    fasta_instance.id = g

            new_fasta_instances.append(fasta_instance)

        return new_fasta_instances

    
    def translate(fasta_file) -> list:

        '''
        Translate nucleotide sequences to protein sequences.

        Parameters:
        - fasta_file (str): Path to the FASTA CDS file.

        Returns:
        - list: List of FASTA instances with translated protein sequences.
        
        '''
        
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

            
    def extract_subsequences (fasta_file, gene_list, from_file, start_position="", end_position="", extract = True, ncbi_tsa_submission = False) -> list:

        '''
        Extracts or removes sequence(s) and subsequence(s) from a fasta file

        Parameters:
          - fasta_file (str): Path to the FASTA file.
          - gene_list (str): Name of file with provided sequence(s).
          - from_file(bool): True if the sequences' names, start and end positions are provided in a file. False if sequence name, start and position are provided as a string.
          - start_position(int): Start position in the sequence. Default: empty_string
          - end_position(int): End position in the sequence. Default: empty_string
          - extract(bool): Extract or remove subsequence(s). Default: extract=True
 
        Returns:
          - list: List of sequences based on options to extract or remove entire sequence(s) and subsequence(s).

        '''
        
        fasta_instances = FASTA.fasta_parser(fasta_file)
        matches = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
        intact_sequences = []
        intact_sequences_r = []
        subsequences = []
        start_pos = int()
        end_pos = int()
        
        for fasta_instance in fasta_instances: ## Find matches between the fasta file and the provided gene list, which is either a file or a string, and add them into a dictionary
            if os.path.isfile(gene_list):
                for gene in isfile(gene_list):
                    if len(gene.split('\t')) == 3:
                        if re.search(gene.split('\t')[0], fasta_instance.id):
                            matches[fasta_instance.id][gene.split('\t')[1]][gene.split('\t')[2]] = fasta_instance.seq
                    elif len(gene.split('\t')) == 1:
                        if re.search(gene, fasta_instance.id):
                            matches[fasta_instance.id] = fasta_instance.seq
                    else:
                        print (f'Please the format of your input {gene_list} file.')

            elif isinstance(gene_list, str): ## Provided gene list is a string, NOT a file

                if start_position == "": ## If no start and/or end positions are provided, return the entire sequence
                    start_position = 1
                    
                if end_position == "":
                    end_position = len(sequence)

                if re.search(gene_list, fasta_instance.id):
                    matches[fasta_instance.id][start_position][end_position] = fasta_instance.seq

            else:
                print ('Provided file is neither a FILE nor a STRING. Please check your input.')
                return fasta_instances

                    
        if len(matches) > 0: ## If matches are found, extract or remove entire sequence(s) or subsequence(s) based on user-provided specifications
            
            intact_sequences = [fasta_instance for fasta_instance in fasta_instances if fasta_instance.id not in matches.keys()] ## Creates a list with all sequences NOT matching the gene list ONCE. Isolates the sequences we want to remove from the file.
                    
            for matching_seq, positions in matches.items(): ## Create dictionaries to save the fasta ID, start, end positions and sequence of matching sequences.
                if isinstance(positions, dict):                   
                    for start_positions, end_positions in positions.items():
                        start_position = int(start_positions)

                        for end_positions, sequences in end_positions.items():
                            end_position = int(end_positions)
                            sequence = str(sequences)                          
        
                if start_position and end_position and sequence: 
                    
                    if 1 <= start_position <= len(sequence) and 1 <= end_position <= len(sequence):
                        
                        if extract:  ## Options to extract sequence(s) and subsequence(s)
                            subsequence = sequence[start_position:end_position]
                            subsequences.append( FASTA(matching_seq, subsequence) )
                                         
                        else: ## Options to remove sequence(s) and subsequence(s)
                            if start_position > 1: 
                                start_pos = start_position - 1
                
                            if end_position < len(sequence):
                                end_pos = end_position + 1
                                
                            if not start_pos == 1 and not end_pos == len(sequence): ## If start and end position are in the middle of the sequence, collapse the subsequence and merge the flanking regions.
                                subsequence = sequence[:start_pos] + sequence[end_pos:]
                                intact_sequences_r.append( FASTA(matching_seq, subsequence) )

                            elif not start_pos == 1 and end_pos == len(sequence): ## If start is in the middle and end position is last position in sequence, keep only up to the start
                                subsequence = sequence[:start_pos]
                                intact_sequences_r.append( FASTA(matching_seq, subsequence) ) 

                            elif start_pos == 1 and not end_pos == len(sequence): ## If start position is 1 and end position is in the middle, keep sequence from start to the provided end position
                                subsequence = sequence[1:end_pos]
                                intact_sequences_r.append( FASTA(matching_seq, subsequence) )
                                    
                    else:
                        print(f"Requested range {start_position} - {end_position} not present in {fasta_instance.id}")
                else:
                    print (f"Start and end positions not provided for {matching_seq}.")
        else:
            print ('No matching sequences found in fasta file.')
            return fasta_instances
            
            
        ## Size filtering here, if ncbi_tsa_submission is enabled to avoid doing multiple times in the loop
        extracted_sequences = []
        kept_sequences = []              
        
        if subsequences: ## Extract option enabled and subsequence(s) were extracted
            for fasta_sequence in subsequences:
                if ncbi_tsa_submission == True:
                    if len(fasta_sequence.seq) >= 200:
                        extracted_sequences.append(fasta_sequence)

                else:
                    extracted_sequences.append(fasta_sequence)    
                        
            return subsequences


        elif len(intact_sequences_r) > 0 and not subsequences: ## Remove option enabled and subsequence(s) were removed           
            intact_sequences_r = intact_sequences_r + intact_sequences

            for fasta_sequence in intact_sequences_r:
                if ncbi_tsa_submission == True:
                    if len(fasta_sequence.seq) >= 200:
                        kept_sequences.append(fasta_sequence)
                else:
                    kept_sequences.append(fasta_sequence)
    
            return kept_sequences

        
        elif intact_sequences and not subsequences and not intact_sequences_r: ## Remove option enabled for entire sequence(s)
            for fasta_sequence in intact_sequences:
                if ncbi_tsa_submission == True:
                    if len(fasta_sequence.seq) >= 200:
                        kept_sequences.append(fasta_sequence)
                else:
                    kept_sequences.append(fasta_sequence)
            
            return kept_sequences



    def check_position (fasta_file, gene_list, position_number, length = 0) -> str:

        '''
        Returns the sequence from a fasta file, based on position number.

        Parameters:
          - fasta_file (str): Path to the FASTA file.
          - gene_list (str): Name of the sequence.
          - position_number (int): Position number in the sequence.
          - length (int): Optional length parameter.

        Returns:
          - str: Subsequence based on the provided position and length.

        '''
        
        fasta_instances = FASTA.fasta_parser(fasta_file)

        sequence = None
        sequences = []

        for fasta_instance in fasta_instances:
            if re.search(gene_list, fasta_instance.id):
                sequence = fasta_instance
                break

        if sequence:
            position_number -= 1 ## counting starts by default at zero
            
            if not position_number > len(sequence.seq) :

                if length > 0:
                    end_position = position_number + length
                    sequences.append(sequence.seq[position_number:end_position])
                    return ''.join(sequences)

                elif length < 0:
                    start_position = position_number - length
                    sequences.append(sequence.seq[start_position:position_number])
                    return ''.join(sequences)
              
                else:
                    return sequence.seq[position_number + length]

            else:
                print ( f"Requested position number not present in {sequence.id}" )
                return None
        else:
            print ( f"Requested sequence {gene_list} does not exist in {fasta_file}" )
            return None
        

        
## Implementation as a main script ##

def parse_arguments():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Library for efficiently manipulating fasta files.')

    ## options to parse and convert multi-line fasta to one-line fasta
    parser.add_argument('-f','--fasta', type=str, help='FASTA file')
    parser.add_argument('-o','--one_line', action="store_true", help='Option to extract FASTA sequences.')

    ## option to compute size of fasta sequences
    parser.add_argument('-s','--size', action="store_true", help='Option to print sizes of FASTA sequences.')

    ## options to extract or remove genes
    parser.add_argument('-obo','--one_by_one', action="store_true", help='Print all fasta sequences in individual files one-by-one.')
    parser.add_argument('-e','--extract', action="store_true", help='Option to extract FASTA sequences.')
    parser.add_argument('-r','--remove', action="store_true", help='Option to remove FASTA sequences.')
    parser.add_argument('-g','--gene_list', type=str, help='User-provided list of sequence IDs to extract, remove or replace with new in FASTA file.')
    parser.add_argument('-st','--start_position', type=str, required = False, help='Start position to extract FASTA sequences.')
    parser.add_argument('-end','--end_position', type=str, required = False, help='End position to extract FASTA sequences.')
    parser.add_argument('-ff','--from_file', action="store_true", help='Option to get start and end positions from file.')
    parser.add_argument('-ncbi','--ncbi_tsa_submission', action="store_true", help='Performcs checks for NCBI TSA submissions.')

    ## replace sequence names; original gene ids should be associated with new names
    parser.add_argument('-nn', '--new_names', action="store_true", help='Option to replace names/IDs of FASTA sequences.')

    ## option to convert CDS to PEP 
    parser.add_argument('-tr','--cds2pep', action="store_true", help='Option to print sizes of FASTA sequences.')

    ## get sequence based on position number
    parser.add_argument('-pos','--position', action="store_true", help='Option to return requested position number of FASTA sequence.')
    parser.add_argument('-gene_list','--sequence_name', type=str, help='Name of sequence to return requested position from a number of FASTA sequence.')
    parser.add_argument('-num','--number', type=int, help='Position number to requested sequence from a FASTA sequence.')
    parser.add_argument('-len','--length', type=int, help='Length to add to position number to requested sequence from a FASTA sequence. e.g. 280 + length')

    args = parser.parse_args()

    if not any(vars(args).values()):
        parser.print_help()
        print("Error: No arguments provided.")

    return parser.parse_args()


def main():
    
    parser = argparse.ArgumentParser(description='Library for efficiently manipulating fasta files.')
    args = parse_arguments()

    threads = []

    if args.fasta:

        inp = re.sub (".aa$|.fa$|.faa$|.fna$|.fasta$|.1l$","",args.fasta)

        if args.one_line:
            with open(f"{inp}.fasta.1l", "w") as f:
                for out in FASTA.fasta_parser(args.fasta):
                    print ( out, file = f )

        elif args.one_by_one:
            FASTA.output_one_by_one(args.fasta)

        elif args.size:
            with open(f"{inp}.fasta.sizes", "w") as f:
                for out in FASTA.fasta_sizer(args.fasta):
                    print ( out, file = f )

        elif args.extract:
            if args.gene_list:
                if args.from_file:
                    if args.ncbi_tsa_submission:
                        with open(f"{inp}.extracted.from_{args.gene_list}.ncbi_tsa_submission.fasta", "w") as f:
                            for out in FASTA.extract_subsequences(fasta_file = args.fasta, gene_list = args.gene_list, ncbi_tsa_submission = True, from_file = True, extract = True):
                                print ( out, file = f )
                    else:
                        with open(f"{inp}.extracted.from_{args.gene_list}.fasta", "w") as f:
                            for out in FASTA.extract_subsequences(fasta_file = args.fasta, gene_list = args.gene_list, from_file = True, extract = True):
                                print ( out, file = f )
                else:
                    with open(f"{inp}.extracted.fasta", "w") as f:
                        if args.ncbi_tsa_submission:
                            with open(f"{inp}.extracted.from_{args.gene_list}.ncbi_tsa_submission.fasta", "w") as f:
                                for out in FASTA.extract_subsequences(fasta_file = args.fasta, gene_list = args.gene_list, start_position = args.start_position, end_position = args.end_position, ncbi_tsa_submission = True, from_file = False, extract = True):
                                    print ( out, file = f )
                        else:
                            with open(f"{inp}.extracted.from_{args.gene_list}.fasta", "w") as f:
                                for out in FASTA.extract_subsequences(fasta_file = args.fasta, gene_list = args.gene_list, start_position = args.start_position, end_position = args.end_position, from_file = False, extract = True):
                                    print ( out, file = f )
        elif args.remove:
            if args.gene_list:
                if args.from_file:
                    if args.ncbi_tsa_submission:
                        with open(f"{inp}.removed.from_{args.gene_list}.ncbi_tsa_submission.fasta", "w") as f:
                            for out in FASTA.extract_subsequences(fasta_file = args.fasta, gene_list = args.gene_list, ncbi_tsa_submission = True, from_file = True, extract = False):
                                print ( out, file = f )

                    else:
                        with open(f"{inp}.removed.from_{args.gene_list}.fasta", "w") as f:
                            for out in FASTA.extract_subsequences(fasta_file = args.fasta, gene_list = args.gene_list, from_file = True, extract = False):
                                print ( out, file = f )
                else:
                    with open(f"{inp}.removed.fasta", "w") as f:
                        if args.ncbi_tsa_submission:
                            with open(f"{inp}.removed.{args.gene_list}.ncbi_tsa_submission.fasta", "w") as f:
                                for out in FASTA.extract_subsequences(fasta_file = args.fasta, gene_list = args.gene_list, start_position = args.start_position, end_position = args.end_position, ncbi_tsa_submission = True, from_file = False, extract = False):
                                    print ( out, file = f )
                        else:
                            with open(f"{inp}.removed.{args.gene_list}.fasta", "w") as f:
                                for out in FASTA.extract_subsequences(fasta_file = args.fasta, gene_list = args.gene_list, start_position = args.start_position, end_position = args.end_position, from_file = False, extract = False):
                                    print ( out, file = f )
            else:
                print ( "Please provide a gene list.")
                
        elif args.new_names:
            with open(f"{inp}.new_names.fasta", "w") as f:
                for out in FASTA.replace_names(args.fasta, args.gene_list):
                    print ( out, file = f )
                    
        elif args.cds2pep:
            with open(f"{inp}.pep.fasta", "w") as f:
                for out in FASTA.translate(args.fasta):
                    print ( out, file = f )

        elif args.position:
            if args.sequence_name:
                if args.number >=0:
                    if FASTA.check_position ( args.fasta, args.sequence_name, args.number, args.length if args.length else 0 ):
                        print ( "Requested site number", args.number, "in", args.sequence_name, "is", args.end==" " )
                        print ( FASTA.check_position ( args.fasta, args.sequence_name, args.number, args.length if args.length else 0 ) )
                else:
                    print ("Please provide a position number to return.")
            else:
                print ( "Please provide the name of the sequence to search for returning the position number.")        
                
        else:
            print("Please select a potential FASTA operation.")

    else:
        print ( "Please provide the input fasta file.")
        

if __name__ == "__main__":
    main()
