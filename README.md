# FastaLib
Custom Python library for parsing and manipulating fasta files.

Includes functions for:
  1. Converting multi-line fasta to one-line fasta
  2. Extracting sequences of user-specified genes
  3. Removing sequences of user-specified genes
  4. Printing length of fasta sequences
  5. Translate CDS to peptide (PEP) sequences
  6. Extracting or removing user-selected sequence(s) and subsequence(s)
  7. Replace sequence IDs
  8. QC of sequence IDs -- checks for duplicate IDs and for IDs longer than 50 characters (not allowed for making blast databases)
