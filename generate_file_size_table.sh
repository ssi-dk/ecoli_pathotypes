ls -l fasta_files/ > l_output_of_fasta_files.txt;
awk '{print $5,$9}' l_output_of_fasta_files.txt > fasta_file_sizes.txt
