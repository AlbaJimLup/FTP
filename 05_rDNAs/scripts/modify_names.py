
def modify_fasta(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    with open(output_file, 'w') as f:
        for line in lines:
            if line.startswith('>'):
                parts = line.strip().split(':')
		#change seq header name
                f.write('>YAO_'+parts[0][-2:]+":"+parts[1]+'\n')
            else:
		#keep seq as it is 
                f.write(line)

# Provide input and output file paths
input_file = 'YAO_45S_hp_unique.fa'
output_file = 'YAO_45S_hp_unique_names.fa'

# Modify FASTA file and save to output file
modify_fasta(input_file, output_file)

