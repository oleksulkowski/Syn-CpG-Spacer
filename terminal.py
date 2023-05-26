import objects
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys, os
from datetime import datetime


def main():
    def load_file():
        while True:
            gene_file_name = input("Enter the name of the Fasta file with the sequence to be recoded: ")
            if chr(27) in gene_file_name:
                sys.exit()
            elif gene_file_name.strip() == "":
                print("Please enter a file name.")
            else:
                if not gene_file_name.endswith(".fasta"):
                    gene_file_name += ".fasta"
                try:
                    with open(gene_file_name):
                        pass
                except FileNotFoundError:
                    print("File not found. Please enter a valid file name.")
                else:
                    break
        # Read the original gene file
        for seq_record in SeqIO.parse(gene_file_name, "fasta"):
            original_sequence = seq_record.seq
            return original_sequence, gene_file_name

    def prompt_for_integer(prompt_message):
        while True:
            input_str = input(prompt_message).strip()
            if chr(27) in input_str:
                sys.exit()
            if input_str.lower() == "o":
                return 'o'
            try:
                input_int = int(input_str)
                if input_int < 0:
                    raise ValueError
            except ValueError:
                print("Please enter a non-negative integer.")
            else:
                return input_int

    def prompt_for_specific_numbers(prompt_message, target_nums):
        while True:
            input_str = input(prompt_message).strip()
            if chr(27) in input_str:
                sys.exit()
            try:
                input_int = int(input_str)
                if input_int not in target_nums:
                    raise ValueError
            except ValueError:
                print(f"Please enter one of the following numbers: {target_nums}.")
            else:
                return input_int

    def get_input_variables():
        minimum_CpG_gap = None
        desired_CpG_gap = None
        
        packaging_signal_length_beginning = prompt_for_integer("Enter the length of the packaging signal to be protected in the beginning of the sequence: ")
        packaging_signal_length_end = prompt_for_integer("Enter the length of the packaging signal to be protected in the end of the sequence: ")
        print("")
        print("If you know the minimum CpG gap you would like to apply, enter 1")
        print("If you have a desired average CpG gap and would like the script to calculate an optimal minimum CpG gap, enter 2")
        print("If you would like to create a sequence with the most consistent CpG spacing (with lowest SD of gap), enter 3")
        gap_method = prompt_for_specific_numbers("Choose: ", [1, 2, 3])

        if gap_method == 1:
            minimum_CpG_gap = prompt_for_integer("Enter the minimum length of a gap between an existing CpG and one to be added: ")
        elif gap_method == 2:
            desired_CpG_gap = prompt_for_integer("Enter the desired average gap: ")
        

        print("")
        return packaging_signal_length_beginning, packaging_signal_length_end, gap_method, minimum_CpG_gap, desired_CpG_gap


    def print_troubleshoot_details():
        gene.read_codons()
        gene.print_changeable_CpG()
    
    def perform_A_mutations(gene):
        print(f"The new average CpG gap is: {gene.calculate_average_gap()}")
        while True:
            ask_A_rich = input("Do you want to make the sequence A-rich? (Y/N): ")
            if chr(27) in ask_A_rich:
                sys.exit()
            if ask_A_rich.lower() == "y":
                gene.mutate_A_rich()
                break
            elif ask_A_rich.lower() == "n":
                break
    
    def collect_statistics(gene):
        original_CpG_count = gene.count_CpGs(gene.original_sequence)
        final_CpG_count = gene.count_CpGs(gene.new_sequence)
        original_average_gap = gene.calculate_average_gap("original")
        final_average_gap = gene.calculate_average_gap()
        print('')
        A_abundance_change = gene.calculate_A_abundance_change()
        print('')
        return original_CpG_count, final_CpG_count, original_average_gap, final_average_gap, A_abundance_change

    def get_output_id(gene_file_name):
        output_id = input("Enter an ID for the output file: ")
        if chr(27) in output_id:
            sys.exit()

        if gene_file_name.endswith(".fasta"):
            gene_file_name = gene_file_name[:-6]
        return output_id

    # Saves the gene in a new file
    def save_new_gene(sequence, gene_file_name, id, final_CpG_count, original_CpG_count, final_average_gap, original_average_gap, A_abundance_change, packaging_signal_length_beginning, packaging_signal_length_end, minimum_CpG_gap):
        my_object = Seq(str(sequence))

        description = f"{gene_file_name}. There are {final_CpG_count} CpG's here vs {original_CpG_count} CpG's in the original sequence. The new average CpG distance is {final_average_gap} (down from {original_average_gap}). A-nucleotides abundance increased by {A_abundance_change}%. Protected original {packaging_signal_length_beginning} nucleotides and final {packaging_signal_length_end}. Minimum new CpG gap was set as {minimum_CpG_gap}"

        record = SeqRecord(my_object,
                            id=id,
                            name=gene_file_name,
                            description=description,
                            annotations={"molecule_type": "DNA"})    
        
        # Get current date as string
        current_date = datetime.now().strftime('%d-%m-%Y')

        # Create a new folder named with the current date if it does not exist
        if not os.path.exists(f'outputs/{current_date}'):
            os.makedirs(f'outputs/{current_date}')

        return SeqIO.write(record, f"outputs/{current_date}/{gene_file_name}-recoded-{id}.fasta", "fasta")

    print("At any prompt, you can exit the program by pressing Esc followed by Enter")
    print("")
    original_sequence, gene_file_name = load_file()
    packaging_signal_length_beginning, packaging_signal_length_end, gap_method, minimum_CpG_gap, desired_CpG_gap = get_input_variables()

    gene = objects.Gene(original_sequence, packaging_signal_length_beginning, packaging_signal_length_end, gap_method, minimum_CpG_gap, desired_CpG_gap)

    # Uncomment if you want to print the original gene sequence in terminal and to print mutable CpG-enriching codons    
    #print_troubleshoot_details()

    perform_A_mutations(gene)
    gene.check_synonymity()
    gene.check_packaging_signals()
    original_CpG_count, final_CpG_count, original_average_gap, final_average_gap, A_abundance_change = collect_statistics(gene)

    save_new_gene(gene.new_sequence, gene_file_name, get_output_id(gene_file_name), final_CpG_count, original_CpG_count, final_average_gap, original_average_gap, A_abundance_change, packaging_signal_length_beginning, packaging_signal_length_end, minimum_CpG_gap)

if __name__ == "__main__":
    main()