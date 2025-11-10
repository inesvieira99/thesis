mport sys
import getopt
import random
import os

# Step 1: Load FASTA file and extract sequence
def load_fasta(file_path):
    if not os.path.exists(file_path):
        print(f"Error: The file {file_path} does not exist.")
        sys.exit(1)
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        header = lines[0].strip()  # Fasta header
        sequence = ''.join(line.strip() for line in lines[1:])
    return header, sequence

# Step 2: Apply random modifications
def modify_sequence(sequence, modification_percentage):
    sequence = list(sequence)  # Convert to list for easier manipulation
    num_modifications = int(len(sequence) * (modification_percentage / 100))
    
    for _ in range(num_modifications):
        index = random.randint(0, len(sequence) - 1)
        current_base = sequence[index]
        new_base = random.choice([b for b in 'AGTC' if b != current_base])
        sequence[index] = new_base
    
    return ''.join(sequence)

# Step 3: Apply fixed modification
def apply_fixed_modification(sequence, fixed_percentage):
    sequence = list(sequence)  # Convert to list for easier manipulation
    num_modifications = int(len(sequence) * (fixed_percentage / 100))
    
    for _ in range(num_modifications):
        index = random.randint(0, len(sequence) - 1)
        current_base = sequence[index]
        new_base = random.choice([b for b in 'AGTC' if b != current_base])
        sequence[index] = new_base
    
    return ''.join(sequence)

# Step 4: Save the modified sequence as a new FASTA file
def save_fasta(header, sequence, output_path):
    try:
        with open(output_path, 'w') as file:
            file.write(f"{header}\n")
            for i in range(0, len(sequence), 60):  # Write in lines of 60 characters
                file.write(f"{sequence[i:i+60]}\n")
    except Exception as e:
        print(f"Error: Could not save the file {output_path}. {e}")
        sys.exit(1)

# Main function to modify scaffold
def modify_scaffold(input_fasta, user_modification_percentage, fixed_modification, random_output_fasta, final_output_fasta):
    # Load the input scaffold
    header, sequence = load_fasta(input_fasta)
    
    # Add a random variation of ± 0.5% to the modification percentage
    random_variation = random.uniform(-0.5, 0.5)
    modification_percentage = round(user_modification_percentage + random_variation, 1)  # Round to 1 decimal place
    
    print(f"Applied random modification percentage: {modification_percentage}%")
    
    # Apply random modification
    random_modified_sequence = modify_sequence(sequence, modification_percentage)
    
    # Save the random modified sequence to a new FASTA file
    save_fasta(header, random_modified_sequence, random_output_fasta)
    print(f"Randomly modified scaffold saved to {random_output_fasta}")
    
    # Apply the fixed modification to the random modified sequence
    final_modified_sequence = apply_fixed_modification(random_modified_sequence, fixed_modification)
    
    # Save the final modified sequence to a second FASTA file
    save_fasta(header, final_modified_sequence, final_output_fasta)
    print(f"Final modified scaffold (with additional {fixed_modification}% modifications) saved to {final_output_fasta}")

# Command-line handling using getopt
def main(argv):
    input_fasta = ''
    user_modification_percentage = 0
    fixed_modification = 0
    random_output_fasta = ''
    final_output_fasta = ''

    try:
        opts, args = getopt.getopt(argv, "i:d:p:o:f:", ["input=", "divergence=", "polimorfism=", "outdiv=", "outpoli="])
    except getopt.GetoptError as err:
        print(err)
        print("Usage: python3 poli.py -input <input_fasta> -divergence <user_modification_percentage> -polimorfism <fixed_modification> -outdiv <random_output_fasta>
        sys.exit(2)

    # Parse options and arguments
    for opt, arg in opts:
        if opt in ("-i", "--input"):
            input_fasta = arg
        elif opt in ("-d", "--divergence"):
            try:
                user_modification_percentage = float(arg)  # Parse divergence as float
            except ValueError:
                print(f"Error: Invalid divergence percentage value '{arg}'. Must be a number.")
                sys.exit(2)
        elif opt in ("-p", "--polimorfism"):
            try:
                fixed_modification = float(arg.replace(",", "."))  # Handle decimal commas
            except ValueError:
                print(f"Error: Invalid polimorfism value '{arg}'. Must be a number.")
                sys.exit(2)
        elif opt in ("-o", "--outdiv"):
                elif opt in ("-o", "--outdiv"):
            random_output_fasta = arg
        elif opt in ("-f", "--outpoli"):
            final_output_fasta = arg

    # Validate that all required arguments are provided
    if not input_fasta:
        print("Error: Missing input FASTA file (-input).")
        sys.exit(2)
    if user_modification_percentage <= 0:
        print("Error: Missing or invalid divergence percentage (-divergence). Must be a positive number.")
        sys.exit(2)
    if fixed_modification <= 0:
        print("Error: Missing or invalid fixed modification percentage (-polimorfism). Must be a positive number.")
        sys.exit(2)
    if not random_output_fasta:
        print("Error: Missing output file for divergence (-outdiv).")
        sys.exit(2)
    if not final_output_fasta:
        print("Error: Missing output file for polimorfism (-outpoli).")
        sys.exit(2)

    # Call the scaffold modification function
    modify_scaffold(input_fasta, user_modification_percentage, fixed_modification, random_output_fasta, final_output_fasta)

# Entry point
if __name__ == "__main__":
    main(sys.argv[1:])