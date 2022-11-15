# BINF_Assignment2.py
# Jesse Wolf 0830233

# Importing the argv and regular expressions modules.
from sys import argv
import re 

# Setting fasta_input to our first argument and enzymes_input to the second argument that we receive from the command line.
fasta_input=argv[1]
enzyme_input=argv[2]

# Creating a function so we can put line breaks in our output.
def line_break():
    print ("-"*63)

# Checking to see if our fasta file is in the right format, if not, we exit and ask the user to input a fasta file.
if ".fa" in fasta_input:
    print ("File is of the right type, proceeding with program!")
    line_break()
else:
    exit ("The sequence file provided cannot be read, please input a fasta file.")

# Creating a function when we get an error for both our fasta and enzyme input.
def fasta_absent_error():
        exit("File " + fasta_input + " is not available, please input a valid file.")
def enzyme_absent_error():
        exit("File " + enzyme_input + " is not available, please input a valid file.")

# Using a try and except block for each of our fasta and enzyme input to check if our input fasta and enzyme files are available in the current directory.
try:
    open (fasta_input, "r")
except FileNotFoundError as e:
    fasta_absent_error()

try:
    open (enzyme_input, "r") 
except FileNotFoundError as e:
    enzyme_absent_error()

# Taking our input files and opening them.
fasta_file = open (fasta_input)
enzyme_file = open (enzyme_input)

# Creating a list of lines for each of the fasta file and enzymes .txt file.
fasta_lines_list = fasta_file.readlines()
enzyme_lines_list = enzyme_file.readlines()

# Closing both our fasta and enzyme file.
fasta_file.close()
enzyme_file.close()

# Checking to see if our enzmes file has enzymes present from which to cut with.
if  len(enzyme_lines_list) <1:
    exit ("No enzymes could be read in the file provided.")

# Getting our sequence name and actual sequence data and creating new variables for both.
fasta_sequence_name = fasta_lines_list[0].lstrip(">")
fasta_sequence = fasta_lines_list[1]

# Printing the fasta and enzyme files that are being used, separated by a line break, and then the name of the sequence and how long it is.
print ("Restriction enzyme analysis of sequence from file " + fasta_input + ".")
print ("Cutting with enzymes found in file " + enzyme_input + ".") 
line_break()
print ("Sequence name: " + fasta_sequence_name + "Sequence is " + str(len(fasta_sequence)) + " bases long.")

# Creating our dictionary of enzymes and their cut sites.
enzyme_dict={}
for i in enzyme_lines_list:
    # Stripping the newline characters from each element in our enzyme list.
    i=i.rstrip("\n")
    # Splitting the values from each element by a semicolon.
    enzyme_split=i.split(";")
    # Adding keys and corresponding values to our dictionary called enzyme_dict.
    enzyme_dict[enzyme_split[0]]=enzyme_split[1]

# Creating a for loop for each enzyme in our enzyme dictionary to cut our sequence data.
for enzyme in enzyme_dict:
    # Splitting each enzyme at the cut site denoted by the % symbol.
    enzyme_split = enzyme_dict[enzyme].split("%")
    # Creating a new variable for our enzymes without the % symbol, using the .join command.
    enzyme_whole = "".join(enzyme_split)
    # Creating a new variable for our enzymes without the % symbol and introducing a line break at the cut site, using the .join command.
    enzyme_digested = "\n".join(enzyme_split)
    # Printing a line break between enzymes.
    line_break()
    # Using re.findall to find cut sites within our fasta sequence using each enzyme.
    if re.findall(enzyme_whole, fasta_sequence):
        # Creating a variable that stores the number of cut sites for each enzyme within our fasta sequence.
        num_cut_sites = len(re.findall(enzyme_whole, fasta_sequence))
        # Creating a variable representing the fragments that each enzyme cuts the sequence into using re.finditer.
        fragments = re.finditer(enzyme_whole, fasta_sequence)
        # Printing a statement with the number of cut sites for each enzyme and where they cut and replacing the % sign with ^ to reflect assignment output requirements.
        print ("There are " + str(num_cut_sites) + " cut sites for " + enzyme + ", cutting at " + str(enzyme_dict[enzyme].replace("%", "^")) + ".")
        # Printing a statement with the number of fragments resulting from each enzyme. We need to add 1 to the cut sites variable as for each cut, there are n+1 fragments produced.
        print ("There are " + str(num_cut_sites+1) + " fragments:" +"\n")
        # Creating a variable for each fragment, replacing each instance within our sequence where an enzyme cut site is found with enzyme_digested, which has a line break at the cut site. This breaks the sequence into n parts, where n is the number of fragments (or the number of cuts +1). We use the split function to break the fragments into a list where the separator is a line break.
        fragments = (str(fasta_sequence).replace(enzyme_whole, enzyme_digested)).split("\n")
        # Creating a variable for the start of sequence positions which we iterate over below.
        seq_start_position = 1
        # Setting up a for loop to iterate over each fragment produced by our restriction enzymes.
        for frag in fragments:
            # For each fragment, calculate the length and print it.
                fragment_lengths = len(frag)
                print ("length: " + str(fragment_lengths))
                # Using a step length of 60, for each value in the range from 0 to the length of a given fragment, employ the below nested loop.
                for x in range (0,len(frag), 60):
                    # Using a step length of 10, for each value from the above for loop, in the range from 0 to the 60th base pair within a given fragment, employ the following.
                    for i in range (0,len([frag[x:x+60]]),10):
                        # Join any fragments up to 60 base pairs in length separated by a newline.
                        frag_break60 = "\n".join([frag[x:x+60]])
                        # Take the fragments separated into 60 base pair lines and re-join them at 10 base pair intervals separated by a space, using regex.
                        frag_final = " ".join(re.findall(r".{1,10}",frag_break60))
                        # Print the start of a given 60 base pair chunk (an integer) and the corresponding 60 base pair fragment.
                        print (str(seq_start_position) + "\t" + frag_final)
                        # Add the length of a fragment to our existing seq_start_position variable and assign the new value back to the same variable to keep a running count between fragments.
                        seq_start_position += len(frag[x:x+60])
        # Printing an extra line break to match the assignment format.
        print ("\n")
    # If there are no matches for a given enzyme in our sequence, we print that there are no sites found for that enzyme and an extra line break.
    else: 
            print ("There are no sites for " + enzyme + ".")
            print ("\n")