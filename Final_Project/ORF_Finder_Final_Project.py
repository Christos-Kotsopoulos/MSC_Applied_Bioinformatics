import re

output_file = open("./Results_standard_ORF.txt", 'w') ## Open a file that will be used later for the storage of the ORFs
reverse_table = str.maketrans({"A" : "U", "G" : "C", "U" : "A", "C" : "G"}) ## A table that will be used in order to create the reverse complement ##
reverse_output_file = open("./Results_reverse_ORF.txt", 'w') ## Open a file that will be used later for the storage of the ORFs


## A good example of a test sequence will be probvided in the below lines. However keep in mind that the ORF this program provides starts with the first AUG codon it meets and ends at the first stop codon it meets.  ##

## Sample sequence:

## CGCTACGTCTTACGCTGGAGCTCTCATGGATCGGTTCGGTAGGGCTCGATCACATCGCTAACCAT  ##
## CGCUTCGUCUUTCGCUGGTGCUCUCTUGGTUCGGUUCGGUTGGGCUCGTUCTCTUCGCUTGCCTU  ##


Translation_table = {
    'UUU': 'F',     'CUU': 'L',     'AUU': 'I',     'GUU': 'V',
    'UUC': 'F',     'CUC': 'L',     'AUC': 'I',     'GUC': 'V',
    'UUA': 'L',     'CUA': 'L',     'AUA': 'I',     'GUA': 'V',
    'UUG': 'L',     'CUG': 'L',     'AUG': 'M',     'GUG': 'V',
    'UCU': 'S',     'CCU': 'P',     'ACU': 'T',     'GCU': 'A',
    'UCC': 'S',     'CCC': 'P',     'ACC': 'T',     'GCC': 'A',
    'UCA': 'S',     'CCA': 'P',     'ACA': 'T',     'GCA': 'A',
    'UCG': 'S',     'CCG': 'P',     'ACG': 'T',     'GCG': 'A',
    'UAU': 'Y',     'CAU': 'H',     'AAU': 'N',     'GAU': 'D',
    'UAC': 'Y',     'CAC': 'H',     'AAC': 'N',     'GAC': 'D',
    'UAA': '--Stop',  'CAA': 'Q',     'AAA': 'K',     'GAA': 'E',
    'UAG': '--Stop',  'CAG': 'Q',     'AAG': 'K',     'GAG': 'E',
    'UGU': 'C',     'CGU': 'R',     'AGU': 'S',     'GGU': 'G',
    'UGC': 'C',     'CGC': 'R',     'AGC': 'S',     'GGC': 'G',
    'UGA': '--Stop',  'CGA': 'R',     'AGA': 'R',     'GGA': 'G',
    'UGG': 'W',     'CGA': 'R',     'AGG': 'R',     'GGG': 'G', 'CGG' : 'R'
}


## Initialising the variables that will be used later ##
num_of_AUG = 0
num_of_UAG = 0
num_of_UGA = 0
num_of_UAA = 0
mscc = 0
reverse_num_of_AUG = 0
reverse_num_of_UAG = 0
reverse_num_of_UGA = 0
reverse_num_of_UAA = 0
reverse_mscc = 0
reverse_counter = 0

## Just some terminal visuals to look cool ##
def print_msg_box(msg, indent=1, width=None, title=None):
    
    lines = msg.split('\n')
    space = " " * indent
    if not width:
        width = max(map(len, lines))
    box = f'╔{"═" * (width + indent * 2)}╗\n'  
    if title:
        box += f'║{space}{title:<{width}}{space}║\n'  
        box += f'║{space}{"-" * len(title):<{width}}{space}║\n'  
    box += ''.join([f'║{space}{line:<{width}}{space}║\n' for line in lines])
    box += f'╚{"═" * (width + indent * 2)}╝'  
    print(box)

info = "This a program used to find simple bacterial ORF's in DNA or RNA sequences. It treats the first AUG codon found from either direction as the start codon and the other AUG's as M residues \n" \
        "It also provided the user with their final ORF, a list where their sequnce's codons are stored and finally the possible protein sequence from the translation of provided sequence \n" \
        "Futhermore it provides the choice for the user to check if the reverse complement of the sequnce her used as input contains an ORF\n" \
        "Finally the program stores the finding in a .txt file in the current directory of the user\n" \
        "To use this program you will need the re library and python3 downloaded  in order to exit you can type -exit- at any time \n\n\n" \
       "To download the re library run:\n" \
       "pip install re\n" \

print_msg_box(msg=info, indent=2, title='Information about this program:')


while True :

    ## Initialising the variables at the start of every loop so they will not mess up the calculations ##
    num_of_AUG = 0
    num_of_UAG = 0
    num_of_UGA = 0  ## Variables for the classic direction of the sequence ##
    num_of_UAA = 0
    mscc = 0

    reverse_num_of_AUG = 0
    reverse_num_of_UAG = 0
    reverse_num_of_UGA = 0 ## Variables for the reverse complement ##
    reverse_num_of_UAA = 0
    reverse_mscc = 0
    reverse_counter = 0
 

    print_msg_box(msg="Write the sequence that you want to be analysed\n", indent=2)
    seq_input = input(" >>").upper()
     
    if seq_input == "EXIT":
        exit()

    seq_input = seq_input.replace("T", "U").upper()

    print("\n")
    ## Check for starting codon ##
    for i in range(0,len(seq_input), 1):
        if seq_input[i:i+3] == "AUG":
            num_of_AUG = num_of_AUG + 1

    if num_of_AUG > 1:
        seq_input_for_aug = seq_input
        seq_input_for_aug = re.split('(AUG)', seq_input_for_aug) ## Getting rid of the sequence infront the AUG codon using regex ##  
        popped_item = seq_input_for_aug.pop(0) 
        sub_seq_str = ''.join(seq_input_for_aug) 
        popped_item_3 = seq_input_for_aug.pop(0) ## Getting rid of the AUG codon to prevent any wrongfully stop codons like A<UAA<AA  ##
        sub_seq_str = ''.join(seq_input_for_aug)
    
        ## Counting for the possible stop codons and identifying them ##
        for i in range(0,len(sub_seq_str), 3):
            if sub_seq_str[i:i+3] == "UAA":
                num_of_UAA = num_of_UAA + 1
            elif  sub_seq_str[i:i+3] == "UAG":
                num_of_UAG = num_of_UAG + 1
            elif  sub_seq_str[i:i+3] == "UGA":
                num_of_UGA = num_of_UGA + 1
    
    mscc = num_of_UAG + num_of_UAA + num_of_UGA
    
    ## Checking for the absence of start codon ##          
    if num_of_AUG == 0:
        print("There is no start codon in your sequence!, please try again")
        

        ## Checking the reverse complement if the user desires in the absence of an AUG codon ##
        for l in range(0,1):


            
            print_msg_box(msg="Do you also want to check for ORF in the reverse complement?(Yes/No)\n", indent=2)
            reverse_complement_choise = input(">>").replace("T", "U").upper()


            if reverse_complement_choise == "YES":

                reverse_complement = seq_input[::-1].translate(reverse_table) ## Reversing the input sequence so that we'll get the reverse complement ##


                ## Check for starting codon in the reverse sequence ##
                for k in range(0,len(seq_input), 1):
                    if reverse_complement[k:k+3] == "AUG":
                        reverse_num_of_AUG = num_of_AUG + 1

                if reverse_num_of_AUG > 1:
                    reverse_seq_input_for_aug = reverse_complement
                    reverse_seq_input_for_aug = re.split('(AUG)', reverse_seq_input_for_aug) ## Getting rid of the sequence infront the AUG codon using regex ##  
                    reverse_popped_item = reverse_seq_input_for_aug.pop(0) 
                    reverse_sub_seq_str = ''.join(reverse_seq_input_for_aug) 
                    reverse_popped_item_3 = reverse_seq_input_for_aug.pop(0) ## Getting rid of the AUG codon to prevent any wrongfully stop codons like A<UAA<AA  ##
                    reverse_sub_seq_str = ''.join(reverse_seq_input_for_aug)

                    ## Counting for the possible stop codons and identifying them in the reverse sequence ##
                    for i in range(0,len(reverse_sub_seq_str), 3):
                        if reverse_sub_seq_str[i:i+3] == "UAA":
                            reverse_num_of_UAA = num_of_UAA + 1
                        elif  reverse_sub_seq_str[i:i+3] == "UAG":
                            reverse_num_of_UAG = num_of_UAG + 1
                        elif  reverse_sub_seq_str[i:i+3] == "UGA":
                            reverse_num_of_UGA = num_of_UGA + 1

                reverse_mscc = reverse_num_of_UAG + reverse_num_of_UAA + reverse_num_of_UGA
               
                ## Checking for the absence of start codon in the reverse sequence ## 
                if reverse_num_of_AUG == 0:
                    print("There is no start codon in your sequence in neither of the complements!, please try again")
                    continue
                
                ## Checking for the absence of stop codons in the reverse sequence ##
                if reverse_mscc == 0 :
                    print("There is no stop codon, in your sequence in neither of the complements! please try again!")
                    continue

                 ## Creating a list that contains three parts the ORF from start codon to the last codon before the possible stop codon, the stop codon and the sequence after the stop codon in the reverse sequence ## 
                test_list = []
                for i in range(0,len(reverse_sub_seq_str), 3):

                    if (reverse_sub_seq_str[i:i+3] != "UAA") and (reverse_sub_seq_str[i:i+3] != "UGA") and (reverse_sub_seq_str[i:i+3] != "UAG"):
                        test_list.append(reverse_sub_seq_str[i:i+3])
                    elif (reverse_sub_seq_str[i:i+3] == "UAA") or (reverse_sub_seq_str[i:i+3] == "UGA") or (reverse_sub_seq_str[i:i+3] == "UAG"):

                        if (reverse_sub_seq_str[i:i+3] == "UAA") :
                            test_list.append("UAA")
                            reverse_sub_seq_str = ''.join(test_list)
                        elif (reverse_sub_seq_str[i:i+3] == "UGA") :
                            test_list.append("UGA")
                            reverse_sub_seq_str = ''.join(test_list)
                        elif (reverse_sub_seq_str[i:i+3] == "UAG"):
                            test_list.append("UAG")
                            reverse_sub_seq_str = ''.join(test_list)
           

                reverse_sub_seq_str = "AUG"+ reverse_sub_seq_str ## Adding the AUD codon that was removed at the start ##

        
                reverse_ORF_calculation_check= len(reverse_sub_seq_str)%3
                if (ORF_calculation_check == 0): ## Checking the triple codon step ##
                    print("Your ORF is",reverse_sub_seq_str, "\n" )
                else :
                    print("No triplets were spotted in your sequence")
                ## Finding, storing and printing the codons in ORF in the reverse sequence ##
                reverse_codons = []
                for z in range(0,len(reverse_sub_seq_str),3):
                    reverse_codons.append(reverse_sub_seq_str[z:z+3])
                print("The codons in your reverse complement are:", reverse_codons)

                ## Translating the codons into a.a residues in the reverse sequence ##
                reverse_protein = []
                for n in range(0,len(reverse_sub_seq_str),3):

                    if reverse_sub_seq_str[n:n+3] in Translation_table:
                        reverse_protein.append(Translation_table.get(reverse_sub_seq_str[n:n+3]))

                reverse_protein_str = ''.join(reverse_protein)
                print("The possible sequence of your protein is :", reverse_protein_str, "\n")

                ## Writing the file that was opened at the start of the program in order to store the input sequence, the ORF and the possible protein sequence of the reverse complement ##
                if reverse_counter == 1 :
                    reverse_output_file.write(f"Starting sequence: {reverse_complement}\n")
                    reverse_output_file.write(f"ORF: {reverse_sub_seq_str}\n")
                    reverse_output_file.write(f"Possible protein sequence:{reverse_protein_str} \n")
                    reverse_output_file.write("........\n" )
                    print("Your results have been written at Results_reverse_ORF.txt at your current directory \n")   

            else :
                print("As you wish!")

        continue

        
    ## Check for the absence of stop codons ##
    if mscc == 0 :
        print("There is no stop codon, in your sequence please try again!")
        ## Checking the reverse complement if the user desires in the absence of a stop codon ##
        for l in range(0,1):
            print_msg_box(msg="Do you also want to check for ORF in the reverse complement?(Yes/No)\n", indent=2) 
            reverse_complement_choise = input(">>").replace("T", "U").upper()


            if reverse_complement_choise == "YES":

                reverse_complement = seq_input[::-1].translate(reverse_table) ## Reversing the input sequence so that we'll get the reverse complement ##

                
                ## Check for starting codon in the reverse sequence ##
                for k in range(0,len(seq_input), 1):
                    if reverse_complement[k:k+3] == "AUG":
                        reverse_num_of_AUG = num_of_AUG + 1

                if reverse_num_of_AUG > 1:

                    reverse_seq_input_for_aug = reverse_complement
                    reverse_seq_input_for_aug = re.split('(AUG)', reverse_seq_input_for_aug) ## Getting rid of the sequence infront the AUG codon using regex ##  
                    reverse_popped_item = reverse_seq_input_for_aug.pop(0) 
                    reverse_sub_seq_str = ''.join(reverse_seq_input_for_aug) 
                    reverse_popped_item_3 = reverse_seq_input_for_aug.pop(0) ## Getting rid of the AUG codon to prevent any wrongfully stop codons like A<UAA<AA  ##
                    reverse_sub_seq_str = ''.join(reverse_seq_input_for_aug)
                
                    ## Counting for the possible stop codons and identifying them in the reverse sequence ##
                    for i in range(0,len(reverse_sub_seq_str), 3):
                        if reverse_sub_seq_str[i:i+3] == "UAA":
                            reverse_num_of_UAA = num_of_UAA + 1
                        elif  reverse_sub_seq_str[i:i+3] == "UAG":
                            reverse_num_of_UAG = num_of_UAG + 1
                        elif  reverse_sub_seq_str[i:i+3] == "UGA":
                            reverse_num_of_UGA = num_of_UGA + 1

                reverse_mscc = reverse_num_of_UAG + reverse_num_of_UAA + reverse_num_of_UGA


                ## Checking for the absence of start codon in the reverse sequence ## 
                if reverse_num_of_AUG == 0:
                    print("There is no start codon in your sequence in neither of the complements!, please try again")
                    continue

                ## Checking for the absence of stop codons in the reverse sequence ##
                if reverse_mscc == 0 :
                    print("There is no stop codon, in your sequence in neither of the complements! please try again!")
                    continue

        

                ## Creating a list that contains three parts the ORF from start codon to the last codon before the possible stop codon, the stop codon and the sequence after the stop codon in the reverse sequence ## 
                test_list = []
                for i in range(0,len(reverse_sub_seq_str), 3):

                    if (reverse_sub_seq_str[i:i+3] != "UAA") and (reverse_sub_seq_str[i:i+3] != "UGA") and (reverse_sub_seq_str[i:i+3] != "UAG"):
                        test_list.append(reverse_sub_seq_str[i:i+3])
                    elif (reverse_sub_seq_str[i:i+3] == "UAA") or (reverse_sub_seq_str[i:i+3] == "UGA") or (reverse_sub_seq_str[i:i+3] == "UAG"):

                        if (reverse_sub_seq_str[i:i+3] == "UAA") :
                            test_list.append("UAA")
                            reverse_sub_seq_str = ''.join(test_list)
                        elif (reverse_sub_seq_str[i:i+3] == "UGA") :
                            test_list.append("UGA")
                            reverse_sub_seq_str = ''.join(test_list)
                        elif (reverse_sub_seq_str[i:i+3] == "UAG"):
                            test_list.append("UAG")
                            reverse_sub_seq_str = ''.join(test_list)
           
               
                reverse_sub_seq_str = "AUG"+ reverse_sub_seq_str ## Adding the AUG codon that was removed at the start ##

    
                reverse_ORF_calculation_check= len(reverse_sub_seq_str)%3
                if (reverse_ORF_calculation_check == 0): ## Checking the triple codon step ##
                    print("Your ORF is",reverse_sub_seq_str, "\n" )
                    reverse_counter = reverse_counter +1
                    
                else :
                    print("No triplets were spotted in your sequence") 

                
                ## Finding, storing and printing the codons in ORF in the reverse sequence ##
                reverse_codons = []
                for z in range(0,len(reverse_sub_seq_str),3):
                    reverse_codons.append(reverse_sub_seq_str[z:z+3])
                print("The codons in your reverse complement are:", reverse_codons)

                ## Translating the codons into a.a residues in the reverse sequence ##
                reverse_protein = []
                for n in range(0,len(reverse_sub_seq_str),3):

                    if reverse_sub_seq_str[n:n+3] in Translation_table:
                        reverse_protein.append(Translation_table.get(reverse_sub_seq_str[n:n+3]))

                reverse_protein_str = ''.join(reverse_protein)
                print("The possible sequence of your protein is :", reverse_protein_str, "\n")

                ## Writing the file that was opened at the start of the program in order to store the input sequence, the ORF and the possible protein sequence of the reverse complement ##
                if reverse_counter == 1 :
                    reverse_output_file.write(f"Starting sequence: {reverse_complement}\n")
                    reverse_output_file.write(f"ORF: {reverse_sub_seq_str}\n")
                    reverse_output_file.write(f"Possible protein sequence:{reverse_protein_str} \n")
                    reverse_output_file.write("........\n" )
                    print("Your results have been written at Results_reverse_ORF.txt at your current directory \n")  

            else :
                print("As you wish!")

        continue
        
            
    ## Creating a list that contains three parts the ORF from start codon to the last codon before the possible stop codon, the stop codon and the sequence after the stop codon ##
    test_list = []
    for i in range(0,len(sub_seq_str), 3):

        if (sub_seq_str[i:i+3] != "UAA") and (sub_seq_str[i:i+3] != "UGA") and (sub_seq_str[i:i+3] != "UAG"):
            test_list.append(sub_seq_str[i:i+3])
        elif (sub_seq_str[i:i+3] == "UAA") or (sub_seq_str[i:i+3] == "UGA") or (sub_seq_str[i:i+3] == "UAG"):

                if (sub_seq_str[i:i+3] == "UAA") :
                    test_list.append("UAA")
                    sub_seq_str = ''.join(test_list)
                elif (sub_seq_str[i:i+3] == "UGA") :
                    test_list.append("UGA")
                    sub_seq_str = ''.join(test_list)
                elif (sub_seq_str[i:i+3] == "UAG"):
                    test_list.append("UAG")
                    sub_seq_str = ''.join(test_list)
           

    sub_seq_str = "AUG"+ sub_seq_str ## Adding the AUG codon that was removed at the start ##

    ORF_calculation_check = len(sub_seq_str)%3
    if (ORF_calculation_check == 0): ## Checking the triple codon step ##
        print("Your ORF is",sub_seq_str, "\n" )
    else :
        print("No triplets were spotted in your sequence")


    ## Finding, storing and printing the codons in ORF ##
    codons = []
    for z in range(0,len(sub_seq_str),3):
        codons.append(sub_seq_str[z:z+3])

    print("Your codons are:", codons)


    ## Translating the codons into a.a residues ##
    protein = []
    for n in range(0,len(sub_seq_str),3):
        
        if sub_seq_str[n:n+3] in Translation_table:
            protein.append(Translation_table.get(sub_seq_str[n:n+3]))
            
    protein_str = ''.join(protein) 
    print("The possible sequence of your protein is :", protein_str, "\n")

    ## Writing the file that was opened at the start of the program in order to store the input sequence, the ORF and the possible protein sequence of the standard complement ##
    output_file.write(f"Starting sequence: {seq_input}\n")
    output_file.write(f"ORF: {sub_seq_str}\n")
    output_file.write(f"Possible protein sequence:{protein_str} \n")
    output_file.write("........\n" )
    print("Your results have been written at Results_standard_ORF.txt at your current directory \n") 


    
    ## Checking the reverse complement if the user desires regardless the absence of an AUG or a stop codon ##
    for l in range(0,1):
            
            print_msg_box(msg="Do you also want to check for ORF in the reverse complement?(Yes/No)\n", indent=2)
            reverse_complement_choise = input(">>").replace("T", "U").upper()

            if reverse_complement_choise == "YES":

                reverse_complement = seq_input[::-1].translate(reverse_table) ## Reversing the input sequence so that we'll get the reverse complement ##
                

                ## Check for starting codon in the reverse sequence ##
                for k in range(0,len(seq_input), 1):
                    if reverse_complement[k:k+3] == "AUG":
                        reverse_num_of_AUG = num_of_AUG + 1
            
                if reverse_num_of_AUG > 1:
                    reverse_seq_input_for_aug = reverse_complement
                    reverse_seq_input_for_aug = re.split('(AUG)', reverse_seq_input_for_aug) ## Getting rid of the sequence infront the AUG codon using regex ##  
                    reverse_popped_item = reverse_seq_input_for_aug.pop(0) 
                    reverse_sub_seq_str = ''.join(reverse_seq_input_for_aug) 
                    reverse_popped_item_3 = reverse_seq_input_for_aug.pop(0) 
                    reverse_sub_seq_str = ''.join(reverse_seq_input_for_aug)

                    ## Counting for the possible stop codons and identifying them in the reverse sequence ##
                    for i in range(0,len(reverse_sub_seq_str), 3):
                        if reverse_sub_seq_str[i:i+3] == "UAA":
                            reverse_num_of_UAA = num_of_UAA + 1
                        elif  reverse_sub_seq_str[i:i+3] == "UAG":
                            reverse_num_of_UAG = num_of_UAG + 1
                        elif  reverse_sub_seq_str[i:i+3] == "UGA":
                            reverse_num_of_UGA = num_of_UGA + 1

                reverse_mscc = reverse_num_of_UAG + reverse_num_of_UAA + reverse_num_of_UGA
               
                
                if reverse_num_of_AUG == 0:
                    print("There is no start codon in your sequence in neither of the complements!, please try again") ## Checking for the absence of start codon in the reverse sequence ##
                    continue
                if reverse_mscc == 0 :
                    print("There is no stop codon, in your sequence in neither of the complements! please try again!") ## Checking for the absence of stop codons in the reverse sequence ##
                    continue

                 ## Creating a list that contains three parts the ORF from start codon to the last codon before the possible stop codon, the stop codon and the sequence after the stop codon in the reverse sequence ## 
                test_list = []
                for i in range(0,len(reverse_sub_seq_str), 3):

                    if (reverse_sub_seq_str[i:i+3] != "UAA") and (reverse_sub_seq_str[i:i+3] != "UGA") and (reverse_sub_seq_str[i:i+3] != "UAG"):
                        test_list.append(reverse_sub_seq_str[i:i+3])
                    elif (reverse_sub_seq_str[i:i+3] == "UAA") or (reverse_sub_seq_str[i:i+3] == "UGA") or (reverse_sub_seq_str[i:i+3] == "UAG"):

                        if (reverse_sub_seq_str[i:i+3] == "UAA") :
                            test_list.append("UAA")
                            reverse_sub_seq_str = ''.join(test_list)
                        elif (reverse_sub_seq_str[i:i+3] == "UGA") :
                            test_list.append("UGA")
                            reverse_sub_seq_str = ''.join(test_list)
                        elif (reverse_sub_seq_str[i:i+3] == "UAG"):
                            test_list.append("UAG")
                            reverse_sub_seq_str = ''.join(test_list)
        
                reverse_sub_seq_str = "AUG"+ reverse_sub_seq_str
        
                reverse_ORF_calculation_check= len(reverse_sub_seq_str)%3
                if (reverse_ORF_calculation_check == 0): ## Checking the triple codon step ##
                    print("Your ORF is",reverse_sub_seq_str, "\n" )
                    reverse_counter = reverse_counter +1 
                else :
                    print("No triplets were spotted in your sequence")

                
                ## Finding, storing and printing the codons in ORF in the reverse sequence ##
                reverse_codons = []
                for z in range(0,len(reverse_sub_seq_str),3):
                    reverse_codons.append(reverse_sub_seq_str[z:z+3])
                print("The codons in your reverse complement are:", reverse_codons)

                ## Translating the codons into a.a residues in the reverse sequence ##
                reverse_protein = []
                for n in range(0,len(reverse_sub_seq_str),3):

                    if reverse_sub_seq_str[n:n+3] in Translation_table:
                        reverse_protein.append(Translation_table.get(reverse_sub_seq_str[n:n+3]))
                    else :
                        reverse_protein.append('unknown')

                reverse_protein_str = ''.join(reverse_protein)
                print("The possible sequence of your protein is :", reverse_protein_str, "\n")

                ## Writing the file that was opened at the start of the program in order to store the input sequence, the ORF and the possible protein sequence of the reverse complement ##
                if reverse_counter == 1 :
                    reverse_output_file.write(f"Starting sequence: {reverse_complement}\n")
                    reverse_output_file.write(f"ORF: {reverse_sub_seq_str}\n")
                    reverse_output_file.write(f"Possible protein sequence:{reverse_protein_str} \n")
                    reverse_output_file.write("........\n" ) 
                    print("Your results have been written at Results_reverse_ORF.txt at your current directory \n") 

            else :
                print("As you wish!")
    

                
                


    

