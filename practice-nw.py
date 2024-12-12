import requests, sys
from Bio import SeqIO

#Establishing the API to ENSMBL
'''rest.ensembl.org/sequence/id/ENSG00000139618'''
ensemblID = input(str("Enter ENSEMBL Sequence ID: "))
server = "https://rest.ensembl.org"
extension = f"/sequence/id/{ensemblID}"

ensemblResponse = requests.get(server + extension, headers = {"Content-Type" : "text/x-fasta"})

if ensemblResponse.ok:
    # Get FASTA content as a string
    fasta_content = ensemblResponse.text
    print(f"Fetched FASTA:\n{fasta_content}")

    # Extract the sequence (ignoring the header)
    fasta_lines = fasta_content.strip().splitlines()
    header = fasta_lines[0]  # First line is the header
    sequence = ''.join(fasta_lines[1:])  # Remaining lines are the sequence
    print("Sequence Extracted.")

    # Save to a file
    with open("Seq_ENSEMBL.fasta", "w") as file:
        file.write(fasta_content)  # Save as is
    print("Sequence written to file.")
if not ensemblResponse.ok:
  ensemblResponse.raise_for_status()
  sys.exit()


#Going to write the fasta sequence to a file since I think it will 
#be easier on the machine and it is easier for me to code right now 
# (something to scale up to is to not have to save the sequence to a file)
#Input the sequences, will add function to parse the input and make sure that only ATCGU are in the sequence
'''
seq1 = str(input('Input sequence 1: ')).upper()
seq2 = str(input('Input sequence 2: ')).upper()
'''
##Now I will add the code to import two files and read the fasta file and store that seq in a variable
'''currently using two test files and not the BRCA2 and variant file'''

for seq1_record in SeqIO.parse("Seq_ENSEMBL.fasta", "fasta"):
    seq1 = seq1_record.seq
    print("Sequence one:" , seq1[0:10])

for seq2_record in SeqIO.parse("Variant700202-cDNA.fasta", "fasta"):
    seq2 = seq2_record.seq
    print("Sequence two:", seq2[0:10])
##I thought I had to use .read not .parse because Biopython manual says it is good for one entry

mismatch = -1
gap = -1
match = 1

##Function to make a matrix of zeros
def zeros(rows, cols):
    #Define an empty list
    matZero = []
    #Set up rows in the matrix
    for x in range(rows):
        #For each row add an empty list
        #this appends an empty array for each x in rows
        matZero.append([])
        #Set up the columns in each row
        for y in range(cols):
            #add zero to each y in cols in each row
            #-1 index is always the last index in an array, it will add 0s to the the last index of the original adding arrays and also the 0s in.
            matZero[-1].append(0)
    #return the matrix of zeros
    #print(matZero)
    return matZero

#zeros(3,5) #testing zeros function

#Return the score between any two bases in alignment
def matchScore(match1, match2):
    if match1 == match2:
        return match
    elif match1 == '-' or match2 == '-':
        return gap
    else:
        return mismatch
    
##Function that fills out the matrix of scores
def NeedlemanWunsch(seq1, seq2, gap = -1, match = 1):
    #Length of two sequence
    s1 = len(seq1) #n
    s2 = len(seq2) #m

    #Generate the matrix of zeros to stores the scores
    #+1 to create the 0s col and rows
    
    scoreMatrix = zeros(s2 + 1, s1 + 1) 
    ##Had an issue where I couldnt align two sequences of different lenghts and I found the issue to be here
    #I originally had it s1+1,s2+1 and it is supposed to be the opposite way. Idk why but when I made this change it fixed the issue.

    #Calculate score table
    #First col and first row of zeros
    #Columns
    for i in range(0, s2 + 1):
        scoreMatrix[i][0] = gap * i

    #Rows
    for j in range(0, s1 + 1):
        scoreMatrix[0][j] = gap * j

    for i in range(1, s2 + 1):
        for j in range(1, s1 + 1):

            #Calculating the score by checking the top, the left, and the diagonal squares
            
            match = scoreMatrix[i - 1][j - 1] + matchScore(seq1[j-1], seq2[i-1])                  #Diagonal one from [i][j] position at that time
            insertGap = scoreMatrix[i - 1][j] + gap                                               #Up one of [i][j] position at that time
            misDelete = scoreMatrix[i][j - 1] + gap                                               #To the left of [i][j] position at that time

            #Now going to recording the maximum score from the three possibilities calculated
            scoreMatrix[i][j] = max(match, insertGap, misDelete)
    #return scoreMatrix

    #Traceback from Wikipedia pseudocode
    #Creating a variable to store each alignment
    alignA = ""
    alignB = ""

    ##Storing the length of both sequences
    i = len(seq2) 
    j = len(seq1)

    ##While loop to trace where we are in the matrix during the traceback
    ##The traceback starts in the lower right square of the matrix and moves backwards throughout the matrix
    ##Starts in bottom right corner and compares the value with the three possible sources to see which it comes from
    #If seq1 and seq2 are aligned they match and it is +1 and the same with the rest of the stuff essentially

    while i > 0 and j > 0:
        #checking that it is a match, if it is then it appends to alignA/B and then jump to the diagnol value after
        if scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1] + matchScore(seq1[j-1], seq2[i-1]):                           #Diagonal step
            alignA += seq1[j - 1]
            alignB += seq2[i - 1]
            i -= 1 #move forward
            j -= 1
        
        #Checking of the index [i][j] and one above are equal with and to update i and j to correspond to that cell
        elif scoreMatrix[i][j] == scoreMatrix[i - 1][j] + gap:                                                          #Step above the [i][j]
            alignA += seq1[j - 1]
            alignB += '-'
            i -= 1
        
        elif scoreMatrix[i][j] == scoreMatrix[i][j - 1] + gap:                                                          #Step to the left of [i][j]
            alignA += '-'
            alignB += seq2[i - 1]
            j -= 1

    #WIll append to alignA or B depending on the result of j or i (seq 1 or seq2) and finish tracing up to the top left cell
    ##Appending occurs also before this line (ln 106 to 123)
    while j > 0:
        alignA += seq1[j-1]
        alignB += '-'
        j -= 1
    
    while i > 0:
        alignA += '-'
        alignB += seq2[i-1]
        i -= 1

    #Reversing the sequences as they were flipped from the traceback
    alignA = alignA[::-1]
    alignB = alignB[::-1]

    matchString = ""
    for i in range (len(alignA)):
        if alignA[i] == alignB[i]:
            matchString += '|'
        elif alignA[i] != alignB[i]:
            if (alignA[i] == '-' or alignB == '-'):
                matchString += " "
            else:
                matchString += "*"
    alignScore = 0
    for i in range(len(matchString)):
        if matchString[i] == '|':
            alignScore += 1
        elif (matchString[i] == '*' or matchString[i] == " "):
            alignScore += -1
    
    return alignA, alignB, matchString, alignScore
output1, output2,matchString, alignScore = NeedlemanWunsch(seq1, seq2, gap = -1, match = 1)  
##Storing the outputs of seq1 and seq2 and other variables and printing the results
print(output1 + '\n' + matchString + '\n' + output2)
print("Alignment Score: ", alignScore)


#To write the sequences to a .txt file to look like EMBOSS
def split_sequence(seq, width):     #sequence is either output1 or output2 or the match string
    return [seq[i:i+width] for i in range(0, len(seq), width)]

#This will split the sequences and alignment markers into chunks with the width of the line being 60
lineWidth = 60
chunks1 = split_sequence(output1, lineWidth)
chunks2 = split_sequence(output2, lineWidth)
match_chunks = split_sequence(matchString, lineWidth)

#Writes the alignment to the text file
output_file = "Variance700202-Alignment.txt"
with open(output_file, "w") as file:
    for c1, match, c2 in zip(chunks1, match_chunks, chunks2):
        file.write(f"{c1}\n")
        file.write(f"{match}\n")
        file.write(f"{c2}\n\n")  #Line for readability

print(f"Alignment written to {output_file}")


