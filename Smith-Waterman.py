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

for seq1_record in SeqIO.parse("Seq_ENSEMBL.fasta", "fasta"):
    seq1 = seq1_record.seq
    print("Sequence one:" , seq1[0:10])

for seq2_record in SeqIO.parse("Variant700202-cDNA.fasta", "fasta"):
    seq2 = seq2_record.seq
    print("Sequence two:", seq2[0:10])
##I thought I had to use .read not .parse because Biopython manual says it is good for one entry


''' 
seq1 = str(input('Input sequence 1: ')).upper()
seq2 = str(input('Input sequence 2: ')).upper()

#Similar scoring to NW but it doesnt allow for negative scores, so the neg is set to 0
match = 1                       #When there is a match
mismatch = -1                   #I am not sure what to put as some sources say different things
gap = -1
maxScore = 0
#Want to add a feature that has the different PAM and BLOSUM paramters stored and use that as the substitution matrix instead
##Scores taken from Wikipedia 
'''

##Function to make a matrix of zeros
def zeros(rows, cols):
    #Define an empty list
    matZero = []
    #Set up rows in the matrix
    for x in range(rows + 1):
        #For each row add an empty list
        #this appends an empty array for each x in rows
        matZero.append([])
        #Set up the columns in each row
        for y in range(cols + 1):
            #add zero to each y in cols in each row
            #-1 index is always the last index in an array, it will add 0s to the the last index of the original adding arrays and also the 0s in.
            matZero[-1].append(0)
    #return the matrix of zeros
    #print(matZero)
    return matZero

def MatZero(rows, cols):
    tempMat = []
    for i in range(rows + 1):
        tempMat.append([])
        for j in range(cols + 1):
            tempMat[-1].append(0)   #-1 index is always the last index in an array, it will add 0s to the the last index of the original adding arrays and also the 0s in.
    return tempMat

def SmithWaterman(seq1, seq2, mismatch = -1, gap = -1, match = 1, maxScore = 0):
    #Create the matrix
    s1 = len(seq1)                  #Sequence one is the vertical sequence to the left
    s2 = len(seq2)                  #Sequence two is teh horizontal sequence on top
    
    tempMat = []
    for i in range(s1 + 1):
        tempMat.append([])
        for j in range(s2 + 1):
            tempMat[-1].append(0)   #-1 index is always the last index in an array, it will add 0s to the the last index of the original adding arrays and also the 0s in.
    #print(tempMat)
    #scoreMatrix = MatZero(s1, s2)
    #Filling in the matrix based off of the algorithm
    print("Matrix made.")
    '''
    for i in range(1, s1 + 1):
        for j in range(1, s2 + 1):
            if seq1[i - 1] == seq2[j - 1]:
                tempMat[i][j] = max(tempMat[i][j - 1] +  gap, tempMat[i - 1][j] + gap, tempMat[i - 1][j - 1] + match, 0)
            else:
                tempMat[i][j] = max(tempMat[i][j - 1] + gap, tempMat[i - 1][j] + gap, tempMat[i - 1][j - 1] + mismatch, 0)
    '''
    print("Starting to fill matrix.")
    for i in range(1, s1 + 1):
        for j in range(1, s2 + 1):
            if seq1[i - 1] == seq2[j - 1]:
                scoreMatch = tempMat[i - 1][j - 1] + match
            else:
                scoreMatch = tempMat[i - 1][j - 1] + mismatch
            
            scoreInsert = tempMat[i][j - 1] + gap
            scoreDelete = tempMat[i - 1][j] + gap

            tempMat[i][j] = max(0, scoreMatch, scoreInsert, scoreDelete)
            if tempMat[i][j] >= maxScore:
                maxScore = tempMat[i][j]
                mismatch, gap = (i, j)
    print("Matrix filled.")
    '''
    This way of filling in the sub matrix was not as good and it was not as accurate as the one above.
    When comparing the same sequence with the same smaller segement it gave a wacky result
    ''' 
           
    #Traceback starts at the element with the highest score
    #Finding the element with the highest score
    alignA = ""
    alignB = ""
    maximum = 0
    for row in range(1, s1 + 1):
        for column in range(1, s2 + 1):
            if maximum < tempMat[row][column]:
                maximum = tempMat[row][column]           #max value is stored
                i = row                                     #Storing the row 
                j = column                                  #Storing column
    print("Starting traceback.")
    #Backtracing
    '''
    Based on the source of each score recursively until 0 is encountered
    Segments that have the highest similarity score based on the given scoring system is generated in this processed
    to obtain the second best local alignment, apply the traceback process starting at the second highest score outside the trace of the best alignment
    Very important is that no negative score is assigned in the scorying system which enables local alignment
    '''
    while tempMat[i][j] != 0:
        if seq1[i - 1] == seq2[j - 1]:
            alignA += seq1[i - 1]
            alignB += seq2[j - 1]
            i -= 1
            j -= 1
        #This is for if the sequences dont match
        elif seq1[i - 1] != seq2[j - 1]:
            tempList = [tempMat[i - 1][j - 1], tempMat[i - 1][j], tempMat[i][j - 1]]
            if max(tempList) == tempList[0]:
                alignA += seq1[i - 1]
                alignB += seq2[j - 1]
                i -= 1
                j -= 1
            #If the maximum value is the 1st indexed position, topvalue
            elif max(tempList) == tempList[1]:
                alignA += seq1[i - 1]
                alignB += '-'
                i -= 1
            #If the max value is the second indexed position, leftvalue
            elif max(tempList) == tempList[-1]:
                alignA += '-'
                alignB += seq2[j - 1]
                j -= 1
        else:
            print('Error. Exit.')
            i = 0
            j = 0
    print("Backtrace complete")
    #Reverse the strings
    alignA = alignA[::-1]
    alignB = alignB[::-1]
    print("Sequences reversed")
    #Storing match  and such  scores and symbols
    matchString = ""
    for i in range (len(alignA)):
        if alignA[i] == alignB[i]:
            matchString += '|'
        elif alignA[i] != alignB[i]:
            if (alignA[i] == '-' or alignB == '-'):
                matchString += " "
            else:
                matchString += "*"
    print("Match scores noted.")
    ##Calculating alignment scores of the local alignment
    alignScore = 0
    for i in range(len(matchString)):
        if matchString[i] == '|':
            alignScore += 1
        elif (matchString[i] == '*' or matchString[i] == " "):
            alignScore += -1
    print("Alignment score calculated.")
    return alignA, alignB, matchString, alignScore
#Printing out the results
output1, output2, matchString, alignScore = SmithWaterman(seq1, seq2, mismatch = -1, gap = -1, match = 1, maxScore = 0)
#print(output1 + '\n' + output2)
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