#Input the sequences, will add function to parse the input and make sure that only ATCGU are in the sequence
seq1 = str(input('Input sequence 1: '))
seq1 = seq1.upper()
seq2 = str(input('Input sequence 2: '))
seq2 = seq2.upper()


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

mismatch = -1
gap = -1
match = 1

#Return the score between any two bases in alignment
def matchScore(match1, match2):
    if match1 == match2:
        return match
    elif match1 == '-' or match2 == '-':
        return gap
    else:
        return mismatch
    
##Function that fills out the matrix of scores
def NeedlemanWunsch(seq1, seq2):
    #Length of two sequence
    s1 = len(seq1) #n
    s2 = len(seq2) #m

    #Generate the matrix of zeros to stores the scores
    #+1 to create the 0s col and rows
    
    scoreMatrix = zeros(s1 + 1, s2 + 1)

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
            
            match = scoreMatrix[i - 1][j - 1] + matchScore(seq1[j-1], seq2[i-1])
            insertGap = scoreMatrix[i - 1][j] + gap
            misDelete = scoreMatrix[i][j - 1] + gap

            #Now going to recording the maximum score from the three possibilities calculated
            scoreMatrix[i][j] = max(match, insertGap, misDelete)
    #return scoreMatrix

    #Traceback from Wikipedia pseudocode
    #Creating a variable to store each alignment
    alignA = ""
    alignB = ""

    ##Storing the length of both sequences
    i = len(seq2) #m
    j = len(seq1) #n

    ##While loop to trace where we are in the matrix during the traceback
    ##The traceback starts in the lower right square of the matrix and moves backwards throughout the matrix
    ##Starts in bottom right corner and compares the value with the three possible sources to see which it comes from
    #If seq1 and seq2 are aligned they match and it is +1 and the same with the rest of the stuff essentially

    while i > 0 and j > 0:
        #checking that it is a match, if it is then it appends to alignA/B and then jump to the diagnol value after
        if scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1] + matchScore(seq1[j-1], seq2[i-1]):
            alignA += seq1[j - 1]
            alignB += seq2[i - 1]
            i -= 1 #move forward
            j -= 1
        
        elif scoreMatrix[i][j] == scoreMatrix[i - 1][j] + gap:
            alignA += seq1[j - 1]
            alignB += '-'
            i -= 1
        
        else: #scoreMatrix[i][j] == scoreMatrix[i][j - 1] + gap:
            alignA += '-'
            alignB += seq2[i - 1]
            j -= 1
        
    while j > 0:
        alignA += seq1[j-1]
        alignB += '-'
        j -= 1
    
    while i > 0:
        alignA += '-'
        alignB += seq2[i-1]
        i -= 1

    alignA = alignA[::-1]
    alignB = alignB[::-1]

    return alignA, alignB

output1, output2 = NeedlemanWunsch(seq1, seq2)
print(output1 + "\n" + output2)
