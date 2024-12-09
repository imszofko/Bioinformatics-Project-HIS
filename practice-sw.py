#Input the sequences, will add function to parse the input and make sure that only ATCGU are in the sequence
import os 
import re

seq1 = str(input('Input sequence 1: ')).upper()
seq2 = str(input('Input sequence 2: ')).upper()
"""
#Similar scoring to NW but it doesnt allow for negative scores, so the neg is set to 0
match = 1                       #When there is a match
mismatch = -1                   #I am not sure what to put as some sources say different things
gap = -1
maxScore = 0
#Want to add a feature that has the different PAM and BLOSUM paramters stored and use that as the substitution matrix instead
##Scores taken from Wikipedia 
"""   
##Function to make a matrix of zeros
'''
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
''' 
#initMat = []
'''
def MatZero(rows, cols):
    tempMat = []
    for i in range(rows + 1):
        tempMat.append([])
        for j in range(cols + 1):
            tempMat[-1].append(0)   #-1 index is always the last index in an array, it will add 0s to the the last index of the original adding arrays and also the 0s in.
    return tempMat
'''
def SmithWaterman(seq1, seq2, mismatch = -1, gap = -1, match = 1, maxScore = 0):
    #Create the matrix
    s1 = len(seq1)                  #Sequence one is the vertical sequence to the left
    s2 = len(seq2)                  #Sequence two is teh horizontal sequence on top
    
    tempMat = []
    for i in range(s1 + 1):
        tempMat.append([])
        for j in range(s2 + 1):
            tempMat[-1].append(0)   #-1 index is always the last index in an array, it will add 0s to the the last index of the original adding arrays and also the 0s in.
    print(tempMat)
    #scoreMatrix = MatZero(s1, s2)
    #Filling in the matrix based off of the algorithm
    '''
    for i in range(1, s1 + 1):
        for j in range(1, s2 + 1):
            if seq1[i - 1] == seq2[j - 1]:
                scoreMatrix[i][j] = max(scoreMatrix[i][j - 1] +  gap, scoreMatrix[i - 1][j] + gap, scoreMatrix[i - 1][j - 1] + match, 0)
            else:
                scoreMatrix[i][j] = max(scoreMatrix[i][j - 1] + gap, scoreMatrix[i - 1][j] + gap, scoreMatrix[i - 1][j - 1] + mismatch, 0)
    '''
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
    #Reverse the strings
    alignA = alignA[::-1]
    alignB = alignB[::-1]
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
    
    ##Calculating alignment scores of the local alignment
    alignScore = 0
    for i in range(len(matchString)):
        if matchString[i] == '|':
            alignScore += 1
        elif (matchString[i] == '*' or matchString[i] == " "):
            alignScore += -1
    
    return alignA, alignB, matchString, alignScore
#Printing out the results
output1, output2, matchString, alignScore = SmithWaterman(seq1, seq2, mismatch = -1, gap = -1, match = 1, maxScore = 0)
#print(output1 + '\n' + output2)
print(output1 + '\n' + matchString + '\n' + output2)
print("Alignment Score: ", alignScore)