# Word report accuracy - a script that scores word report data
#
# Written in Python 2 by Patrick McClure in collaboration with Heidi S. Økland 
# at the MRC Cognition and Brain Sciences Unit, April/May 2016. 
# Modified September 2017 by HSØ (nothing in the main code, though!).
#
# This is a script that will compare two sets of written text and compute the 
# number and percent of matching words between the two sets. 
#
# It was originally written to compare the number of matching words in two 
# sentences - one as uttered by a participant (subject) in an experiment 
# ("transcription"), and the other being the sentence they tried to repeat 
# back ("target"). Participants said "don't know" if they weren't able to 
# repeat any words back, hence there is a portion of the code below where the 
# score is set to 0 whenever the transcription says "don't know".
# 
# PLEASE NOTE: The script is written in such a way that it will score matching
# words as correct regardless of word order. A potential problem is that this 
# approach will give the same score to a sentence with jumbled words as to a 
# sentence with correctly ordered words. Given that subjects usually tend to
# a) not produce randomly ordered sentences and b) report words in more or less
# the correct order, this scoring algorithm is the best we could come up with 
# in order to automate the scoring.
#
# To run this script, you'll need:
# 1) A data folder (e.g. ..\Data\ in which all the subject-specific 
# transcriptions are stored in subfolders (..\Data\01\, ..\Data\02\, etc). 
# Put this script in \Data\! This makes it easier to run it.
# 2) A text file with the target sentences. It's probably best to copy them 
# from an excel sheet into e.g. Notepad such that there is one sentence per 
# line. The path for this file is hard-coded below in "targetPath".
# 3) Text files with transcriptions that are named in a consistent way (again,
# you'd probably want to copy from excel) such as 'Transcription_SUBJ01.txt',
# 'Transcription_SUBJ02.txt', etc. This means we can easily loop over subjects
# by defining the default name ("transcriptionPath" below) and simply add the
# subject number to it, plus '.txt'.
# 4) A text file with some sort of "sentence IDs", e.g. S001, S002, .., with 
# one ID per line (again, copy from excel). The path for this file is hard-
# coded below in "idPath".
#
# The script outputs to the command line (using the function print) as well as
# to a text file. The text file contains the following info:
# Sentence ID, number of words in target, number of words in transcription,
# number of matched words (i.e. word report score), the proportion of matched 
# words correct, the transcription, and the target sentences.


from __future__ import division
import numpy as np

import os
from os import path

# Define the working directory (which is the same as the folder that this script lives in)
directory = path.dirname(path.realpath(__file__))

# Define the filename for the text file which includes the original text to 
# which we want to compare (henceforth "target")
sentFile = 'TargetSentences_AVSpeechMEG.txt'

# Define the sentence ID text file filename 
idFile = 'SentenceIDs_AVSpeechMEG.txt'

# Define the default name for the text file which contains the transcription
transFile = 'Transcription_SUBJ' # e.g. 'Transcription_SUBJ' so that we can just add '01', '02', etc. as we loop over subjects

# Define the default name for the output text file
outFile = 'WordReportScores_SUBJ'

# List of the subject numbers to score (also names of the folders)
subjects = ['sub-01' 'sub-02'] # e.g. ['01', '02' ..]

# Define file suffix
suffix = '.txt'

# Loop over subjects/data sets
for s in subjects:
    
    # Define output and transcription files (with full paths) for this subject
    outFile_s = os.path.join(directory, s, outFile + s + suffix)
    transFile_s = os.path.join(directory, s, transFile + s + suffix)
                           
    # Open the output file to save the word report accuracy scores
    file = open(outFile_s,'w')

    # Write output file headers
    file.write('Sentence ID\tNumber of words in target sentence\tNumber of words in transcription\tNumber of correct words\tWords correct\tTranscription\tTarget\n')

    # Open transcription, target and sentence ID files
    # (each is a text file where each piece of text is separated by a \n, i.e. new line)
    with open(transFile_s) as transcriptionfile, open(sentFile) as targetfile, open(idFile) as idfile: 
        
         # Loop over sentences
        for x, y, z in zip(transcriptionfile, targetfile, idfile):
            # make transcription lower case and print it            
            x = x.strip().lower() 
            x_print = x
            z = z.strip().split() 
            print(z)
            print(x)
   
            # If transcription says 'don't know': 
            if x == 'don\'t know':
                # set score to 0
                score = 0
                # set number of words in transcript to 0
                x = '' 
                # make target sentence lower-case and print it
                y = y.strip().lower()
                y_print = y
                print(y)
                # split target sentence into words
                y = y.split()
                
            # Otherwise count the number of matching words in the transcription
            else:
                # split transcription into words
                x = x.split()
                # make target sentence lower-case and print it
                y = y.strip().lower()
                y_print = y
                print(y)
                # split target sentence into words
                y = y.split()
                
                # Find the matching words in the transcription (x) and the target (y).
                # ifMatched is a binary array where 0s correspond to unmatched transcription words and 1s correspond to matched transcription words
                ifMatched = np.zeros((len(x),), dtype=np.int)

                # For each word in the target
                for yword in y:
                    i = -1
                    # For each word in the transciption
                    for xword in x:
                        i = i + 1
                        
                        # If the target and transcription word match AND the transcription word has not been matched already (i.e. hasn't been "seen" before)
                        if yword == xword and ifMatched[i] == 0:
                            # ..then label the transcription word as matched (= count it as a correct word)
                            ifMatched[i] = 1
                            break

                # Calculate the number of matched words for this sentence
                score = sum(ifMatched)
            
            # Print output
            print(z[0] + '\t' + str(len(y)) + '\t' + str(len(x)) + '\t' + str(score) + '\t' + str(float(score)/float(len(y))*100)+'%' )
            file.write(z[0]  + '\t' + str(len(y)) + '\t' + str(len(x)) + '\t' + str(score) + '\t' + str(float(score)/float(len(y))*100)+'%' + '\t' + x_print + '\t' + y_print + '\n') # output to text file
    # Close the text file
    file.close()#file.write(z[0] + '\t' + str(len(y)) + '\t' + str(score) + '\t' + str(round(score/len(y)*100))+'%' + '\n')# + '\t' + x_print + '\t' + y_print + '\n') # output to command line