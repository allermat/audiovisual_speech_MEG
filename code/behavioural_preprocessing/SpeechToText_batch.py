# -*- coding: utf-8 -*-
"""
@author: Connor Quinn

This script takes a list of folders (participants). For each folder it tries
to recognise the words in the .wav files and outputs the text results.

USE:
You should save this script in a folder, e.g. '/speech_data/heidi_script.py'
You should also have subfolders for each of your participants:
/speech_data/heidi_script.py
/speech_data/01
/speech_data/02
/speech_data/03

Probably easiest if the folder (participant) names are numerical, but no
problem if not. Within each participant folder you should have the .wav files:
/speech_data/01/file1.wav
/speech_data/01/file2.wav
/speech_data/01/file3.wav

The only thing you should have to change in this script is the list of
folders to look in (folder_list). The elements in this list must match the names you have
given to the participants' folders; the '01' in folder_list corresponds to the
folder /speech_data/01

"""


# Here we specify which libraries to import.
import speech_recognition as sr
import os
from os import path

folder_list = ["sub-01" "sub-02"] # change this to match your folder structure

'''
Before running any code first set up the functions that we will use:
You shouldn't need to change any of the functions.
'''

def recog_wav(wav_file):
    '''
    Reads in one wav file at a time and passes out the resulting text.
    This is the critical piece of code
    '''
    r = sr.Recognizer()
    with sr.WavFile(str(wav_file)) as source:
        audio = r.record(source) # read the entire WAV file
    try:
        # for testing purposes, we're just using the default API key
        # to use another API key, use `r.recognize_google(audio, key="GOOGLE_SPEECH_RECOGNITION_API_KEY")`
        # instead of `r.recognize_google(audio)`
        return r.recognize_google(audio)
        # r.recognize_ibm(audio, username=IBM_USERNAME, password=IBM_PASSWORD))
    except sr.UnknownValueError:
        print "Google Speech Recognition could not understand audio"
    except sr.RequestError as e:
        print "Could not request results from Google Speech Recognition service; {0}".format(e)

def recog_files(list, file_out_path):
    '''
    Read in a list of .wav files and output the text to a txt.file
    This function calls the 'recog_files' function for each .wav files in the list.
    '''
#    for wav_file in list:
#        out_text = recog_wav(wav_file)
#        with open(str(file_out_path), 'a') as out_doc:
#            out_doc.write('\t'.join([wav_file, out_text]).encode('utf-8'))
#            #out_doc.write('\t'.join(map(str, [wav_file, out_text]))) # this line sometimes produces an error!
#            out_doc.write('\n')
#            out_doc.close()
            
    for wav_file in list:
        try:        
            out_text = recog_wav(wav_file)
            with open(str(file_out_path), 'a') as out_doc:
                out_doc.write('\t'.join([wav_file, out_text]).encode('utf-8'))
                out_doc.write('\n')
                out_doc.close()
        # sometimes the encoding of the Google Speech Recognition output fails, 
        # so we add an exception to keep the script running:      
        except:
            with open(str(file_out_path), 'a') as out_doc:
                out_doc.write('\t'.join([wav_file, 'TRANSCRIPTION FAILED']).encode('utf-8'))
                out_doc.write('\n')
                out_doc.close() 

def set_out_file(file_out_path):
    '''
    Create an output file and set the column headings
    '''
    f_out = open(file_out_path, 'a')
    f_out.write('\t'.join(['Filename', 'Automatic transcription',]))
    f_out.write('\n')
    f_out.close()


#All of the action happens below here.
# Gets the path to whereever you are running this script
DIR = path.dirname(path.realpath(__file__))

for subj in folder_list: # for each participant
    SUBJ_DIR = os.path.join(DIR,  subj) # find their folder
    print SUBJ_DIR    
    # set the path for the output file
    OUT_PATH = os.path.join(DIR, 'SpeechToText')
    if not os.path.exists(OUT_PATH):
        os.makedirs(OUT_PATH)
    file_out_path = os.path.join(OUT_PATH, 'SUBJ{}_SpeechToText.txt'.format(subj))
    # create the output file
    set_out_file(file_out_path)
    # make a list of all .wav files in that folder
    wav_list = [os.path.join(SUBJ_DIR, file) for file in os.listdir(SUBJ_DIR) if file.endswith(".wav")]
    # recognise each file in that list and output results
    recog_files(wav_list, file_out_path)






