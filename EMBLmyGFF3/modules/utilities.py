#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""

"""

import logging
import sys
import re
import curses.ascii

def print_overwritable(text):
    sys.stderr.write(text)
    sys.stderr.flush()
    sys.stderr.write( ("%c" % curses.ascii.BS) * len(text) ) 

#Wrap is 75 for all line (OC   =5 characters + 75 = 80) except FT
#For FT the wrap is 59 (FT   =5 + 16 characters (feature type + blanc) + 59 characters = 80)
# no_wrap to allow split lines
# splitW to allow split Word
# split_char character to use for split instead of split by word. 
def multiline(prefix, data, sep=";", suffix="", indentprefix = 3, featureType="", wrap=75 ,no_wrap=None, splitW="", split_char=""):
    """
    Creates a multiline entry.

    If data is a list, the entries are listed with "sep" as separator
    """
    #logging.error("prefix: %s featureType: %s data: %s" % (prefix, featureType, data) )
    
    if no_wrap: # equivalent to no wrap when split line deactivated
        wrap = 1000000
    output=""

    # List Case
    previousChunck=""
    if type(data) == type([]):
        #logging.error("list case")
        for i, item in enumerate(data):
            if item:
                currentChunck = item + previousChunck
                # If item of the list is too long we have to split it as well
                if len(currentChunck) > wrap:
                    output,lastLine = _splitStringMultiline(output, currentChunck, wrap, splitW, split_char)
                    previousChunck="\n"+lastLine

                else:
                    previousChunck=currentChunck

                #Now add separator between chuncks
                if (i+1) != len(data): # avoid last iter
                    if len(previousChunck) >= wrap : # >= Because when previousChunck is last line and is wrapSize char length, adding the \n will give string longer than wrapSize
                        output+=previousChunck+"\n"
                        previousChunck=" %s" % sep
                    else:
                        previousChunck+=" %s" % sep
                  
        #In oder to avoid last line longer than expected when adding the suffix
        if len(previousChunck)+len(suffix) > wrap:
            output+=previousChunck+"\n"+suffix
            #logging.error("previousChunckL: %i suffixL %i wrapL %i" % ( len(previousChunck), len(suffix), wrap ) )

        else:
            output+=previousChunck+suffix

    # String case
    else:
        #logging.error("string case")
        output,lastLine = _splitStringMultiline(output, data, wrap, splitW, split_char)

        #add suffix if needed_
        if suffix:
            if len(lastLine)+len(suffix) > wrap:
                # Compile output and lasline as it should
                if not output:
                    output+=lastLine+"\n"+suffix
                else:
                    output+="\n"+lastLine+"\n"+suffix
            else:
                if not output:
                    output+=lastLine+suffix
                else:
                    output+="\n"+lastLine+suffix
        else:
            # Compile output and lasline as it should
            if not output:
                output+=lastLine
            else:
                output+="\n"+lastLine

    #Last step: add prefix and middle at each line
    cleanOutput=""
    if output:

        listLine= output.split("\n")
        for i, line in enumerate(listLine):
            if i == 0:
                if prefix == "FT":
                    cleanOutput += "%s%s" % (prefix, " "*indentprefix) + "%s" % "{:16}".format(featureType) + line
                else:
                    cleanOutput += "%s%s" % (prefix, " "*indentprefix) + line
            else:
                if prefix == "FT":
                    cleanOutput += "\n%s%s" % (prefix, " "*indentprefix) + "%s" % " "*16 + line
                else:
                    cleanOutput += "\n%s%s" % (prefix, " "*indentprefix) + line
    else:
        cleanOutput += "%s%s" % (prefix, " "*indentprefix) #the "+sep" is a trick to keep the final cleaning within the return working properly

    return "\n" + cleanOutput

# This method allow to wrap a string at a size of wrapSize taking care of quote
# It return back the result in different part: the last line and everything before if exists.
def _splitStringMultiline(output, data, wrap, splitW, split_char):
    lastLine=""
    string = " ".join(data.split("\n"))

    roundl=0
    while string:
        roundl+=1
        if roundl == 1: #Within the round 1 the indentation already exists
            if len(string) <= wrap:
                lastLine = string
                string = string[len(string):]
            else: # len(string) > wrapSize:
                splitLoc = _splitWordsMax(string,wrap,splitW,split_char)
                line = string[:splitLoc]
                string = string[len(line):]
                string=string.strip() # remove white space
                output +=line

        else: #Not the first round
            if len(string) > wrap:
                splitLoc = _splitWordsMax(string,wrap,splitW,split_char)
                line = string[:splitLoc]
                string = string[len(line):]
                string=string.strip() # remove white space
                output +="\n"+line
            else: #it the last round
                lastLine +=string
                string = string[len(string):] #needed to stop the while

    return output,lastLine

def _splitWordsMax(string, valueMax, splitW, split_char):
    position=0
    positionBefore=0
    words=[]

    if split_char:
        words = _splitkeepsep(string, split_char)
        #logging.error(words)
    else:
        words = string.split()

    newString=words.pop(0)
    position = len(newString)
    #logging.error("position = %i" % position)

    if position >= valueMax:
        return valueMax

    while position <= valueMax :
        positionBefore=position
        if split_char:
            newString += words.pop(0)
        else:
            newString += " "+words.pop(0)
        position = len(newString)

    return positionBefore

def _splitkeepsep(s, sep):
    return reduce(lambda acc, elem: acc[:-1] + [acc[-1] + elem] if elem == sep else acc + [elem], re.split("(%s)" % re.escape(sep), s), [])