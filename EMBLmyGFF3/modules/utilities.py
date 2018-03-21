#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""

"""


import logging
import sys
import curses.ascii

def print_overwritable(text):
    sys.stderr.write(text)
    sys.stderr.flush()
    sys.stderr.write( ("%c" % curses.ascii.BS) * len(text) ) 

def multiline(prefix, data, sep=";", suffix="", indentprefix = 3, quoted = False, featureType="", wrap=75):
    """
    Creates a multiline entry.

    If data is a list, the entries are listed with "sep" as separator, if data is
    a string, it's quoted over as many lines as needed.
    """
    output=""

    #particular case when RT come empty. We must print ; wihtout quotes
    if(prefix == "RT" and data == ";"):
        output = "%s%s" % (prefix, " "*indentprefix)
        output += str(data)
        return "\n" + output + suffix

    # List Case
    previousChunck=""
    if type(data) == type([]):
        for i, item in enumerate(data):
            if item:
                currentChunck = item + previousChunck
                # If item of the list is too long we have to split it as well
                if len(currentChunck) > wrap:
                    output,lastLine = _splitStringMultiline(output, currentChunck, quoted, wrap)
                    previousChunck="\n"+lastLine

                else:
                    previousChunck=currentChunck

                #Now add separator between chuncks
                if len(previousChunck) >= wrap : # >= Because when previousChunck is last line and is wrapSize char length, adding the \n will give string longer than wrapSize
                    output+=previousChunck+"\n"
                    previousChunck="%s " % sep
                else:
                    previousChunck+="%s " % sep

        output+=previousChunck

    # String case
    else:
        output,lastLine = _splitStringMultiline(output, data, quoted, wrap)
        if len(lastLine) == wrap:
            output+=lastLine+"\n"+sep
        else:
            output+="\n"+lastLine+sep

    #Check if we have output. If not we have to avoid the strip at the end
    doNotStrip=False
    if not output:
        doNotStrip = True

    #Last step: add prefix and middle at each line
    cleanOutput=""
    if output:
        logging.error("output avant print final: %s" % ( output))
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

    if doNotStrip: # Because is only i.e >KW    <
        return "\n" + cleanOutput + suffix
    else:
        return "\n" + cleanOutput.strip().strip(sep) + suffix

# This method allow to wrap a string at a size of wrapSize taking care of quote
# It return back the result in different part: the last line and everything before if exists.
def _splitStringMultiline(output, data, quoted, wrap):
    lastLine=""
    string = " ".join(data.split("\n"))
    output += "\"" if quoted else ""

    roundl=0
    while string:
        roundl+=1
        if roundl == 1: #Within the round 1 the indentation already exists
            if quoted:
                if len(string) + 2 <= wrap: #peculiar case quotes plus string exactly wrapSize
                    lastLine += "\""
                    lastLine = string
                    string = string[len(string):]
                else:# len(string) + 1 > wrap: # + 1 quote
                    splitLoc = _splitWordsMax(string,wrap)
                    line = string[:splitLoc]
                    string = string[len(line):]
                    string=string.strip() # remove white space
                    output += "\""
                    output +=line
            else:
                if len(string) <= wrap:
                    logging.error("round1 if")
                    lastLine = string
                    string = string[len(string):]
                else: # len(string) > wrapSize:
                    logging.error("round1 else")
                    splitLoc = _splitWordsMax(string,wrap)
                    line = string[:splitLoc]
                    string = string[len(line):]
                    string=string.strip() # remove white space
                    output +=line

        else: #Not the first round
            if quoted:
                if len(string)+1 > wrap:
                    splitLoc = _splitWordsMax(string,wrap)
                    line = string[:splitLoc]
                    string = string[len(line):]
                    string=string.strip() # remove white space
                    output +="\n"+line
                else: #it the last round
                    lastLine += string
                    string = string[len(string):] #needed to stop the while
            else:
                if len(string) > wrap:
                    logging.error("roundX if")
                    splitLoc = _splitWordsMax(string,wrap)
                    line = string[:splitLoc]
                    string = string[len(line):]
                    string=string.strip() # remove white space
                    output +="\n"+line
                else: #it the last round
                    logging.error("roundX else - add %s" % string)
                    lastLine +=string
                    string = string[len(string):] #needed to stop the while

    lastLine += "\"" if quoted else ""

    return output,lastLine

def _splitWordsMax(string, valueMax):
    position=0
    positionBefore=0

    words = string.split()
    #If the string was one word longer than the longer value, we have to split the word
    #if len(words) == 1:
    #    words = string.split()
    newString=words.pop(0)
    position = len(newString)
    logging.error("word %s Position %i" % ( newString, position))

    if position >= valueMax:
        logging.error("return value max %i" % ( valueMax))
        return valueMax

    while position <= valueMax :
        positionBefore=position
        newString += " "+words.pop(0)
        position = len(newString)

    return positionBefore