#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
import sys
from IPython.core.debugger import Tracer

"""
Location examples (http://www.insdc.org/files/feature_table.html)

The following is a list of common location descriptors with their meanings: 

Location                  Description   

467                       Points to a single base in the presented sequence 

340..565                  Points to a continuous range of bases bounded by and
                          including the starting and ending bases

<345..500                 Indicates that the exact lower boundary point of a feature
                          is unknown.  The location begins at some  base previous to
                          the first base specified (which need not be contained in 
                          the presented sequence) and continues to and includes the 
                          ending base 

<1..888                   The feature starts before the first sequenced base and 
                          continues to and includes base 888

1..>888                   The feature starts at the first sequenced base and 
                          continues beyond base 888

102.110                   Indicates that the exact location is unknown but that it is 
                          one of the bases between bases 102 and 110, inclusive

123^124                   Points to a site between bases 123 and 124

join(12..78,134..202)     Regions 12 to 78 and 134 to 202 should be joined to form 
                          one contiguous sequence


complement(34..126)       Start at the base complementary to 126 and finish at the 
                          base complementary to base 34 (the feature is on the strand 
                          complementary to the presented strand)


complement(join(2691..4571,4918..5163))
                          Joins regions 2691 to 4571 and 4918 to 5163, then 
                          complements the joined segments (the feature is on the 
                          strand complementary to the presented strand) 

join(complement(4918..5163),complement(2691..4571))
                          Complements regions 4918 to 5163 and 2691 to 4571, then 
                          joins the complemented segments (the feature is on the 
                          strand complementary to the presented strand)
  
J00194.1:100..202         Points to bases 100 to 202, inclusive, in the entry (in 
                          this database) with primary accession number 'J00194'
 
join(1..100,J00194.1:100..202)
                          Joins region 1..100 of the existing entry with the region
                          100..202 of remote entry J00194

   /\
  /| \  We do not deal with features that behave like join(complement(4918..5163),complement(2691..4571)). 
 / .  \ If a feature is on the negative strand we will always deal like complement(join(2691..4571,4918..5163))
/______\
"""

class Span(object):
    
    def __init__(self, start , end = None):
        self.complement = False
        self.complement = False
        self.start = None
        self.end = None
        self.join = ".."  # TODO: support other join operators!!
        self.accession = ""
        self.start = start
        self.end = end
    
    def __repr__(self):
        output = ""
#        if self.start > self.end and not self.complement:
#            self.complement = True
        
#        start = self.end if self.complement else self.start
        start = self.start    
#        end = self.start if self.complement else self.end
        end = self.end
#        start_mod = self.end_modifier if self.complement else self.start_modifier
#        end_mod = self.start_modifier if self.complement else self.end_modifier
        if start != None:
#           output += start_mod if start_mod else ""
            output += "%s" % start
        if end:
            output += self.join
#            output += end_mod if end_mod else ""
            output += "%s" % end
        if self.accession:
            output = "%s:%s" % (self.accession, output)
 #       if self.complement:
 #           output = "complement(%s)" % output
        return output
      
    def set_accession(self, value):
        self.accession = value
    
    def set_any(self):
        self.join = "."
    
    def set_between(self):
        self.join = "^"
        if self.end - self.start != 1:
            print "WARNING: Setting operator 'between' for non-adjoining start-end values"
    
    def set_complement(self):
         self.complement = True

class Location(object):
    """
    The location descriptor can be one of the following: 
    (a) a single base number, e.g. "467"
    (b) a site between two indicated adjoining bases, e.g. "123^124"
    (c) a single base chosen from within a specified range of bases (not allowed for new
        entries), e.g. "102.110"
    (d) the base numbers delimiting a sequence span, e.g. "240..565"
    (e) a remote entry identifier followed by a local location descriptor
        (i.e., a-d), e.g. "J00194.1:100..202"
    
    
    multiple locations can be joined as "join(12..78,134..202)".
    locations on the complementary strand can be presented as "complement(34..126)"
    'fuzzy' ends can be represented as "<1..>888"
    """
    
    def __init__(self, feature, *args, **kwargs):
        
        self.spans = []
        #sys.stderr.write("deux et demi part: '%s' \n" % feature)
        for arg in args:
            nb_part=len(arg.parts)
            if hasattr(arg, "start"): # Deal with BioPython FeatureLocation basically
                count_part=0
                for part in arg.parts:

                    count_part+=1
                    # We assume the position arrive sorted !!!!!! croissant
                    start=part.start+1    
                    # there is no start, it is positive strand and the first exon
                    if('has_start' in feature.qualifiers):
                        if feature.qualifiers['has_start'] == "no" and  count_part == 1 and hasattr(part, "strand") and part.strand > 0:
                            start="<"+str(start)
                    if('has_stop' in feature.qualifiers):
                        if feature.qualifiers['has_stop'] == "no" and  count_part == 1 and hasattr(part, "strand") and part.strand < 0:
                            start="<"+str(start)
                    span = Span(start)

                    #Add the complement information
                    if hasattr(part, "strand") and part.strand < 0:
                        #Tracer()()
                        span.set_complement()

                    if hasattr(part, "end"):
                        end=part.end
                        span.end = end
                        if('has_stop' in feature.qualifiers):
                            if feature.qualifiers['has_stop'] == "no" and  count_part == nb_part and hasattr(part, "strand") and part.strand > 0:
                                end=">"+str(end)
#                            sys.stderr.write("There is end: '%s' \n" % part)
                        if('has_start' in feature.qualifiers):
                            if feature.qualifiers['has_start'] == "no" and  count_part == nb_part and hasattr(part, "strand") and part.strand < 0:
                                end=">"+str(end)
                        span.end = end   

                    if hasattr(part, "ref") and hasattr(part, "ref_db") and part.ref and part.ref_db:
                        span.accession = "%s:%s" % (part.ref_db, part.ref)
                    self.spans += [span]




            elif not self.spans:
                self.spans += [Span(arg)]
            elif arg == "any":
                self.spans[-1].set_any()
            elif arg == "between":
                self.spans[-1].set_between()
            elif self.spans[-1].end != None:
                self.spans += [Span(arg)]
            else:
                self.spans[-1].set_end(arg)
    
    def __repr__(self):
        output = ""
        if not self.spans:
            return "."
        if len(self.spans) > 1:
            output = "join(%s)" % (",".join(map(str, self.spans)))
        else:
            output = "".join(map(str, self.spans))
        if self.spans[0].complement:       #added JD - We just check the firt location if complement is true we add the information
            output = "complement(%s)" % output
        return output




if __name__ == '__main__':
    
    print " ----- Feature Location implementation for EMBL format ----- "
    print
    print "Location(\"<23\", \"700\")               =>", Location("<23", "700")
    print "Location(\"<23\", 700, \">332\", \"<931\") =>", Location("<23", 700, ">332", "<931")
    print "Location(\">2321\", 700, \"332\")        =>", Location(">2321", 700, "332")
    print "Location(\">2321\", \"single\", \"332\")   =>", Location(">2321", "single", "332")
    print "Location(\">23\", 166, \"any\")          =>", Location(">23", 166, "any")
    print "Location(\"165\", 166, \"between\")      =>", Location("165", 166, "between")
    
    try:
        from Bio.SeqFeature import FeatureLocation, BeforePosition, AfterPosition
        print "location = FeatureLocation(BeforePosition(5), AfterPosition(33))"
        location = FeatureLocation(BeforePosition(5), AfterPosition(33))
        print "Location( location, 100, \"<50\")      =>", Location( location, 100, "<50")
    except:
        pass
    