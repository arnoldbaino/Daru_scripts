#!/usr/bin/env python
'''
Created on 6Sep.,2016

@author: 
'''

import argparse, re, os

def main():
    '''Do something'''
    print "running the main routine"
    
    '''Get command line arguments'''
    parser = argparse.ArgumentParser()
    parser.add_argument("bases", help= "Maximum number of bases separating SNPS to be excluded")
    parser.add_argument("vcf", help = "A VCF file")
    parser.add_argument("filtered_vcf", help = "The new filtered VCF file")
    args = parser.parse_args()

    variants = {}
    if os.path.isfile(args.filtered_vcf):
        os.remove(args.filtered_vcf)
    
    with open(args.vcf, "r") as read_file:
        for line in read_file:
            if re.search("^#", line) == None:
                variants[line.split('\t', 2)[1]] = line.rstrip()  
            else:
                with open(args.filtered_vcf, "a") as fpt:
                    fpt.write(line)
    
    first =""
    second = ""
    initial = True
    for variant in sorted(variants, key=int):
        print "working with " + variant
        if first == "":
            first = variant
            print "added " + variant + " to first"
        elif second == "":
            # ASSERT: ONLY GET HERE IF FIRST HAS GOOD OR NO LOWER BOUND
            if (int(variant) - int(first)) > int(args.bases):
                if initial:
                    # INITIAL ENTRY AND GOOD UPPER BOUND
                    print "writing initial " + first
                    with open(args.filtered_vcf, "a") as fpt:
                        fpt.write(variants[first] + "\n")
                    first = variant
                    print "added  " + variant + " to first"
                    initial = False
                else:
                    # FIRST HAS GOOD LOWER AND UPPER BOUND
                    print "writing first " + first
                    with open(args.filtered_vcf, "a") as fpt:
                        fpt.write(variants[first] + "\n")
                    print "added " + variant + " to first"
                    # NEW FIRST HAS GOOD LOWER BOUND
                    first = variant
            else:
                # FIRST HAD GOOD LOWER BOUND BUT BAD UPPER
                print "first and variant too close " + first + " " + variant
                # SECOND HAS BAD LOWER BOUND
                second = variant
                print "added  " + variant + " to second"
                initial = False
        else:
            # ASSERT FIRST HAS BAD UPPER BOUND AND UNKNOWN LOWER
            # ASSERT SECOND HAS HAS BAD LOWER BOUND
            if (int(variant) - int(second)) > int(args.bases):
                print "variant has good lower bound"
                first = variant
                print "added " + variant + " to first"
                second = ""
                print "reset second"
            else:
                # FIRST HAS BAD UPPER AND UNKNOWN LOWER
                # SECOND HAS BAD LOWER AND UPPER BOUND
                # VARIANT HAS BAD LOWER BOUND
                print "have second and variant has bad lower bound"
                first = second
                print "added second " + second + " to first"
                second = variant
                print "added " + variant + " to second"
    
    # ASSERT LAST CASE NEEDS WRITING IF GOOD LOWER BOUND
    # ASSERT SECOND IS SET IF HAS BAD LOWER BOUND 
    if second == "":
        # ASSERT: ONLY GET HERE IF FIRST HAS GOOD OR NO LOWER BOUND
        print "writing first " + first
        with open(args.filtered_vcf, "a") as fpt:
            fpt.write(variants[first] + "\n")      

if __name__ == '__main__':
    main()
