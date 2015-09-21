#!/usr/local/anaconda/bin/python

import sys
import os
import shutil
import argparse

STAGE2ARTS = {'reg': [('registered_trx.mat','overwrite'),('ledregistrationimage.png','new'),('registrationdata.mat','new'),('registrationdata.txt','new'),('registrationimage.png','new')] }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("exp",help="full path to FlyBubble experiment")
    parser.add_argument("stage",help="eg reg")
    parser.add_argument("arcdir",help="archive subdir")
    parser.add_argument("--pre",action="store_true")
    parser.add_argument("--post",action="store_true")
    args = parser.parse_args()
        
    if not os.path.exists(args.exp):
        sys.exit('ERROR: Experiment %s not found' % args.exp)
    arcdirfull = os.path.join(args.exp,args.arcdir)
    if not os.path.isdir(arcdirfull):
        os.mkdir(arcdirfull)

    if args.pre:
        preArc(args.exp,args.stage,args.arcdir)
    if args.post:
        postArc(args.exp,args.stage,args.arcdir)


def preArc(exp,stage,arcdir):
    arts = STAGE2ARTS[stage]
    for (afile,amode) in arts:
        if amode=='overwrite':
            afullsrc = os.path.join(exp,afile)
            if os.path.isfile(afullsrc):
                afulldst = os.path.join(exp,arcdir,afile+'.pre')
                print "copying " + afullsrc + " to " + afulldst
                shutil.copy2(afullsrc,afulldst)
            else:
                print "warning: file " + afullsrc + " doesn't exist"

def postArc(exp,stage,arcdir):
    arts = STAGE2ARTS[stage]
    for (afile,amode) in arts:
        afullsrc = os.path.join(exp,afile)
        if os.path.isfile(afullsrc):
            afulldst = os.path.join(exp,arcdir,afile)
            if amode=='overwrite':
                print "copying " + afullsrc + " to " + afulldst
                shutil.copy2(afullsrc,afulldst)
            else:
                print "moving " + afullsrc + " to " + afulldst
                shutil.move(afullsrc,afulldst)
        else:
            print "warning: file " + afullsrc + " doesn't exist"            
            
if __name__=="__main__":
    main()
