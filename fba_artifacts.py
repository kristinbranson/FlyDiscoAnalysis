#!/usr/local/anaconda/bin/python

import sys
import os
import shutil
import argparse
import re

'''
Use cases for this script:
1. You want to delete all FBA output artifacts starting from stage S.
2. You are planning to rerun stage S, but you want to keep the existing results for comparison. Since rerunning stage S invalidates all downstream stages, you want to archive all FBA output artifacts starting from stage S.

A note regarding registered_trx.mat. This file is created in the reg stage, then overwritten during sex and wtco. In the sex stage, the trx.sex field is populated and sexclassifierinfo variable is appended; if these fields exist already, they are ignored I believe and created fresh. For wtco, the trx.theta field is updated, and here I believe but am less certain that the existing trx.theta data is ignored/irrelevant to the computation. If the existing/incoming trx.theta data is used, then there will be history-dependence wrt the wtco stage; the output artifacts of wtco do not depend only on the most recent wtco run. If this is true, a full/true re-running of wtco really requires a re-running of reg (and consequently sex as well). (An alternative is to archive registered_trx.mat before any run of wtco.)
'''

# dependency graph: DEPGRAPH[stage] gives stages that are immediately dependent on stage
DEPGRAPH = {
    'ctrax': ['reg'],
    'reg': ['led','sex','wtco'],
    'led': [],
    'sex': ['pff'],
    'wtco': ['pff'],
    'pff': [] }

ORIGARTS = ['automatic_checks_incoming_info.mat',
            'automatic_checks_incoming_results.txt',
            'cx__BARCODE.csv',            
            'LOGS',
            'Log.txt',
            'movie.ufmf',
            'protocol.mat',
            'QuickStats.png',
            'QuickStats.txt',
            'StimulusTimingLog.txt',
            'SUCCESS']
ORIGARTSPAT = ['^Metadata\.xml.*$',
              '^FlyBowlDataCaptureParams_.*$']
            
STAGE2ARTS = { 
    'ctrax': [('ctrax_diagnostics.txt','new'),
              ('ctrax_results.mat','new'),
              ('movie.ufmf.ann','new')],
    'led': [('indicatordata.mat','new'),
            ('indicator_log.txt','new')],
    'reg': [('registered_trx.mat','new'),
            ('ledregistrationimage.png','new'),
            ('registrationdata.mat','new'),
            ('registrationdata.txt','new'),
            ('registrationimage.png','new')],
    'sex': [('sexclassifier_diagnostics.txt','new'),
            ('sexclassifier.mat','new'),
            ('registered_trx.mat','overwrite')],
    'wtco': [('wingtracking_results.mat','new'),
             ('registered_trx.mat','overwrite')],
    'pff': [('perframe','new'),
            ('perframefeatures_info.mat','new')]}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("exp",help="full path to FlyBubble experiment")
    parser.add_argument("stage",help="eg ctrax, reg, led, sex, wtco, pff")
    parser.add_argument("arcdir",help="archive subdir")
    parser.add_argument('-d','--dryrun', action='store_true', 
                        help="only print the changes that would be made, leaving the file system untouched")
    parser.add_argument('--LOGS', action='store_true', 
                        help="archive LOGS dir")
    #parser.add_argument("--pre",action="store_true") # goofy/related to 'overwrite' thing, not sure useful yet
    #parser.add_argument("--post",action="store_true") # goofy again, etc
    args = parser.parse_args()
        
    if not os.path.exists(args.exp):
        sys.exit('ERROR: Experiment %s not found' % args.exp)
    if not STAGE2ARTS.has_key(args.stage):
        sys.exit('ERROR: invalid stage %s' % args.stage)

    arcdirfull = os.path.join(args.exp,args.arcdir)
    if not os.path.isdir(arcdirfull):
        if args.dryrun:
            print "os.mkdir(\'{0}\')".format(arcdirfull)
        else:
            os.mkdir(arcdirfull)

#    if args.pre:              
#        preArc(args.exp,args.stage,args.arcdir)
#    if args.post:
#        postArc(args.exp,args.stage,args.arcdir)

    # Compile stages to be archived
    stages = [args.stage]
    iStage = 0
    while iStage<len(stages):
        newstages = DEPGRAPH[stages[iStage]]
        for s in newstages:
            if s not in stages:
                stages.append(s)
        iStage += 1
    print "Stages to archive:"
    for s in stages:
        print s

    # Compile artifacts to archive
    artsmv = []
    for s in stages:
        arts = STAGE2ARTS[s]
        arts = [el[0] for el in arts if el[1]!='overwrite']
        for a in arts:
            if a not in artsmv and os.path.exists(os.path.join(args.exp,a)):
                artsmv.append(a)
    if args.LOGS and os.path.exists(os.path.join(args.exp,'LOGS')):
        artsmv.append('LOGS')

    # move artifacts
    for a in artsmv:
        afull = os.path.join(args.exp,a)
        if args.dryrun:
            print "shutil.move(\'{0}\',\'{1}\')".format(afull,arcdirfull)
        else:
            shutil.move(afull,arcdirfull)

    # check remainder
    unknownFiles = []
    for f in os.listdir(args.exp):
        fMatchesPat = any([re.match(pat,f) for pat in ORIGARTSPAT])
        if (args.dryrun and f not in ORIGARTS and not fMatchesPat and f not in artsmv) or \
           (not args.dryrun and f not in ORIGARTS and not fMatchesPat and f!=args.arcdir):
            unknownFiles.append(f)
    if unknownFiles:
        print args.exp + ': Unknown remainder files: [ ' + ', '.join(unknownFiles) + ' ]'


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
