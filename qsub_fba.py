#!/usr/local/anaconda/bin/python

from __future__ import print_function
import sys
import os
import stat
import argparse
import subprocess
import datetime
import re

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--elist",help="file containing list of experiments to process")
    parser.add_argument("--bindir",default="/groups/flyprojects/home/leea30/git/fba.build/bubble/current")
    parser.add_argument("--setdir",default="/groups/flyprojects/home/leea30/git/fba.flybubble/settings")
    parser.add_argument("-ap","--anlsprot",default="current_bubble")
    parser.add_argument("--account",default="bransonk",help="account to charge")
    parser.add_argument("--outsubdir",help="optional output subdir")
    parser.add_argument("--outsubdirauto",action="store_true")
    parser.add_argument("--aci",action="store_true")
    parser.add_argument("--reg",action="store_true")
    parser.add_argument("--sex",action="store_true")
    parser.add_argument("--wgt",action="store_true")
    parser.add_argument("--pff",action="store_true")
    parser.add_argument("--dec",action="store_true")
    parser.add_argument("--jdt",action="store_true")
    parser.add_argument("--mov",action="store_true")
    parser.add_argument("exps",nargs="*",help="full path to experiments to process")

    args = parser.parse_args()
    
    STAGES = ['aci','reg','sex','wgt','pff','dec','jdt','mov']
    stagemap = dict();
    for s in STAGES:
        val = getattr(args,s)
        delattr(args,s)
        stagemap[s] = val
    args.stageon = stagemap        
    
    stage2bin = dict()
    stage2bin['aci'] = "run_FlyBubbleAutomaticChecks_Incoming.sh"
    stage2bin['reg'] = "run_FlyBubbleRegisterTrx.sh"
    stage2bin['sex'] = "run_FlyBubbleClassifySex.sh"
    stage2bin['wgt'] = "run_FlyBubbleTrackWings.sh";
    stage2bin['pff'] = "run_FlyBubbleComputePerFrameFeatures.sh"
    stage2bin['dec'] = "run_FlyBubbleDectectIndicatorLedOnOff.sh"
    stage2bin['jdt'] = "run_FlyBubbleJAABADetect.sh"
    stage2bin['mov'] = "run_FlyBubbleMakeCtraxResultsMovie.sh"
    
    # misc other args maybe settable in future
    args.KEYWORD = "fba"; # used for log/sh filenames, sge job name
    args.MCR = "/groups/branson/bransonlab/projects/olympiad/MCR/v717"
    args.USERNAME = subprocess.check_output("whoami").strip()
    args.TMP_ROOT_DIR = "/scratch/" + args.USERNAME
    args.MCR_CACHE_ROOT = args.TMP_ROOT_DIR + "/mcr_cache_root"
    args.stage2bin = stage2bin
    args.QSUBARGS = "-pe batch 1 -j y -b y -cwd"
    
    # read/check explist
    if args.elist is None:
        exps = args.exps
    else: 
        if args.exps:
            print("elist and exps both supplied; using elist")
        with open(args.elist) as f:
            exps = f.read().splitlines()
    badexps = [e for e in exps if not os.path.exists(e)]
    if badexps:
        print("Bad experiments found:")
        for e in badexps: 
            print(e)
    exps = [e for e in exps if os.path.exists(e)]
    nexps = len(exps)

    # output subdir
    # if outsubdir but no osdauto, then check if exists but just do it
    # if outsubdir and osdautonum, use dt
    # if not outsubdir, then outdir is expdir and use LOGS
    nowstr = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
    if args.outsubdir:
        args.outsubdir = "run." + args.outsubdir
        if args.outsubdirauto:
            args.outsubdir = "{0:s}{1:s}".format(args.outsubdir,nowstr)
        else:
            tf = map(lambda x: os.path.isdir(os.path.join(x,args.outsubdir)),exps)
            nexist = sum(tf)
            if nexist>0:
                print("{0:d} outsubdirs '{1:s}' already exist.".format(nexist,
                                                                       args.outsubdir))

    
    # create template script
    templatescript = "./qsub_template.{0:s}.sh".format(nowstr)
    gencode(templatescript,"<myexp>","<jobid>",args)

    # summarize for user, proceed y/n?
    argsdisp = vars(args).copy()
    del argsdisp['stage2bin']
    del argsdisp['MCR_CACHE_ROOT']
    del argsdisp['TMP_ROOT_DIR']
    del argsdisp['MCR']
    del argsdisp['KEYWORD']
    
    pprintdict(argsdisp)
    print("To process {0:d} exps.".format(nexps))
    resp = raw_input("Template script written to {0:s}. Proceed? y/[n]".format(templatescript))
    if not resp=="y":
        exit()

    DTPAT = '[0-9]{8,8}T[0-9]{6,6}'
    for exp in exps:
        exp.rstrip('/')
        expbase = os.path.basename(exp)

        # jobid
        match = re.search(DTPAT,expbase)
        nowstr = datetime.datetime.now().strftime("%Y%m%dT%H%M%S%f")
        nowstr = nowstr[:-3] # keep only milliseconds
        jobid = match.group() + "-" + nowstr if match else nowstr
        jobid = args.KEYWORD + "-" + jobid
        print(jobid)

        # set outdir, logdir: full paths for output/logs
        # At this point, if outsubdirauto was used then outsubdir has been set
        if args.outsubdir:            
            args.outdirfull = os.path.join(exp,args.outsubdir)
            args.logdirfull = args.outdirfull
            if not os.path.isdir(args.outdirfull):
                os.mkdir(args.outdirfull)
        else:
            args.outdirfull = None
            args.logdirfull = os.path.join(exp,'LOGS')

        # generate code
        shfile = os.path.join(args.logdirfull,"{0:s}.sh".format(jobid))
        logfile = os.path.join(args.logdirfull,"{0:s}.log".format(jobid))
        gencode(shfile,exp,jobid,args)

        # submit 
        qargs = "-A {0:s} -o {1:s} -N {2:s} {3:s} {4:s}".format(args.account,
                     logfile,jobid,args.QSUBARGS,shfile)
        qsubcmd = "qsub " + qargs
        print(qsubcmd)
        subprocess.call(qsubcmd,shell=True)

    exit()

def gencodehelp(f,expdir,args,stage):
    if args.stageon[stage]:
        if args.outsubdir:
            print(args.bindir + "/" + args.stage2bin[stage],args.MCR,
                  expdir,"settingsdir",args.setdir,"analysis_protocol",
                  args.anlsprot,"outdir",args.outsubdir,file=f)
        else:
            print(args.bindir + "/" + args.stage2bin[stage],args.MCR,
                  expdir,"settingsdir",args.setdir,"analysis_protocol",
                  args.anlsprot,file=f)
    

def gencode(fname,expdir,jobid,args):
    f = open(fname,'w')
    print("#!/bin/bash",file=f)
    print("",file=f)
    print("source ~/.bashrc",file=f)
    print("umask 002",file=f)
    print("export MCR_CACHE_ROOT="+args.MCR_CACHE_ROOT + "." + jobid,file=f)
    print("echo $MCR_CACHE_ROOT",file=f)

    print("",file=f)
    gencodehelp(f,expdir,args,'aci')
    gencodehelp(f,expdir,args,'reg')
    gencodehelp(f,expdir,args,'sex')
    gencodehelp(f,expdir,args,'wgt')
    gencodehelp(f,expdir,args,'pff')
    gencodehelp(f,expdir,args,'dec')
    gencodehelp(f,expdir,args,'jdt')
    gencodehelp(f,expdir,args,'mov')
    print("",file=f)

    print("rm -rf",args.MCR_CACHE_ROOT+"."+jobid,file=f)
    f.close()
    os.chmod(fname,stat.S_IRUSR|stat.S_IXUSR|stat.S_IRGRP|stat.S_IXGRP|stat.S_IROTH);

def pprintdict(d, indent=0):
   for key, value in sorted(d.items()):
      print('\t' * indent + str(key))
      if isinstance(value, dict):
         pprintdict(value, indent+1)
      else:
         print('\t' * (indent+1) + str(value))

if __name__=="__main__":
    main()

