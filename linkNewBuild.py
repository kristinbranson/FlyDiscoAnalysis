#!/usr/local/anaconda/bin/python

import sys
import os
import shutil
import argparse
import subprocess

STAGES = [
    'AutomaticChecks_Incoming',
    'ClassifySex',
    'ComputePerFrameFeatures',
    'DectectIndicatorLedOnOff',
    'JAABADetect',
    'RegisterTrx',
    'TrackWings',
    'MakeCtraxResultsMovie'
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--builddir",default="/groups/branson/home/leea30/git/fba.build/bubble",help="root of build directory")
    parser.add_argument("--newbuild",default="",help="Optional, build dir (rel to builddir) where links will be created")
    #parser.add_argument("dirTgt",help="build dir (rel to dirBase) containing old/existing bins")
    args = parser.parse_args()
        
    if not os.path.exists(args.builddir):
        sys.exit('ERROR: Dir %s not found' % args.builddir)

    # figure out newbuild
    if not args.newbuild:
        builddircontents = os.listdir(args.builddir)
        builds = [el for el in builddircontents if el[0:3]=='201' and len(el)==8]
        builds.sort()
        if not builds:
            sys.exit("ERROR: No builds found in builddir.")
        else:
            print "Found most recent build: " + builds[-1]            
            args.newbuild = builds[-1]

    # figure out current
    currentptr = os.path.join(args.builddir,'current')
    currentptrlink = os.readlink(currentptr)
    print "the current ptr: " + currentptr + " points to: " + currentptrlink

    if not os.path.exists(os.path.join(args.builddir,args.newbuild)):
        sys.exit('ERROR: Dir %s not found' % args.newbuild)
    if not os.path.exists(os.path.join(args.builddir,currentptrlink)):
        sys.exit('ERROR: Dir %s not found' % currentptrlink)

        
    # link binaries that don't exist to last known version (current)
    for stage in STAGES:
        binfile = "FlyBubble" + stage
        runfile = "run_FlyBubble" + stage + ".sh"
        linkIfDNE(args.builddir,args.newbuild,currentptrlink,binfile)
        linkIfDNE(args.builddir,args.newbuild,currentptrlink,runfile)
        
    # update current ptr
    if os.path.exists(currentptr):
        os.remove(currentptr)
    cmd = "ln -s " + args.newbuild + " " + currentptr
    print(cmd)
    subprocess.call(cmd,shell=True)


def linkIfDNE(dirbase,dirnew,dirtgt,fname):
    fnew = os.path.join(dirbase,dirnew,fname)
    if os.path.exists(fnew):
        print "Exists: " + fnew
    else:
        reltgt = os.path.join("..",dirtgt,fname)
        cmd = "ln -s " + reltgt + " " + os.path.join(dirbase,dirnew,fname)
        print(cmd)
        subprocess.call(cmd,shell=True)

if __name__=="__main__":
    main()
