#!/bin/env python
#
# usage: save.py [options] session
#
# This is a simple script to cycle through a set of simulation files.
# It saves the current set in the sub-directory "Runs" using a sequence
# number stored in the file ".sequence.$rea".  Then, it increments the
# sequence file and copies the field file to the restart file.
#
# $Id$

import os
import stat
import sys
import string

usage = "usage: save.py [options] session"
help  = "options:\n" \
        "-rea     save the .rea file\n" \
        "-chk     use the checkpoint file the field file is empty\n"\
        "-reset   force the sequence number back to 1"

# default suffix lists

clist = [ 'avg' ]
zlist = [ 'his', 'mea', 'sta', 'cp', 'log' ]

# flags

verbose = 0
chk     = 0
reset   = 0
rundir  = "./Runs"
session = None

# Class for storing all of the file names associated with a session

class Session:
    def __init__(self,base):
        self.base = base

        # associated text files
        self.rea  = base + '.rea'
        self.his  = base + '.his'
        self.mea  = base + '.mea'
        self.sta  = base + '.sta'
        self.cp   = base + '.cp'
        self.log  = base + '.log'

        # associated binary files
        self.fld  = base + '.fld'
        self.rst  = base + '.rst'
        self.chk  = base + '.chk'
        self.avg  = base + '.avg'

        # sequencing file
        self.seq  = '.sequence.' + base
        return

    def __repr__(self):
        return self.base
    
# ------------------------------------------------------------------------- 

def exists(path):
    return os.path.exists(path)

def isempty(path):
    st = os.stat(path)
    return st[stat.ST_SIZE] == 0

def init_sequence(session):
    seq = 1
    if exists(session.seq):
        fp  = open(session.seq,"r")
        seq = string.atoi(fp.read())
        fp.close()
    return seq

def roll_sequence(session, seq):
    fp = open(session.seq,"w")
    fp.write("%d\n" % (seq+1))
    fp.close()
    return

def copy(path,seq):
    src = path
    dst = "%s/%s.%03d" % (rundir,path,seq)
    os.system("cp %s %s" % (src,dst))
    if verbose: print "copying:", src, dst
    return

def compress(path,seq):
    src = path
    dst = "%s/%s.%03d" % (rundir,path,seq)
    os.system("cp %s %s ; gzip %s" % (src,dst,dst))
    if verbose: print "zipping:", src, dst+'.gz'
    return
    
def parse_args(argv):
    global verbose, chk, reset
    global session
    global zlist
    
    for c in argv:
        if c[0] == '-':
            if   c == '-v':
                verbose = 1
            elif c == '-rea':
                zlist.append('rea')
            elif c == '-chk':
                chk = 1
            elif c == '-reset':
                reset = 1
            elif c == '-h':
                print usage
                print help
                sys.exit(0)
            else:
                print "save: unknown option --", c[1:]
        else:
            session = Session(c)
    return

# ----------------------------------------------------------------------------

def Save(base):
    global session
    
    if session == None:
        session = Session(base)

    if not exists(session.fld):
        print "save: no field file!"
        sys.exit(1)
        
    # get a sequence number
    if reset == 1:
        seq = 1
    else:
        seq = init_sequence(session)

    # create the run directory?
    if not exists(rundir): os.mkdir(rundir)

    # process field file
    if isempty(session.fld):
        if chk:
            print "Empty field file -- using the checkpoint file"
            copy(session.chk,seq)
        else:
            print "ERROR: empty field file"
            sys.exit(1)
    else:
        copy(session.fld,seq)

    # process remaining files
    for suffix in clist:
        file = session.base+'.'+suffix
        if exists(file):
            copy(file)
            
    for suffix in zlist:
        file = session.base+'.'+suffix
        if exists(file):
            compress(file,seq)

    # update the restart file
    os.system("cp %s %s" % (session.fld, session.rst))

    # update the sequence file
    roll_sequence(session,seq)
    
    return None

# ----------------------------------------------------------------------------

if __name__ == '__main__':
    parse_args(sys.argv[1:])
    
    # prequisites
    if session == None:
        print usage
        sys.exit(1)

    Save(session.base)
