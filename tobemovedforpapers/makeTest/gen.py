import argparse
import os
os.chdir("/home/teocollin/gitcode/diderot/tests/makeTest")
parser = argparse.ArgumentParser(description='Make a Diderot program wrap setup')
parser.add_argument("--shell", type=str, default="/bin/sh")
parser.add_argument("--python", type=str, default="python3")
parser.add_argument("--pythonIncludes", type=str,
                    default="-I/usr/include/python3.6m -I/usr/include/python3.6m")
parser.add_argument("programName", type=str, default="prog")
parser.add_argument("diderotPath", type=str)  # /home/teocollin/gitcode/diderot
parser.add_argument("teemPath", type=str)  # ugg
parser.add_argument("dir", type=str)  # actual argument

args = parser.parse_args()
targetDir = args.dir
programName  = args.programName
argsp = {'shell': args.shell,
         "python": args.python,
         "pythonIncludes": args.pythonIncludes,
         "programName": args.programName,
         "diderotPath": args.diderotPath,
         "teemPath": args.teemPath,}
os.mkdir(targetDir)
with open("Makefile.in") as f:
    makeFileInsert = f.read()
    newMakeFile = makeFileInsert.format(**argsp)
    with open(targetDir + "/Makefile", "w+") as g:
        g.write(newMakeFile)
with open("runData.py.in") as f:
    runData = f.read()
    newRunFile = runData.format(**argsp)
    with open(targetDir + "/runData.py", "w+") as g:
        g.write(newRunFile)
with open("runStatic.py.in") as f:
    runStat = f.read()
    newRunStat = runStat.format(**argsp)
    with open(targetDir + "/runStatic.py", "w+") as g:
        g.write(newRunStat)

with open("prog.diderot.in") as f:
    did = f.read()
    newDid = did.format(**argsp)
    with open(targetDir + "/{0}.diderot".format(programName), "w+") as g:
        g.write(newDid)
