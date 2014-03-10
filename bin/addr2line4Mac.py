#! /usr/bin/python

import sys
filename=sys.argv[1]
address=sys.argv[2]
import re
from os import environ,path

fullFile=None
if path.exists(filename):
    fullFile=filename

for v in ["PATH","LD_LIBRARY_PATH"]:
    if not fullFile:
        for d in environ[v].split(':'):
            if path.exists(path.join(d,filename)):
                fullFile=path.join(d,filename)
                break

if not fullFile:
    fullFile=filename

answer="??:0"

if path.exists(fullFile):
    import subprocess

    result=subprocess.Popen(["xcrun", "atos",
                             "-o",fullFile,
                             address],
                            stdout=subprocess.PIPE
                        ).communicate()[0]
    match=re.compile('.+ \((.+)\) \((.+)\)').match(result)
    if match:
        answer=match.group(2)+" "+match.group(1)
    else:
        import os
        result=subprocess.Popen(["xcrun", "atos",
                                 "-p",str(os.getppid()),
                                 address],
                                stdout=subprocess.PIPE
                            ).communicate()[0]
        match=re.compile('.+ \((.+)\) \((.+)\)').match(result)
        if match:
            answer=match.group(2)+" "+match.group(1)

print answer,

sys.exit(255)
