import sys
import os.path


def writeToFile(f, t, type, s, l, h, v):
    f.write('%s\t%s\t%s\t%s\t%s\t%f\t%f\n'%(t,type,s, l[0],l[1], h, v))


if len(sys.argv) <  2:
    print "Please enter the log file to parse"
    sys.exit()

dir, logFile = os.path.split(sys.argv[1])
if logFile == '':
    logFile = 'log.txt'

inFile = open(os.path.join(dir,logFile), 'r')
outFile = open(os.path.join(dir,"bzflag2.par"), 'w')
objFile = open(os.path.join(dir,"bzflag2.par.str"), 'w')

state = "present"
speed = 0
heading = 0
location = (0,0)
for s in inFile.readlines():
    tokens = s[:-1].split('\t')
    if tokens[2] == 'MOVINGOBJECT_LINEARSPEED' and tokens[3] == 'PandaEPL_avatar':
        speed = float(tokens[4])
    elif tokens[2] == 'SESS_START':
        objFile.write("%s\n"%tokens[4])
    elif tokens[2] == 'OBJ_LOCATION':
        coordinates = tokens[4][tokens[4].find('(')+1:tokens[4].find(')')]
        subTokens = coordinates.split(',')
        objFile.write("%s\t%s\t%s\n"%(tokens[3],subTokens[0],subTokens[1]))
    elif tokens[2] == 'TRAIN_ARRIVE_V':
        state = "train_v"
        writeToFile(outFile,tokens[0],'pickup',state,location,heading,speed)
    elif tokens[2] == 'TRAIN_ARRIVE_I':
        state = "train_i"
        writeToFile(outFile,tokens[0],'pickup',state,location,heading,speed)
    elif tokens[2] == 'TEST_ARRIVE':
        state = "seek"
        writeToFile(outFile,tokens[0],'pickup',state,location,heading,speed)
    elif tokens[2] == 'TRAIN_START':
        state = "train_v"
    elif tokens[2] == 'TEST_START':
        state = "seek"
    elif tokens[2] == 'CAMERA_POS_HEAD':
        heading = float(tokens[4])%360.0

        p3 = tokens[3]
        coordinates = p3[p3.find('(')+1:p3.find(')')]
        subTokens = coordinates.split(',')
        location = (subTokens[0],subTokens[1])

        writeToFile(outFile,tokens[0],'refresh',state,location,heading,speed)
        
inFile.close()
outFile.close()
objFile.close()
