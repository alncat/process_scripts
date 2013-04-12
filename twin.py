#!/usr/local/bin/python

import sys, getopt, math, re, pdb
import numpy as np

class Atom:
    'A class used to describe atom in proteins'
    natom = 0

    def __init__(self, No, name, res, chain, resNo, coord, aniso):
        self.No = No
        self.name = name
        self.res = res
        self.chain = chain
        self.resNo = resNo
        self.coord = coord
        self.aniso = aniso
        self.natom += 1
    def display(self):
        print self.name, self.res, self.resNo

def readpdb(filename, atoms):
    with open(filename, "r") as f:
        line = f.readline()
        while (line != ''):
            if line[:6]=="ATOM  ":
                break
            else:
                line = f.readline()
        while (line != ''):
            if line[:6]=="ATOM  ":
                no = int(line[6:11])
                name = line[12:16]
                res = line[17:20]
                chain = line[21:22]
                resNo = int(line[22:26])
                coords = line[30:54].split()
                coord = []
                switch = 1
                for i in coords:
                    coord.append(float(i))
            elif line[:6]=="ANISOU" and switch==1:
                switch = 0
                anisos = line[28:70].split()
                aniso = []
                for i in anisos:
                    aniso.append(float(i)*0.0001)
                tmp = []
                init(tmp, 3)
                unpack(tmp, aniso)
                if(np.linalg.det(tmp) > 0.0):
                    atom = Atom(no, name, res, chain, resNo, coord, aniso)
                    atoms.append(atom)
                else:
                    pass
#                if name == ' OXT':
#                    break
            line = f.readline()
    f.closed
    return
def readtls(filename, tls):
    with open(filename, "r") as f:
        line = f.readline()
        while (line != ''):
            if re.findall('tls', line):
                seg = re.findall('[A-Z]|[0-9]\d*', line)
                seg[1] = int(seg[1])
                seg[2] = int(seg[2])
                tls.append(seg)
            line = f.readline()
    f.closed
    return

def writeout(filename, out):
    with open(filename, "w") as f:
        for i in range(0, len(out)):
            if out[i]!=[]:
                for j in range(0, len(out[i])):
                    if(j != len(out[i])-1):
                        s1 = str(out[i][j][0])
                        s2 = str(out[i][j][1])
                        s3 = str(i)
                        f.write(s3+'-'+s1+' '+s2+' ')
                    else:
                        s1 = str(out[i][j][0])
                        s2 = str(out[i][j][1])
                        s3 = str(i)
                        f.write(s3+'-'+s1+' '+s2+'\n')
    f.closed
    return

def distance(a, b):
    r = 0.0
    for i in range(0, 3):
        r += pow(a[i]-b[i], 2)
    return math.sqrt(r)

def dkl(i, j):
    r = -1.5
    for p in range(0, 3):
        r += 0.5*np.log(j[0][p]/i[0][p])
    for k in range(0, 3):
        for l in range(0, 3):
            r += 0.5*(i[0][k]/j[0][l])*(np.dot(i[1][k], j[1][l])**2)
    return r
def unpack(u, aniso):
    u[0][0] = aniso[0]
    u[1][1] = aniso[1]
    u[2][2] = aniso[2]
    u[0][1] = aniso[3]
    u[1][0] = aniso[3]
    u[0][2] = aniso[4]
    u[2][0] = aniso[4]
    u[1][2] = aniso[5]
    u[2][1] = aniso[5]
def init(u, i):
    for j in range(0, i):
        u.append([])
        for k in range(0, i):
            u[j].append([])
    return

def initdkl(u, i):
    for j in range(0, i):
        u.append([])
    return

def main(argv):
    inputpdbfile = ''
    inputtlsfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hp:t:o:", ["ipfile=", "itfile=", "ofile="])
    except getopt.GetoptError:
        print 'test.py -p <inputpdbfile> -t <inputtlsfile> -o <outputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'test.py -ip <inputpdbfile> -it <inputtlsfile> -o <outputfile>'
            sys.exit()
        elif opt in ("-p", "--ipfile"):
            inputpdbfile = arg
        elif opt in ("-t", "--itfile"):
            inputtlsfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
#    print 'Input file is "', inputpdbfile
#    print 'Input tls file is "', inputtlsfile 
#    print 'Output file is "', outputfile
    atomsa = []
    readpdb(inputpdbfile, atomsa)
    atomsb = []
    readpdb(inputtlsfile, atomsb)
    contact = []
#    for i in range(0, len(atoms)):
#        atoms[i].display()
    print len(atomsa), len(atomsb)

    for i in range(0, len(atomsa)):
        for j in range(0, len(atomsb)):
            if atomsa[i].name==atomsb[j].name and atomsa[i].resNo==atomsb[j].resNo and atomsa[i].chain==atomsb[j].chain:
                u = []
                init(u, 3)
                unpack(u, atomsa[i].aniso)
                atomsa[i].eig = np.linalg.eig(u)
                unpack(u, atomsb[j].aniso)
                atomsb[j].eig = np.linalg.eig(u)
#                kldist = dkl(atomsa[i].eig, atomsb[j].eig) + dkl(atomsb[j].eig, atomsa[i].eig)
                kldist = dkl(atomsa[i].eig, atomsb[j].eig)
#                kldist /= 2.0
                contact.append(kldist)

    dklav = 0.0
    print len(contact)
    for i in range(0, len(contact)):
        dklav += contact[i]
    dklav /= len(contact)
    print dklav

    for i in range(0, len(atomsb)):
        atomsb[i].display()

if __name__ == "__main__":
    main(sys.argv[1:])

