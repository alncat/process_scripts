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
        print self.No, self.name, self.res

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
                switch = 1
                no = int(line[6:11])
                name = line[12:16]
                res = line[17:20]
                chain = line[21:22]
                resNo = int(line[22:26])
                coords = line[30:54].split()
                coord = []
                for i in coords:
                    coord.append(float(i))
            elif line[:6]=="ANISOU" and switch == 1:
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
        r += 0.5*math.log(j[0][p]/i[0][p])
    for k in range(0, 3):
        for l in range(0, 3):
            r += 0.5*(i[0][k]/j[0][l])*pow(np.dot(i[1][k], j[1][l]), 2)
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
    atoms = [ ]
    readpdb(inputpdbfile, atoms)
    contact = []
#    for i in range(0, len(atoms)):
#        atoms[i].display()
    for i in range(0, len(atoms)-1):
        contact.append([])
        for j in range(i+1, len(atoms)):
            di = distance(atoms[i].coord, atoms[j].coord)
            if 0 < di < 5.0:
                contact[i].append([j, di])
    for i in range(0, len(atoms)):
        u = []
        init(u, 3)
        unpack(u, atoms[i].aniso)
        atoms[i].eig = np.linalg.eig(u)
    for i in range(0, len(contact)):
        for j in range(0, len(contact[i])):
            k = contact[i][j][0]
            r = contact[i][j][1]
            kldist = dkl(atoms[i].eig, atoms[k].eig)+dkl(atoms[k].eig, atoms[i].eig) 
            kldist /= 2.0*r
            contact[i][j].append(kldist)
#    for i in range(0, len(contact)):
#        for j in range(0, len(contact[i])):
#            print i, contact[i][j][0], contact[i][j][2], '\n'
    tls = []
    readtls(inputtlsfile, tls)
    if len(tls)==1:
        print "No boundary effect because there is only one tls group"
        exit()

    dklinside = []
    dklboundary = []
    insidestat = []
    initdkl(dklinside, len(tls))
    initdkl(dklboundary, len(tls))
    initdkl(insidestat, len(tls))
    dkli = 0.0
    dklb = []
    switch = 0
    instat = 0
    chain=[]
    for i in range(0, len(tls)):
        if i==0:
            chain.append(tls[i][0])
        else:
            for j in range(0, len(chain)):
                if tls[i][0]!=chain[j] and j==len(chain) - 1:
                    chain.append(tls[i][0]) 
    chainlength = [0]
    lengthtmp=0
    for i in range(0, len(chain)):
        for j in range(0, len(tls)):
            if(tls[j][0]==chain[i]):
                lengthtmp += 1
        chainlength.append(lengthtmp)

#    pdb.set_trace()
    for i in range(0, len(tls)):
        for j in range(0, len(contact)):
            low = atoms[j].resNo - tls[i][1]
            high = atoms[j].resNo - tls[i][2]
            if low < 0:
                pass
            elif low >= 0 and high <= 0 and atoms[j].chain==tls[i][0]:
                for k in range(0, len(contact[j])):
                    pos = contact[j][k][0]
                    for i1 in range(0, len(tls)):
                        if atoms[pos].resNo - tls[i1][1] >= 0 and atoms[pos].resNo - tls[i1][2] <= 0 and atoms[pos].chain==tls[i1][0]:
                            if i1==i:
                                dkli += contact[j][k][2]
                                instat += 1
                                break
                            else:
                                if switch==0:
                                    dklb.append([i1, contact[j][k][2], 1])
                                    switch = 1
                                    break
                                else:
                                    j1 = 0
                                    while j1 < range(0, len(dklb)):
                                        if i1==dklb[j1][0]:
                                            dklb[j1][1] += contact[j][k][2]
                                            dklb[j1][2] += 1
                                            break
                                        elif j1==len(dklb) - 1:
                                            dklb.append([i1, contact[j][k][2], 1])
                                            break
                                        else:
                                            j1 += 1
                                    break
            elif high > 0:
                pass
        dklinside[i] = dkli
        insidestat[i] = instat
        dklboundary[i] = dklb
        dkli = 0.0
        instat = 0
        switch = 0
        dklb = []
        dklinsideav = 0.0
    for i in range(0, len(dklinside)):
        if insidestat[i] != 0:
#            print dklinside[i], insidestat[i], dklinside[i]/insidestat[i]
            dklinside[i] /= insidestat[i]
            dklinsideav += dklinside[i]
    dklinsideav /= len(dklinside)
    print dklinsideav
    dklboundaryav = []
    dklboundaryavtmp = 0.0
    for i in range(0, len(dklboundary)):
        for j in range(0, len(dklboundary[i])):
            dklboundary[i][j].append(dklboundary[i][j][1]/dklboundary[i][j][2])
            dklboundaryavtmp += dklboundary[i][j][3]
            if j == len(dklboundary[i]) - 1:
                dklboundaryavtmp /= len(dklboundary[i])
                dklboundaryav.append(dklboundaryavtmp)
                dklboundaryavtmp = 0.0
    for i in range(0, len(dklboundaryav)):
        dklboundaryavtmp += dklboundaryav[i]
#        print i, dklboundaryav[i]
    dklboundaryavtmp /= len(dklboundary)
    print dklboundaryavtmp
#    for i in range(0, len(dklboundary)):
#        print dklboundary[i]
    bklstat = []
    initdkl(bklstat, len(tls))
    for i in range(0, len(dklinside)-1):
        #        if dklinside[i] >= 1.0E-4:
        for j in range(0, len(dklboundary[i])):
            pos = dklboundary[i][j][0]
            bkltmp = dklboundary[i][j][3]/dklinside[pos]
            bkltmp += dklboundary[i][j][3]/dklinside[i]
            bkltmp /= 2.0
            bklstat[i].append([dklboundary[i][j][0], bkltmp])
    writeout(outputfile, bklstat)
    bklav = 0.0
    bkln = 0
    for i in range(0, len(bklstat)):
        for j in range(0, len(bklstat[i])):
            bklav += bklstat[i][j][1]
            bkln += 1
    if bkln != 0:
        bklav /= bkln

    print bklav
    print bkln
            
    contactsum = 0
    for i in range(0, len(contact)):
        contactsum += len(contact[i])
    print contactsum

    contactsum = 0
    for i in range(0, len(insidestat)):
        contactsum += insidestat[i]
        print dklboundary[i]
        print dklinside[i]
        for j in range(0, len(dklboundary[i])):
            contactsum += dklboundary[i][j][2]
    print contactsum



if __name__ == "__main__":
    main(sys.argv[1:])

