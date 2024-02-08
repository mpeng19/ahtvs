import math
import numpy


def overlaps(MGF, minusHOMO=0, plusLUMO=0):
    with open(MGF,'r') as f:
        NAtoms=int(f.readline().split()[0])
        Atoms = []
        Coords = []

        NOrbs = 0
        for i in range(NAtoms):
            readin=f.readline().split()[0:4]
            Atoms.append( int(readin[0]) )
            Coords.append([float(readin[1]),float(readin[2]),float(readin[3])])

        NElect= sum(Atoms)-2*len([i for i in Atoms if i>2])-8*len([i for i in Atoms if i>10])
        NOrbs=NAtoms+3*len([i for i in Atoms if i>1])
        # Zetas = []
        # for i in range(NAtoms):
        #     readin=map(float,f.readline().split())
        #     Zetas.append(readin)
        # ZetasList = [Zetas[i][j] for i in range(len(Zetas)) for j in [0,1,1,1,2,2,2,2,2] if Zetas[i][j]>0.001  ]


        line = 'aaa'
        stopper=str(NElect/2-minusHOMO-1)+'A'
        #print stopper
        while stopper not in line :
            line=f.readline()
        #     print line
        # print 'HOMO-1 En', line[len(stopper)+3:-1]

        CoefsHOMOm=[]
        for i in range((NOrbs+4)/5):
            readin =  f.readline()
            CoefsHOMOm.append([float(readin[15*n:15*n+15].replace("D","E"))  for n in range(len(readin)/15)])
        CoefsHOMOm=numpy.array([val for subl in CoefsHOMOm for val in subl])


        stopper=str(NElect/2-minusHOMO)+'A'
        #print stopper
        while stopper not in line :
            line=f.readline()
        #     print line
        # print 'HOMO En', line[len(stopper)+3:-1]
        HOMOEn = float(line.split()[-1])

        CoefsHOMO=[]
        for i in range((NOrbs+4)/5):
            readin =  f.readline()
            CoefsHOMO.append([float(readin[15*n:15*n+15].replace("D","E"))  for n in range(len(readin)/15)])
        CoefsHOMO=numpy.array([val for subl in CoefsHOMO for val in subl])

        stopper=str(NElect/2+1+plusLUMO)+'A'
        # print stopper
        while stopper not in line :
            line=f.readline()
        #print 'LUMO En', line[len(stopper)+3:-1]
        LUMOEn = float(line.split()[-1])

        CoefsLUMO = []
        for i in range((NOrbs+4)/5):
            readin =  f.readline()
            CoefsLUMO.append([float(readin[15*n:15*n+15].replace("D","E"))  for n in range(len(readin)/15)])
        CoefsLUMO=numpy.array([val for subl in CoefsLUMO for val in subl])

        stopper=str(NElect/2+2+plusLUMO)+'A'
        # print stopper
        while stopper not in line :
            line=f.readline()
        # print 'LUMO+1 En', line[len(stopper)+3:-1]

        CoefsLUMOp = []
        for i in range((NOrbs+4)/5):
            readin =  f.readline()
            CoefsLUMOp.append([float(readin[15*n:15*n+15].replace("D","E"))  for n in range(len(readin)/15)])
        CoefsLUMOp = numpy.array([val for subl in CoefsLUMOp for val in subl])

    HLoverlap = sum([math.sqrt(CoefsHOMO[i]*CoefsHOMO[i]*CoefsLUMO[i]*CoefsLUMO[i]) for i in range(len(CoefsHOMO)) ])
    HLoverlapSQ = sum([CoefsHOMO[i]*CoefsHOMO[i]*CoefsLUMO[i]*CoefsLUMO[i] for i in range(len(CoefsHOMO)) ])
    return(HOMOEn, LUMOEn , HLoverlap, HLoverlapSQ)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='HOMO-LUMO overlap from MOPAC mgf file')
    parser.add_argument('MGF', type=str, help='mgf file to analyse')
    parser.add_argument('-H', dest='minusHOMO', type=int, help='Occupied Orbital below HOMO, 0 for HOMO', default=0)
    parser.add_argument('-L', dest='plusLUMO', type=int, help='Occupied Orbital above LUMO, 0 for LUMO', default=0)
    args = parser.parse_args()
    minusHOMO=args.minusHOMO
    plusLUMO=args.plusLUMO

    print Overlaps(args.MGF, 1, 1)




        line = 'aaa'
        while line[:8] != ' INVERSE' :
            line=f.readline()
        SMat=[['AA' for i in range(len(ZetasList)) ] for i in range(len(ZetasList)) ]
        awful = [i.strip('\n').replace("D","E") for i in f.readlines()[:-1]]
        better = [float(i[15*n:15*n+15]) for i in awful  for n in range(len(i)/15)  ]



    k=0
    for i  in range(len(CoefsHOMO)):
        for j  in range(i+1):
            SMat[i][j]=SMat[j][i]=better[k]
            k+=1

    SMat=numpy.array(SMat)
    print SMat



    SMat=numpy.dot(SMat,SMat)
    SMat=numpy.linalg.inv(SMat)





    print sum([CoefsHOMO[i]*CoefsHOMO[i] for i in range(len(CoefsHOMO)) ])
    print sum([CoefsLUMO[i]*CoefsLUMO[i] for i in range(len(CoefsLUMO)) ])
    print sum([CoefsHOMO[i]*CoefsLUMO[i] for i in range(len(CoefsHOMO)) ])
    print 'H   L   Overlap', sum([math.sqrt(CoefsHOMO[i]*CoefsHOMO[i]*CoefsLUMO[i]*CoefsLUMO[i]) for i in range(len(CoefsHOMO)) ])
    print 'H-1 L   Overlap', sum([math.sqrt(CoefsHOMOm[i]*CoefsHOMOm[i]*CoefsLUMO[i]*CoefsLUMO[i]) for i in range(len(CoefsHOMO)) ])
    print 'H   L+1 Overlap', sum([math.sqrt(CoefsHOMO[i]*CoefsHOMO[i]*CoefsLUMOp[i]*CoefsLUMOp[i]) for i in range(len(CoefsHOMO)) ])
    print 'H-1 L+1 Overlap', sum([math.sqrt(CoefsHOMOm[i]*CoefsHOMOm[i]*CoefsLUMOp[i]*CoefsLUMOp[i]) for i in range(len(CoefsHOMO)) ])

    print 'H   L   OverlapB', sum([CoefsHOMO[i]*CoefsHOMO[i]*CoefsLUMO[i]*CoefsLUMO[i] for i in range(len(CoefsHOMO)) ])
