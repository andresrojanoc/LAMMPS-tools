#! python3

"""
This program updates LAMMPS' configuration file according to a dump file
"""
from datetime import datetime
import sys, logging, time
start = time.time()

FORMAT = '%(asctime)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.DEBUG,format = FORMAT)

"""
Leave the following line uncommented to avoid Info. messages during the run
"""
#logging.disable(logging.CRITICAL)

class DumpFrame:
    def __init__(self):
        self.timestep = None
        self.numberAtoms = None
        self.boxBounds = None
        self.xlo = None
        self.xhi = None
        self.ylo = None
        self.yhi = None
        self.zlo = None
        self.zhi = None
        self.atomInfo = None
        self.atomData = []

class Atom:
    def __init__(self,id,mol,type,q,x,y,z,nx,ny,nz):
        self.id = id
        self.mol = mol
        self.type = type
        self.q = q
        self.x = x
        self.y = y
        self.z = z
        self.nx = nx
        self.ny = ny
        self.nz = nz

class Velocity:
    def __init__(self, id,vx, vy,vz):
        self.id = id
        self.vx =vx
        self.vy =vy
        self.vz =vz

class Bond:
    def __init__(self,type,bond1,bond2):
        self.type = type
        self.bond1 = bond1
        self.bond2 = bond2

class Angle:
    def __init__(self,type,angle1,angle2,angle3):
        self.type = type
        self.angle1 = angle1
        self.angle2 = angle2
        self.angle3 = angle3

class Dihedral:
    def __init__(self,type,dihedral1,dihedral2,dihedral3,dihedral4):
        self.type = type
        self.dihedral1 = dihedral1
        self.dihedral2 = dihedral2
        self.dihedral3 = dihedral3
        self.dihedral4 = dihedral4

class LmpConf:
    def __init__(self):
        self.xlo = None
        self.xhi = None
        self.ylo = None
        self.yhi = None
        self.zlo = None
        self.yhi = None
        self.nAtoms = None
        self.nBonds = None
        self.nAngles = None
        self.nDihedrals = None
        self.nAtomTypes = None
        self.nBondTypes = None
        self.nAngleTypes = None
        self.nDihedralTypes = None
        self.mass = []
        self.atoms = []
        self.velocities = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.blen = []
        self.bhlf = []
        self.blni = []
        self.molecules = []

    def initialize_box(self):
        self.blen.append(abs(self.xhi-self.xlo))
        self.blen.append(abs(self.yhi-self.ylo))
        self.blen.append(abs(self.zhi-self.zlo))
        for i in range(3):
            self.bhlf.append(self.blen[i]/2.0)
            self.blni.append(1.0/self.blen[i])

    def apply_pbcs(self,dr,dim):
        if dr > self.bhlf[dim]:
            dr = dr - self.blen[dim] * float(int(dr * self.blni[dim] + 0.5))
        elif dr < -self.bhlf[dim]:
            dr = dr - self.blen[dim] * float(int(dr * self.blni[dim] - 0.5))
        return dr

    def read_configuration_file(self,confFile):
        try:
            file = open(confFile,"r")
            line = file.readline()
            while line:
                if "atoms" in line:
                    self.nAtoms = int(line.strip().split()[0])
                if "bonds" in line:
                    self.nBonds = int(line.strip().split()[0])
                if "angles" in line:
                    self.nAngles = int(line.strip().split()[0])
                if "dihedrals" in line:
                    self.nDihedrals = int(line.strip().split()[0])
                if "atom types" in line:
                    self.nAtomTypes = int(line.strip().split()[0])
                if "bond types" in line:
                    self.nBondTypes = int(line.strip().split()[0])
                if "angle types" in line:
                    self.nAngleTypes = int(line.strip().split()[0])
                if "dihedral types" in line:
                    self.nDihedralTypes = int(line.strip().split()[0])
                if "xlo xhi" in line:
                    self.xlo = float(line.strip().split()[0])
                    self.xhi = float(line.strip().split()[1])
                if "ylo yhi" in line:
                    self.ylo = float(line.strip().split()[0])
                    self.yhi = float(line.strip().split()[1])
                if "zlo zhi" in line:
                    self.zlo = float(line.strip().split()[0])
                    self.zhi = float(line.strip().split()[1])
                if "Masses" in line :
                    file.readline()
                    for i in range(self.nAtomTypes):
                        l = file.readline().strip().split()
                        self.mass.append(float(l[1]))
                if "Pair Coeffs" in line:
                    logging.info('This script is not designed to read Pair Coeffs data')
                    sys.exit()
                if "Bond Coeffs" in line:
                    file.readline()
                    for i in range(self.nBondTypes):
                        logging.info('This script is not designed to read Bond Coeffs data')
                        sys.exit()
                if "Angle Coeffs" in line:
                    file.readline()
                    for i in range(self.nAngleTypes):
                        logging.info('This script is not designed to read Angle Coeffs data')
                        sys.exit()
                if "Dihedral Coeffs" in line:
                    file.readline()
                    for i in range(self.nDihedralTypes):
                        logging.info('This script is not designed to read Dihedral Coeffs data')
                        sys.exit()

                if "Atoms" in line:
                    file.readline()
                    for i in range(self.nAtoms):
                        l = file.readline().strip().split()
                        if len(l) == 7:
                            nx,ny,nz = 0,0,0
                            self.atoms.append(Atom(int(l[0]),int(l[1]),int(l[2]),float(l[3]),float(l[4]), float(l[5]),float(l[6]),nx,ny,nz))
                        elif len(l) == 10:
                            self.atoms.append(Atom(int(l[0]),int(l[1]),int(l[2]),float(l[3]),float(l[4]), float(l[5]),float(l[6]),int(l[7]),int(l[8]),int(l[9])))
                        else:
                            logging.info('Inconsistent number of entries in Atoms line %d for molecular atom_style' % i+1)
                            sys.exit()

                if "Velocities" in line:
                    file.readline()
                    for i in range(self.nAtoms):
                        l = file.readline().strip().split()
                        self.velocities.append(Velocity(int(l[0]),float(l[1]),float(l[2]),float(l[3])))


                if "Bonds" in line:
                    file.readline()
                    for i in range(self.nBonds):
                        l = file.readline().strip().split()
                        self.bonds.append(Bond(int(l[1]),int(l[2]),int(l[3])))

                if "Angles" in line:
                    file.readline()
                    for i in range(self.nAngles):
                        l = file.readline().strip().split()
                        self.angles.append(Angle(int(l[1]),int(l[2]),int(l[3]),int(l[4])))

                if "Dihedrals" in line:
                    file.readline()
                    for i in range(self.nDihedrals):
                        l = file.readline().strip().split()
                        self.dihedrals.append(Dihedral(int(l[1]),int(l[2]),int(l[3]),int(l[4]),int(l[5])))
                line = file.readline()

        except IOError:
            logging.info('Could not open %s' % confFile)
            sys.exit()
        self.initialize_box()
        file.close()

    def read_dump_file(self,dumpFile):
        try:
            file = open(dumpFile,"r")
        except IOError:
            logging.info('Could not open %s' % dumpFile)
            sys.exit()
        else:
            localDumps = []
            line = file.readline()
            while line:
                if "ITEM: TIMESTEP" in line:
                    frame = DumpFrame()
                    frame.timestep = int(file.readline().strip().split()[0])
                    self.timestep =  frame.timestep
                if "ITEM: NUMBER OF ATOMS" in line:
                    frame.numberAtoms = int(file.readline().strip().split()[0])
                    number = frame.numberAtoms
                if "ITEM: BOX BOUNDS" in line:
                    line = line.strip().split()
                    b1 = line[-3]
                    b2 = line[-2]
                    b3 = line[-1]
                    frame.boxBounds = "%s %s %s" % (b1,b2,b3)
                    line = file.readline().strip().split()
                    frame.xlo = float(line[0])
                    frame.xhi = float(line[1])
                    line = file.readline().strip().split()
                    frame.ylo = float(line[0])
                    frame.yhi = float(line[1])
                    line = file.readline().strip().split()
                    frame.zlo = float(line[0])
                    frame.zhi = float(line[1])
                if "ITEM: ATOMS" in line:
                    frame.atomInfo = line.strip().replace("ITEM: ATOMS ","").split()
                    for i in range(frame.numberAtoms):
                        frame.atomData.append(file.readline().strip().split())
                    localDumps.append(frame)
                line = file.readline()
            file.close()

            """[-1] is for using the last dump change the value to the desired dump"""
            return localDumps[-1]


    def update_configuration_file(self,dumpData,confFile):
        # First, build a map from dump to conf
        siteMap = {}
        # Sense strand
        for i in range(dumpData.numberAtoms+1): # plus 1 for including the last value
            siteMap[i] = i
        # Reorder the data in dumpData
        info = dumpData.atomInfo
        infoDict = {}
        for i in range(len(info)):
            infoDict[info[i]]=  int(i)
        final = dumpData.numberAtoms
        for i in range (final):
            data = dumpData.atomData[i]
            index = int(data[infoDict["id"]])
            if index in siteMap:
                idx = siteMap[index]  # Get the index in the configuration file (assuming everything is in order)

                """ For further customization activate the if statement and indent the following lines """
               # if float(data[infoDict["type"]]) == 1 :
                self.atoms[idx-1].x = float(data[infoDict["x"]])
                self.atoms[idx-1].y = float(data[infoDict["y"]])
                self.atoms[idx-1].z = float(data[infoDict["z"]])

                """ Note: It is assumed that the velocities and atoms are in the same number """
                self.velocities[idx-1].vx = float(data[infoDict["vx"]])
                self.velocities[idx-1].vy = float(data[infoDict["vy"]])
                self.velocities[idx-1].vz = float(data[infoDict["vz"]])


    def write_configuration_file(self,outFile):
        file = open(outFile,"w")
        file.write("""LAMMPS configuration file created using %s on %s \n\n""" % (sys.argv[0],datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
        file.write("\t%ld atoms\n" % len(self.atoms))
        file.write("\t%ld bonds\n" % self.nBonds)
        file.write("\t%ld angles\n" % self.nAngles)
        file.write("\t%ld dihedrals\n\n" % self.nDihedrals)

        file.write("\t%ld atom types\n" % self.nAtomTypes)
        file.write("\t%ld bond types\n" % self.nBondTypes)
        file.write("\t%ld angle types\n" % self.nAngleTypes)
        file.write("\t%ld dihedral types\n\n" % self.nDihedralTypes)
        file.write("\t%lf\t%lf xlo xhi\n" % (self.xlo,self.xhi))
        file.write("\t%lf\t%lf ylo yhi\n" % (self.ylo,self.yhi))
        file.write("\t%lf\t%lf zlo zhi\n\n" % (self.zlo,self.zhi))

        logging.info('Total number of atoms: %s' % len(self.atoms))
        if self.nBonds > 0:
            logging.info('Total number of bonds: %s' % self.nBonds)
        if self.nAngles > 0:
            logging.info('Total number of angles: %s' % self.nAngles)
        if self.nDihedrals > 0:
            logging.info('Total number of dihedral: %s' % self.nDihedrals)

        logging.info('Atom types: %s' % self.nAtomTypes)
        if self.nBonds > 0:
            logging.info('Bond types: %s' % self.nBondTypes)
        if self.nAngles > 0:
            logging.info('Angle types: %s' % self.nAngleTypes)
        if self.nDihedrals > 0:
            logging.info('Total number of dihedral types: %s' % self.nDihedralTypes)
        x_length = ("(%s - %s)" % (self.xlo,self.xhi))
        y_length = ("(%s - %s)" % (self.ylo,self.yhi))
        z_length = ("(%s - %s)" % (self.ylo,self.zhi))
        logging.info('Simulation box: [%s, %s, %s]' % (x_length,y_length,z_length))

        file.write("Masses\n\n")
        for i in range(self.nAtomTypes):
            file.write("\t%d\t%lf\n" % (i+1,(self.mass[i])))
        file.write("\n")

        file.write("Atoms\n\n")
        for i in range(len(self.atoms)):
            a = self.atoms[i]
            file.write("\t%ld\t%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d\n" % (a.id,a.mol,a.type,a.q,self.atoms[i].x,self.atoms[i].y,self.atoms[i].z,0,0,0))
        file.write("\n")

        file.write("Velocities\n\n")
        for i in range(len(self.atoms)):
            v = self.velocities[i]
            file.write("\t%ld\t%lf\t%lf\t%lf\n" % (v.id,self.velocities[i].vx,self.velocities[i].vy,self.velocities[i].vz))
        file.write("\n")

        if len(self.bonds) > 0:
           file.write("Bonds\n\n")
           for i in range(len(self.bonds)):
               b = self.bonds[i]
               file.write("\t%ld\t%ld\t%ld\t%ld\n" % (i+1,b.type,b.bond1,b.bond2))
           file.write("\n")

        if len(self.angles) > 0:
           file.write("Angles\n\n")
           for i in range(len(self.angles)):
               ang = self.angles[i]
               file.write("\t%ld\t%ld\t%ld\t%ld\t%ld\n" % (i+1,ang.type,ang.angle1,ang.angle2,ang.angle3))
           file.write("\n")

        if len(self.dihedrals) > 0:
           file.write("Dihedrals\n\n")
           for i in range(len(self.dihedrals)):
               d = self.dihedrals[i]
               file.write("\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n" % (i+1,d.type,d.dihedral1,d.dihedral2,d.dihedral3,d.dihedral4))
           file.close()


def get_command_line_args(args):
    if len(args) != 4:
        logging.info('''Error: The number of files is invalid
        \t\t\t Usage: %s <confFile> <dumpFile> <modified confFile> ''' % args[0])
        sys.exit()
    return args

def main():
    args = get_command_line_args(sys.argv)
    confFile = args[1]
    dumpFile = args[2]
    outFile = args[3]
    logging.info('confFile: %s' % confFile)
    logging.info('dumpFile: %s' % dumpFile)
    logging.info('modified confFile: %s' % outFile)
    lmpconf = LmpConf()
    lmpconf.read_configuration_file(confFile)
    dumpData = lmpconf.read_dump_file(dumpFile)
    lmpconf.update_configuration_file(dumpData,confFile)
    lmpconf.write_configuration_file(outFile)

if __name__ == "__main__":
    main()

end = time.time()
logging.info('Total elapsed time: %s seconds' % (end - start))

"""
NOTE: Adjust the Shebang (Header) line according to your OS
Windos: #! python3
Linux: #! /usr/bin/env python3
OSX: #! /usr/bin/python3
"""
