#Python script to convert temp/chunk points to a grid in ttm simulations
#This script assumes that nx, ny, and nz are equal
#This script prints the max and average temperature values of the ionic temperature estimator in each provided timestep

import sys
import os
import statistics

class Grid:

    def __init__(self):
        self.separated_files = []

    def separate(self,input_file,root_name):
        infile = open(input_file,"r")
        line = infile.readline()

        nodes=0

        while line:
            if line.startswith("#"):
                pass
            else:
                l = line.strip().split()
                if l[0]=="0":
                    nodes=l[1]
                if l[1]==nodes:
                    timestep=l[0]
                    outputfile = root_name[0]+"_"+timestep+"_grid_"+root_name[1]
                    outfile = open(outputfile,"w")
                    outfile.write("%s %s\n" % (timestep,nodes))
                    self.separated_files.append(outputfile)
                    for i in range(int(nodes)):
                        chunk = infile.readline().strip().split()
                        outfile.write("%s %s\n" % (chunk[0],chunk[1]))
            line = infile.readline()
        infile.close()

    def write_nodes(self,separated_file,root_name):
        separated_file = open(separated_file,"r")
        line = separated_file.readline()
        header = line.strip().split()
        timestep=header[0]
        nodes=header[1]
        node=round(int(nodes)**(1/3)) #Assuming cubic ionic grid
        output_file=root_name[0]+"_"+timestep+root_name[1]
        file = open(output_file,"w")
        file.write("#Ionic temperature on %sx%sx%s grid at step %s\n" % (node,node,node,timestep))
        x,y,z= 1,1,1
        temperatures=[]
        line = separated_file.readline()
        while line:
            l = line.strip().split()
            write_line = str(x)+" "+str(y)+" "+str(z)+" "+str(l[1])+"\n"
            temperatures.append(float(l[1]))
            file.write(write_line)
            x += 1
            if x > 10:
                x=1
                y+=1
            if y > 10:
                y=1
                z+=1
            line = separated_file.readline()
        separated_file.close()
        temperatures.sort()
        print(timestep,temperatures[-1], statistics.mean(temperatures), min(temperatures))

def get_command_line_args(args):
    if len(args) != 2:
        print('''Error: The number of files is invalid
        Usage: %s <CIT input file>''' % args[0])
        print('bad number of files: %s' % (len(args)))
        sys.exit(1)
    return args

def main():
    args = get_command_line_args(sys.argv)
    input_file = args[1]
    root_name = os.path.splitext(input_file)
    grid = Grid()
    grid.separate(input_file,root_name)
    print("#timestep maxT avgT minT")
    for separated_file in grid.separated_files:
        grid.write_nodes(separated_file,root_name)


if __name__ == "__main__":
    main()
