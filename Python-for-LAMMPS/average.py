#! python3

"""
This program averages LAMMPS' ave/chunk output files
"""
import sys, logging, time
start = time.time()

FORMAT = '%(asctime)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.DEBUG,format = FORMAT)

"""
Leave the following line uncommented to avoid Info. messages during the run
"""
#logging.disable(logging.CRITICAL)

class Profile:
    def __init__(self,*args):
        attributes = args[0]
        attribute = 0
        self.list = []
        while attribute < len(attributes):
            self.list.append(attributes[attribute])
            attribute = attribute + 1
        global number_of_variables
        number_of_variables = attribute

class LmpConf:
    def __init__(self):
        self.step = None
        self.nchunks = None
        self.tchunks = None
        self.chunk = []

    def read_configuration_file(self,confFile,document):
        try:
            file = open(confFile,"r")
            line = file.readline()
            line = file.readline()
            line = file.readline()
            if document == 0:
                chunk_coor_ncount =  str(line).strip("#\n")
                logging.info('The variables in the files: %s' % chunk_coor_ncount )
            line = file.readline()
            try:
                self.step = int(line.strip().split()[0])
                self.nchunks = int(line.strip().split()[1])
                self.tchunks = int(line.strip().split()[2])
            except:
                logging.info('Error: File %s have invalid content' % (document+1))
                sys.exit()
            list_file = []
            while line:
                for i in range(self.nchunks):
                    l = file.readline().strip().split()
                    list = []
                    for j in range (len(l)):
                        list.append(0)
                    try:
                        for k in range (len(l)):
                            if k == 0:
                                list[k] = int(l[0])
                            else:
                                list[k] = float(l[k])
                    except:
                        logging.info("Error: File %s has invalid text format this file is number %s on the input order"  % (confFile, document+1))
                        sys.exit()

                    if len(list) == 0:
                        logging.info('Error: In the file %s one line is empty' % (confFile))
                        sys.exit()
                    list_file.append(Profile(list))
                line = file.readline()
            self.chunk.append(list_file)

        except:
            logging.info("Error: Fail to correctly read File %s this file is number %s on the input order"  % (confFile, document+1))
            sys.exit()
        file.close()

    def write_configuration_file(self,outFile):
        file = open(outFile,"w")
        file.write("\t%ld\t%ld\t%ld\n" % (self.step,self.nchunks,self.tchunks))
        line_file = []
        logging.info('Total number of variables: %s' % number_of_variables)
        logging.info('Total number of bins (chunks): %s' % self.nchunks)
        for nchunk in range(self.nchunks):
            line_file.append([0]*number_of_variables)
        for document in range(documents):
            list = self.chunk[document]
            for chunk in range(self.nchunks):
                a = list[chunk]
                alist = a.list
                variables = len(alist)
                line = []
                for variable in range(variables):
                    line.append(alist[variable])
                    try:
                        line_file[chunk][variable] = line_file[chunk][variable] + alist[variable]
                    except:
                        logging.info('''Error: File number %s has %s or more variables in one of its lines
                        \t Note that the number of variables is set to %s according to file number %s
                        ''' % (document+1,variable+1,number_of_variables,documents))
                        sys.exit()
        for a in range(self.nchunks):
            variable = 0
            for variable in range(variables):
                if variable == 0:
                    line_file[a][variable] = ("%ld" % (line_file[a][variable]/documents))
                elif variable == 1 or variable == 2 :
                    line_file[a][variable] = ("%lf" % (line_file[a][variable]/documents))
                else:
                    line_file[a][variable] = ("%lE" % (line_file[a][variable]/documents))

            writing_line = line_file[a]
            file.writelines(' '.join(writing_line))
            file.write("\n")
        file.close()

def get_command_line_args(args):
    if len(args) <= 2:
        logging.info('''Error: The number of files is invalid
        \t\t\t Usage: %s <LAMMPS reading file(s)> <output file> ''' % args[0])
        sys.exit(1)
    logging.info('The number of files: %s' % (len(args)-2))
    return args

def main():
    args = get_command_line_args(sys.argv)
    global documents
    documents = len(args)-2
    lmpconf = LmpConf()
    for document in range (documents):
        lmpconf.read_configuration_file(args[document+1],document)
    lmpconf.write_configuration_file(args[-1])

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
