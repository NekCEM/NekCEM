import linefeeder




def extract_value(s):
    return s.split()[0]

def extract_int(s):
    return int(extract_value(s))




class BoxFile:
    def __init__(self, filename=None):
        self.ReaName = ""
        self.Dimensions = 2
        self.Before = []
        self.NElements = [4,4]
        self.After = []

	if filename is not None:
	    self.read(filename)

    def read(self, filename):
        inf = file(filename, "r")
        feeder = linefeeder.LineFeeder(inf.readlines())
        inf.close()

        self.ReaName = feeder.get_line()
        self.Dimensions = extract_int(feeder.get_line())
        self.Before = []
        while not feeder.eof() and feeder.peek().strip().lower() != "box":
            self.Before.append(feeder.get_line())
        self.Before.append(feeder.get_line()) # append box line, too.
        self.NElements = [-int(s) 
                for s in feeder.get_line().split()[:self.Dimensions]]
        for els in self.NElements: assert els > 0
        self.After = feeder.get_rest()

    def write(self, filename):
        outf = file(filename, "w")
        outf.write(self.ReaName)
        outf.write("%d\n" % self.Dimensions)
        outf.write("".join(self.Before))
        outf.write(" ".join(["%d" % -i for i in self.NElements])+" (number of elements in x,y,z)\n")
        outf.write("".join(self.After))
        outf.close()






if __name__ == "__main__":
    import sys
    f = BoxFile(sys.argv[1])
    print f.NElements
    f.write("blah.temp")

