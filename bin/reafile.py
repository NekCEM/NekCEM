import linefeeder
import nektools




def extract_value(s):
    return s.split()[0]

def extract_comment(s):
    # HACK!
    return " ".join(s.split()[2:])

def extract_int(s):
    return int(extract_value(s))

def extract_float(s):
    return float(extract_value(s))

def extract_bool(s):
    return extract_value(s).lower() == "t"

def build_line(value, index, comment):
    value = str(value)
    pad = (30-len(value))*" "
    # DO NOT eliminate leading whitespace
    # or else n2to3 will get very mad.
    return "  %s%s%d: %s\n" % (value, pad, index, comment)

class REALineFeeder(linefeeder.LineFeeder):
    def get_int(self):
	return extract_int(self.get_line())

    def get_float(self):
	return extract_float(self.get_line())

class REAFile:
    def __init__(self, filename=None):
	self.Header = ["", "", "", "0 PARAMETERS FOLLOW"]
	self.Dimensions = 0
	self.Parameters = []
	self.PassiveScalars = []
	self.Switches = []
	self.Rest = []

	if filename is not None:
	    self.read(filename)
	else:
	    self.build_name_maps()

    def copy(self):
        result = REAFile()
        result.Header = self.Header[:]
        result.Dimensions = self.Dimensions
        result.Parameters = self.Parameters[:]
        result.PassiveScalars = self.PassiveScalars[:]
        result.Switches = self.Switches[:]
        result.Rest = self.Rest[:]
        result.build_name_maps()
        return result

    def read(self, filename):
	feeder = REALineFeeder(file(filename, "r").readlines())

	self.Header = feeder.get_lines(2)
	self.Dimensions = feeder.get_int()
	parameter_count = feeder.get_int()
	self.Parameters = feeder.get_lines(parameter_count)

	passive_header = feeder.get_line()
	passive_count = extract_int(passive_header)
	self.PassiveScalars = [passive_header] + feeder.get_lines(passive_count)

	switch_count = feeder.get_int()
	self.Switches = feeder.get_lines(switch_count)

	self.Rest = feeder.get_rest()

	self.build_name_maps()

    def write(self, filename):
	outf = file(filename, "w")
	outf.write("".join(self.Header))
	outf.write("%d DIMENSIONAL RUN\n" % self.Dimensions)
	outf.write("%d PARAMETERS FOLLOW\n" % len(self.Parameters))
	outf.write("".join(self.Parameters))
	outf.write("".join(self.PassiveScalars))
	outf.write("%d LOGICAL SWITCHES FOLLOW\n" % len(self.Switches))
	outf.write("".join(self.Switches))

	outf.write("".join(self.Rest))

    def build_name_maps(self):
	def build_map(lines):
	    import re

	    ambiguous = nektools.MySet()
	    map = {}
	    for i, line in nektools.my_enumerate(lines):
		try:
		    name = line.split()[2].lower()
		    if re.match("^[A-Za-z0-9]+$", name):
			if name in map:
			    ambiguous.add(name)
			elif name not in ambiguous:
			    map[name] = i+1
		except IndexError:
		    pass
	    return map

	self.ParameterMap = build_map(self.Parameters)
	self.SwitchMap = build_map(self.Switches)

    def get_line(self, lines, map, which):
	if isinstance(which, str):
	    which = map[which.lower()]
	return lines[which-1]

    def set_value(self, lines, map, which, value):
	if isinstance(which, str):
	    which = map[which.lower()]
	lines[which-1] = build_line(value, which, 
		extract_comment(lines[which-1]))

    def get_switch(self, which):
	return extract_bool(self.get_line(self.Switches, self.SwitchMap, which))

    def set_switch(self, which, value):
        if value == True:
            value = "T"
        elif value == False:
            value = "F"

	return self.set_value(self.Switches, self.SwitchMap, which, value)

    def get_parameter(self, which):
	return extract_float(self.get_line(self.Parameters, self.ParameterMap, which))

    def set_parameter(self, which, value):
	value = str(value)
	return self.set_value(self.Parameters, self.ParameterMap, which, value)



if __name__ == "__main__":
    f = REAFile("../examples/cylwave/cylwave.rea")
    print f.get_switch("IFUPWIND")
    f.set_switch("IFUPWIND", False)
    print f.get_parameter("TMODTYP")
    f.set_parameter("tmodtyp", 17)
    f.write("demo.rea")
