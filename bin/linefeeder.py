class LineFeeder:
    def __init__(self, lines):
	self.Lines = lines
	self.Index = 0

    def get_lines(self, count):
	result = self.Lines[self.Index:self.Index+count]
	if self.Index + count > len(self.Lines):
	    raise RuntimeError, "requesting %d lines at line %d -- only have %d" % \
		    (count, self.Index+1, len(self.Lines)-self.Index)
	self.Index += count
	return result

    def get_line(self):
	return self.get_lines(1)[0]

    def get_rest(self):
	result = self.Lines[self.Index:]
	self.Index = len(self.Lines)
	return result

    def peek(self):
	return self.Lines[self.Index]

    def eof(self):
        return self.Index == len(self.Lines)
