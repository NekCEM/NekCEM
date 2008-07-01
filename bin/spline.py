from __future__ import division



def vandermonde(vector, degree=None):
    """Return a Vandermonde matrix, i.e.

    1 v1 v1^2 ...
    1 v2 v2^2 ...
    ...
    """
    import Numeric as num

    if degree is None:
        degree = len(vector)-1

    mat = num.zeros((len(vector), degree+1), num.Float64)
    for i, v in enumerate(vector):
        for power in range(degree+1):
            mat[i,power] = v**power
    return mat




def fit_polynomial(x_vector, data_vector):
    """Fit a polynomial of degree `degree' at the points 
    in `x_vector' to the `data_vector' by interpolation.
    """

    import Numeric as num
    import LinearAlgebra as la

    degree = len(x_vector)-1

    vdm = vandermonde(x_vector, degree)
    result = la.solve_linear_equations(vdm, num.array(data_vector))
    result = list(result)
    result.reverse()
    return Polynomial(result)




class Polynomial:
    def __init__(self, coefficients):
        # highest-order first
        self.Coefficients = coefficients

    def degree(self):
        return len(self.Coefficients)-1

    def __call__(self, x):
        result = 0.
        for c in self.Coefficients:
            result = result * x + c
        return result

    def as_fortran_contlines(self, var):
        cont = "     $  "
        result = cont \
                +self.degree()*"(" \
                +repr(self.Coefficients[0])+"\n"
        for c in self.Coefficients[1:]:
            result += "%s)*%s+%s\n" % (cont, var, repr(c))
        return result




class Slab:
    def __init__(self, start, end, polynomials, attributes):
        self.Start = start
        self.End = end
        self.Polynomials = polynomials # a list of Polynomial instances
        self.Attributes = attributes # a list of strings




class Spline:
    def __init__(self, filename=None):
        self.clear()
        if filename is not None:
            self.read(filename)

    def clear(self):
        self.Slabs = []

    def read(self, filename):
        import re
        from linefeeder import LineFeeder

        lines = file(filename, "r").readlines()

        # eliminate comments and blank lines
        mylines = []
        for line in lines:
            com_start = line.find("#")
            if com_start != -1:
                line = line[:com_start]
            line = line.strip()
            if len(line):
                mylines.append(line)

        feeder = LineFeeder(mylines)

        def get_grouphead():
            line = feeder.peek()
            match = re.match(r"^\[(\w+)\]$", line)
            if match is None:
                raise RuntimeError, "group head expected, '%s' found" % line
            feeder.get_line()
            return match.group(1)

        def is_grouphead():
            line = feeder.peek()
            match = re.match(r"^\[(\w+)\]$", line)
            return match is not None

        def parse_keyvalue(line):
            reg_match = re.match(r"^(\w+)\=(.*)$", line)
            subscr_match = re.match(r"^(\w+)\[([a-zA-Z0-9,.]*)\]\=(.*)$", line)
            if subscr_match is not None:
                return subscr_match.group(1), subscr_match.group(2), \
                        subscr_match.group(3)
            elif reg_match is not None:
                return reg_match.group(1), None, reg_match.group(2)
            else:
                raise RuntimeError, "key-value pair expected, '%s' found" % line

        def get_group():
            head = get_grouphead()
            group_lines = []
            
            while not feeder.eof() and not is_grouphead():
                group_lines.append(feeder.get_line())
            return head, group_lines

        def get_group_line_dict(group_lines):
            result = {}
            for line in group_lines:
                name, subscript, value = parse_keyvalue(line)
                if subscript is not None:
                    result[name,subscript] = value
                else:
                    result[name] = value
            return result

        def parse_file():
            constants = {}
            composition = []
            current_section = []
            sections = {}

            while not feeder.eof():
                head, group_lines = get_group()
                if head == "newsection":
                    name = get_group_line_dict(group_lines)["name"]
                    current_section = []
                    sections[name] = current_section
                elif head == "composition":
                    composition = group_lines
                elif head == "constants":
                    for line in group_lines:
                        name, subscript, value = parse_keyvalue(line)
                        if subscript is not None:
                            raise RuntimeError, "no subscripts allowed in constant groups"
                        constants[name] = eval(value, constants)
                elif head == "slab":
                    length = None
                    values = {}
                    attributes = {}
                    gl_parsed = [parse_keyvalue(line) for line in group_lines]

                    for name, subscript, value in gl_parsed:
                        if name == "length":
                            length = eval(value, constants)
                        if name == "attribute":
                            attributes[subscript] = value
                        elif name == "value":
                            evaluated = eval(value, constants)
                            if subscript is None:
                                raise RuntimeError, "subscript required in value"
                            subelements = subscript.split(",")
                            var = subelements[0]
                            if len(subelements) == 1:
                                values.setdefault(var, []).append((None, evaluated))
                            elif len(subelements) == 2:
                                parm_point = eval(subelements[1], constants)
                                values.setdefault(var, []).append((parm_point, evaluated))
                            else:
                                raise RuntimeError, "invalid number of subscripts in '%s'" % subscript

                    if length is None:
                        raise RuntimeError, "must specify a length"

                    current_section.append((length, values, attributes))
                elif head == "settings":
                    pass
                else:
                    raise RuntimeError, "invalid group '%s'" % head
            return composition, sections

        def compose_spline(composition, sections):
            def fill_abscissae(abscissae, desired_length, end_value):
                if len(abscissae) == 0:
                    abscissae.append(0.)
                    last_absc = 0.
                else:
                    last_absc = abscissae[-1]

                new_value_count = desired_length - len(abscissae)
                if new_value_count == 0:
                    return

                slope = (end_value-last_absc)/new_value_count
                for i in range(new_value_count):
                    abscissae.append((i+1)*slope+last_absc)

            position = 0
            slabs = []

            for section_name in composition:
                inverted = section_name.startswith("-")
                if inverted:
                    section_name = section_name[1:]
                section = sections[section_name]
                if inverted:
                    section = section[:]
                    section.reverse()

                for slab_length, slab_values, slab_attributes in section:
                    slab_start = position
                    position = slab_end = position + slab_length

                    polys = {}
                    for var, values in slab_values.iteritems():
                        abscissae = []
                        ordinates = []
                        for abscissa, ordinate in values:
                            ordinates.append(ordinate)
                            if abscissa is not None:
                                fill_abscissae(abscissae, len(ordinates), 
                                        abscissa)

                        fill_abscissae(abscissae, len(ordinates), 1.)
                        if inverted:
                            abscissae = [1-x for x in abscissae]
                        abs_map_slope = slab_end-slab_start
                        abscissae = [abs_map_slope*x+slab_start for x in abscissae]
                        poly = fit_polynomial(abscissae, ordinates)
                        for x, y in zip(abscissae, ordinates):
                            assert abs(poly(x)-y) < 1e-5
                        polys[var] = poly

                    slabs.append(Slab(slab_start, slab_end, polys,
                        slab_attributes))
            return slabs

        composition, sections = parse_file()
        self.Slabs = compose_spline(composition, sections)

    def as_fortran_func(self, name, var):
        body = ""
        if_prefix = ""
        for s in self.Slabs:
            body += \
                    "      %sif ((x.ge.%g).and.(x.le.%g)) then\n" \
                    "        %s =\n%s" % (if_prefix, 
                            s.Start, s.End,
                            name,
                            s.Polynomials[var].as_fortran_contlines("x"))
            if_prefix = "else"

        return \
                "      real function %s(x)\n" \
                "      real x\n\n" \
                "%s" \
                "      elseif (x.le.0) then\n" \
                "        %s = %d\n" \
                "      elseif (x.ge.%g) then\n" \
                "        %s = %d\n" \
                "      else\n" \
                "        write(*,*) '%s: how did i get here?'\n" \
                "        stop\n" \
                "      endif\n\n" \
                "      return\n" \
                "      end\n" \
                % (name, body, 
                        name, self.Slabs[0].Polynomials[var](0),
                        self.end(), name, self.Slabs[-1].Polynomials[var](self.end()),
                        name)

    def end(self):
        if not self.Slabs:
            return 0
        else:
            return self.Slabs[-1].End

    def __call__(self, var, x):
        for s in self.Slabs:
            if s.Start <= x <= s.End:
                return s.Polynomials[var](x)

        raise ValueError, "%f out of range for spline" % x

    def plot(self, vars):
        end = self.end()
        n = 1000
        xpts = [end/n*i for i in range(n+1)]

        import Gnuplot
        gp = Gnuplot.Gnuplot()
        datasets = []
        for var in vars:
            ypts = [self(var, x) for x in xpts]
            datasets.append(Gnuplot.Data(xpts, ypts, with="lines",
                title=var))
        gp.plot(*datasets)
        raw_input("Hit Enter: ")


if __name__ == "__main__":
    import sys
    s = Spline(sys.argv[1])
    s.plot(sys.argv[2])
    #print s.as_fortran_func("f", sys.argv[2])
