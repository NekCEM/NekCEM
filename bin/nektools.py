from data_tags import *




def find_root():
    from os.path import join, isfile, abspath

    curdir = "."
    for i in range(10):
	if isfile(curdir+"/src/cem_dg.F"):
	    return abspath(curdir)
	curdir = join(curdir, "..")

    raise RuntimeError, "root of tree not found"




def parse_fortran_float(s):
    # Awful hack. Why can't Fortran be more like the other kids?
    try:
        return float(s.replace("D","E"))
    except ValueError:
        import re
        float_re = re.compile("^([-+]?[0-9]*\.[0-9]+)[DEde]?([+-][0-9]+)$")
        match = float_re.match(s.strip())
        if match is None:
            raise ValueError, "Not even the screwy float handling could deal with '%s'" % s
        return float(match.group(1)) * 10.** int(match.group(2))




def sort_by(list, key_func, ascending=True):
    list.sort(lambda x,y: cmp(key_func(x), key_func(y)))
    if not ascending:
        list.reverse()




def my_enumerate(seq):
    return zip(range(len(seq)), seq)



class MySet:
    """Poor man's sets.Set implementation for Python 2.2 compatiblity.
    """
    def __init__(self, items=[]):
        self.Dict = {}
        for i in items:
            self.Dict[i] = 1

    def __len__(self): return len(self.Dict)
    def __contains__(self, i): return i in self.Dict
    def __iter__(self): return self.Dict.iterkeys()
    def __or__(self, other): 
        result = self.copy()
        result.update(other)
        return result
    def add(self, i): self.Dict[i] = 1
    def copy(self):
        return MySet(self.Dict.keys())
    def update(self, iterable):
        for i in iterable:
            self.Dict[i] = 1



def short_threshold(x):
    # VTK stores data in single precision internally and it reports 
    # errors on numbers smaller than 1e-38 (roughly) rather
    # than doing the reasonable thing and rounding to zero.
    # This routine does the necessary chopping to avoid errors.
    if abs(x) < 1e-30:
	return 0.
    elif abs(x) > 1e+30:
	return x/abs(x)*1e+30
    else:
	return x




def parse_point(line):
    return tuple([
	    short_threshold(parse_fortran_float(i))
	    for i in line.split()])




def add_tuples(ta, tb):
    return tuple([a+b for a,b in zip(ta, tb)])





class VTKDumpGeometry:
    def __init__(self, filename):
	print "geometry..."
	geo_file = file(filename, 'r').readlines()
	first_line_fields = geo_file[0].split()

	self.NElements, \
		self.PointsPerElX, \
		self.PointsPerElY, \
		self.PointsPerElZ = [int(i) for i in first_line_fields]
	self.Points = [parse_point(line) for line in geo_file[1:]]
	assert len(self.Points) == self.PointCount

	self.PointIndexLookup = self.generate_index_lookup()

	print "fine grid..."
        if self.is2d():
            self.Volumes = list(self.generate_quads())
        else:
            self.Volumes = list(self.generate_hexes())

	print "faces..."
        self.Faces = list(self.generate_faces())

	print "cell polygons..."
	self.Cells = list(self.generate_cells())

    def is2d(self):
        return self.PointsPerElZ == 1

    def volumes(self):
        return self.Volumes

    def cells(self):
        return self.Cells

    def generate_index_lookup(self):
	index = 0
	pil = {}
	for e in range(self.NElements):
	    for k in range(self.PointsPerElZ):
		for j in range(self.PointsPerElY):
		    for i in range(self.PointsPerElX):
			pil[i,j,k,e] = index
			index += 1
        return pil

    def generate_quads(self):
	pil = self.PointIndexLookup
	for e in range(self.NElements):
            for j in range(self.PointsPerElY-1):
                for i in range(self.PointsPerElX-1):
                    yield (pil[i  ,j  ,0  ,e],
                            pil[i+1,j  ,0  ,e],
                            pil[i+1,j+1,0  ,e],
                            pil[i  ,j+1,0  ,e],
			    )

    def generate_hexes(self):
	pil = self.PointIndexLookup
	for e in range(self.NElements):
	    for k in range(self.PointsPerElZ-1):
		for j in range(self.PointsPerElY-1):
		    for i in range(self.PointsPerElX-1):
			yield (
			    pil[i  ,j  ,k  ,e],
			    pil[i+1,j  ,k  ,e],
			    pil[i+1,j+1,k  ,e],
			    pil[i  ,j+1,k  ,e],
			    pil[i  ,j  ,k+1,e],
			    pil[i+1,j  ,k+1,e],
			    pil[i+1,j+1,k+1,e],
			    pil[i  ,j+1,k+1,e],
			    )

    def generate_cells(self):
	def get_manhattan_line((i1, j1, k1), (i2, j2, k2), el):
	    pil = self.PointIndexLookup
	    while i1 != i2 or j1 != j2 or k1 != k2:
		yield pil[i1, j1, k1, el]
		if i1 < i2: i1 += 1
		if j1 < j2: j1 += 1
		if k1 < k2: k1 += 1
		if i1 > i2: i1 -= 1
		if j1 > j2: j1 -= 1
		if k1 > k2: k1 -= 1

	def get_side(base, axis1, axis2, el):
	    a = base
	    b = add_tuples(base, axis1)
	    c = add_tuples(add_tuples(base, axis1), axis2)
	    d = add_tuples(base, axis2)
	    return list(get_manhattan_line(a, b, el))+\
		    list(get_manhattan_line(b, c, el))+\
		    list(get_manhattan_line(c, d, el))+\
		    list(get_manhattan_line(d, a, el))

	nx, ny, nz = self.PointsPerElX-1, self.PointsPerElY-1, self.PointsPerElZ-1
	for e in range(self.NElements):
	    yield get_side((0,0,0), (nx,0,0), (0,ny,0), e)
            if not self.is2d():
                yield get_side((0,0,nz), (nx,0,0), (0,ny,0), e)
                yield get_side((0,0,0), (0,ny,0), (0,0,nz), e)
                yield get_side((nx,0,0), (0,ny,0), (0,0,nz), e)
                yield get_side((0,0,0), (nx,0,0), (0,0,nz), e)
                yield get_side((0,ny,0), (nx,0,0), (0,0,nz), e)

    def generate_faces(self):
	pil = self.PointIndexLookup

        class Face:
            def __init__(self, element, face_idx, point_indices, points):
                self.Element = element
                self.FaceIndex = face_idx
                self.PointIndices = point_indices
                self.Points = points

	def get_face_3d(el, face_idx, base, axis1, axis2):
            e_tup = (el,)
            a = pil[base+e_tup]
            b = pil[add_tuples(base, axis1)+e_tup]
            c = pil[add_tuples(add_tuples(base, axis1), axis2)+e_tup]
            d = pil[add_tuples(base, axis2)+e_tup]
            point_indices = [a, b, c, d]
            points = [self.Points[pi] for pi in point_indices]
            return Face(el, face_idx, point_indices, points)

        if self.is2d():
            return
        else:
            nx, ny, nz = self.PointsPerElX-1, self.PointsPerElY-1, self.PointsPerElZ-1
            for e in range(self.NElements):
                # ordering consistent with pff/sym notation (see connect1.f)
                yield get_face_3d(e, 1, (0 ,0,0), (0,ny,0), (0,0,nz))
                yield get_face_3d(e, 2, (nx,0,0), (0,ny,0), (0,0,nz))

                yield get_face_3d(e, 3, (0,0 ,0), (nx,0,0), (0,0,nz))
                yield get_face_3d(e, 4, (0,ny,0), (nx,0,0), (0,0,nz))

                yield get_face_3d(e, 5, (0,0,0 ), (nx,0,0), (0,ny,0))
                yield get_face_3d(e, 6, (0,0,nz), (nx,0,0), (0,ny,0))

    def point_count(self):
	return self.NElements*self.PointsPerElX*self.PointsPerElY*self.PointsPerElZ
    PointCount = property(point_count)




# Run log data ---------------------------------------------------
class TaggedDataTable(object):
    """This is a miniature database table. Each value can have 
    multiple `tags', by which the table can be filtered.

    Below, it is used for error values.
    """
    def __init__(self, dtable = []):
        self.DataTable = dtable[:]

    def union(self, other):
        return TaggedDataTable(self.DataTable + other.DataTable)

    def __len__(self):
        return len(self.DataTable)

    def __add__(self, other):
        return self.union(other)

    def add_row(self, value, tags):
        if not isinstance(tags, dict):
            raise ValueError, "expeceted dictionary as `tags' argument"
        self.DataTable.append((tags, value))

    def add_dependent_tag(self, tag, f):
        for row_tags, row_value in self.DataTable:
            row_tags[tag] = f(row_tags)

    def tag_values(self, tag):
        result = MySet()
        for row_tags, row_value in self.DataTable:
            result.add(row_tags.get(tag, None))
        return list(result)

    def filter(self, tag, tag_value):
        result = []
        for row_tags, row_value in self.DataTable:
            if row_tags.get(tag, None) == tag_value:
                new_row_tags = row_tags.copy()
                del new_row_tags[tag]
                result.append((new_row_tags, row_value))
        return TaggedDataTable(result)

    def iteritems(self, key_tag):
        return [(row_tags[key_tag], row_value) 
          for row_tags, row_value in self.DataTable]

    def values(self):
        return [row_value \
                for row_tags, row_value in self.DataTable]

    def values_sorted_by(self, sort_tag):
        sorted_values = self.DataTable[:]
        sort_by(sorted_values, lambda x: x[0].get(sort_tag, None))
        return [row_value \
                for row_tags, row_value in sorted_values]





def parse_run_log(edat, usr_stem=None, pdegree=None, rea_stem=None, 
        box_resolution=None, dt=None, rea_switches={}, rea_parameters={}):
    import math
    lines = edat.split("\n")[1:]
    errors = TaggedDataTable()
    times = TaggedDataTable()

    common_tags = {}
    if usr_stem is not None: common_tags[TAG_USRSTEM] = usr_stem
    if rea_stem is not None: common_tags[TAG_REASTEM] = rea_stem
    if dt is not None: common_tags[TAG_DT] = dt
    if pdegree is not None: common_tags[TAG_PDEGREE] = pdegree

    if box_resolution is not None:
        common_tags[TAG_RESOLUTION] = box_resolution

    for rsk, rsv in rea_switches.iteritems():
        common_tags[rsk] = rsv

    for rpk, rpv in rea_parameters.iteritems():
        common_tags[rpk] = rpv

    for line in lines:
        if len(line.strip()) == 0:
            continue
        words = line.split()

        if words[0] == "errors":
            errors_here = TaggedDataTable()
            errtype = words[1]
            if len(errtype) > 4:
                # HACK!
                if errtype.startswith("L2"):
                    words[2:2] = [errtype[2:]]
                    errtype = "L2"
                elif errtype.startswith("Linf"):
                    words[2:2] = [errtype[4:]]
                    errtype = "Linf"

            try:
                step = int(words[2])
            except ValueError:
                print "*** Warning: invalid step number: %s" % words[2]
                continue

            idx = 3

            while idx < len(words):
                assert words[idx].endswith(":")
                component = words[idx][:-1]
                error = parse_fortran_float(words[idx+1])

                etags = {TAG_TYPE: errtype,
                        TAG_STEP: step,
                        TAG_COMPONENT: component,
                        }
                etags.update(common_tags)

                errors_here.add_row(error, etags)

                idx += 2
                
            etags = {TAG_TYPE: errtype,
                    TAG_STEP: step,
                    TAG_COMPONENT: "sum",
                    }
            etags.update(common_tags)

            import operator 
            errors_here.add_row(math.sqrt(reduce(operator.add,
                [x**2 for x in errors_here.values()])),
                etags)

            errors += errors_here
        elif words[0] == "time":
            try:
                step = int(words[1])
            except ValueError:
                print "*** Warning: invalid step number: %s" % words[2]
                continue

            ttags = {TAG_TYPE: "wall",
                    TAG_STEP: step,
                    TAG_COMPONENT: "sum",
                    }


            ttags.update(common_tags)
            times.add_row(parse_fortran_float(words[2]), ttags)
        else:
            raise ValueError, "unkown run log line tag '%s'" % words[0]

    return errors, times

