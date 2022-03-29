import numpy as np


class Point:

    def __init__(self, x: float, y: float):
        self.x = x
        self.y = y


    def __repr__(self):
        return f'[{self.x}, {self.y}]'


    def to_np(self):
        return np.asarray([self.x, self.y])


    def segment(self, direction):
        '''
        Return directed segment, which has the direction of @direction, formed by starting at @self 
        that has length of 1

        @direction: float, in [0, 2 * pi)

        Compute the coordinate of the end point b by forming right triangle abc with right angle 
        at vertex c
        '''
        b = None
        acute_angle = None
        signs = None

        if direction <= np.pi / 2:
            signs = [1, 1]
            acute_angle = direction
        elif direction <= np.pi:
            signs = [-1, 1]
            acute_angle = np.pi - direction
        elif direction <= 3 * np.pi / 2:
            signs = [-1, -1]
            acute_angle = direction - np.pi
        else:
            signs = [1, -1]
            acute_angle = 2 * np.pi - direction

        bc_length = np.sin(acute_angle)
        ac_length = np.cos(acute_angle)
        b = Point(self.x + signs[0] * ac_length, self.y + signs[1] * bc_length)

        return Segment(self, b)



class Segment:
    '''directed straight line'''

    def __init__(self, start_point, end_point):
        '''
        @start_point, @end_point: Point
        '''
        self.start_point = start_point
        self.end_point = end_point
        self.line = self.line_formula()


    def direction(self):
        '''
        Return angle between right oriented horizontal ray (0, 1) and @self, in [0, 2 * pi)
        Using a.b = |a||b|cos(theta)
        '''
        s = np.array([self.end_point.x - self.start_point.x, self.end_point.y - self.start_point.y])

        theta = np.arctan2(s[1], s[0])
        
        return (theta / np.pi % 2) * np.pi


    def angle(self, t):
        '''
        Return the angle between @self and @t

        @t: Segment
        '''
        return t.direction() - self.direction()


    def perpendicular(self, p):
        '''
        Return the segment pq perpendicular to @self with source @p and target on line(@self)

        @p: Point

        line(@self = [a, b]): (a.y - b.y)x - (a.x - b.x)y = a.y * b.x - a.x * b.y
        q in line(@self) -> (a.y - b.y)q.x - (a.x - b.x)q.y = a.y * b.x - a.x * b.y
        pq perpendicular to @self -> [p.x - q.x, p.y - q.y].[a.x - b.x, a.y - b.y] = 0
                                or (a.x - b.x)q.x + (a.y - b.y)q.y = (a.x - b.x)p.x + (a.y - b.y)p.y
        '''
        a = self.start_point
        b = self.end_point

        coefficient_matrix = [[self.line[0], self.line[1]], [b.x - a.x, b.y - a.y]]
        ordinate = [-self.line[2], (b.x - a.x) * p.x + (b.y - a.y) * p.y]
        q = np.linalg.solve(coefficient_matrix, ordinate)
        q = Point(q[0], q[1])

        return Segment(p, q)


    def length(self):
        '''
        Return the length of @self
        '''
        return np.sqrt((self.start_point.x - self.end_point.x) ** 2 + 
            (self.start_point.y - self.end_point.y) ** 2)


    def line_formula(self):
        '''
        Return formula of the line passing @self

        line: mx + ny + p = 0
        '''
        a = self.start_point
        b = self.end_point

        m = a.y - b.y
        n = b.x - a.x
        p = a.x * b.y - a.y * b.x 

        return (m, n, p)


    def intersection(self, g):
        '''
        Return the intersector of @self and @g

        @g: Segment
        '''

        intersect = not (abs(self.direction() - g.direction()) < 1e-7 
            or abs(abs(self.direction() - g.direction()) - np.pi) < 1e-7)
        intersector = None

        if intersect:
            coefficient_matrix = [[self.line[0], self.line[1]], [g.line[0], g.line[1]]]
            ordinate = [-self.line[2], -g.line[2]]

            intersector = np.linalg.solve(coefficient_matrix, ordinate)
            intersector = Point(intersector[0], intersector[1])

        return intersector


    def __repr__(self):
        return f'[{self.start_point}, {self.end_point}]'


class Evtype:

    def __init__(self, v, e):
        '''
        @v: Point
        @e: Segement
        '''

        self.v = v # Point
        self.e = e # Segment


    def edge(self):
        return self.e


    def vertex(self):
        return self.v


    def angle(self, p):
        '''
        Return the angle by the edges of 2 antipodal edge-vertex pairs

        @p: Evtype
        '''
        return self.e.angle(p.e)


    def width(self):
        '''
        Return the width of an antipodal pair (e, v), which is the distance of v
        to its orthogonal projection on the supporting line of e
        '''
        return self.e.perpendicular(self.v).length()


    def __repr__(self):
        return f'{self.e}, {self.v}'


class Parallelogram:

    def __init__(self, z1, z2):
        '''
        @z1, @z2: Evtype
        '''
        self.a = None   # Point
        self.b = None   # Point
        self.c = None   # Point
        z1_direction = z1.e.direction()
        z2_direction = z2.e.direction()
        self.a = z1.e.intersection(z2.e)
        self.b = z1.e.intersection(z2.v.segment(z2_direction))
        self.c = z2.v.segment(z2_direction).intersection(z1.v.segment(z1_direction))


    def drawable(self):
        return not (self.a == None or self.b == None or self.c == None)


    def d(self):
        '''
        Return fourth point from 3 other pts of the parallelogram
        The order in clockwise: a - b - c - d
        -> od = oa - bc
        '''
        return Point(self.a.x - self.b.x + self.c.x, self.a.y - self.b.y + self.c.y)


    def angle(self):
        return Segment(self.a, self.b).angle(Segment(self.b, self.c))


    def area(self):
        theta = abs(self.angle()) % np.pi
        return Segment(self.a, self.b).length() * Segment(self.b, self.c).length() \
                * np.sin(theta)


    def __repr__(self):
        return f'{self.a}, {self.b}, {self.c}, {self.d()}'


def antipodal_pairs(vertices):
    '''
    Traverse through every pair of adjacent vertices, find the point which combines to
    the current segment to create an antipodal pair.

    @vertices: list[Point]
        list of vertices of the polygon
    '''
    antipodal_evs = []

    for i in range(len(vertices)):
        v1 = vertices[i]
        v2 = vertices[(i+1) % len(vertices)]

        max_distance = 0
        antipodal_ev = None

        for v in vertices:
            if not ((v.x == v1.x and v.y == v1.y) or (v.x == v2.x and v.y == v2.y)):
                ev = Evtype(v, Segment(v1, v2))
                distance = ev.width()
                if distance > max_distance:
                    max_distance = distance
                    antipodal_ev = ev

        antipodal_evs.append(antipodal_ev)

    return antipodal_evs


def simple_mep(evs):
    '''
    @evs: list[Evtype]
        list of antipodal pairs of the polygon
    '''
    min_area = 1e10
    mep = None
    ev1 = None
    ev2 = None

    for i in range(len(evs) - 1):
        for j in range(i, len(evs)):
            pargram = Parallelogram(evs[i], evs[j])

            if not pargram.drawable():
                continue

            area = pargram.area()
            if area < min_area:
                min_area = area
                mep = pargram
                ev1 = evs[i]
                ev2 = evs[j]

    return mep, ev1, ev2


if __name__ == '__main__':
    pts = [[6, 6], [4, 5], [3, 4], [4, 2], [7, 1], [8, 2], [7, 5]]
    convex_polygon = [Point(pt[0], pt[1]) for pt in pts]

    antipodal_evs = antipodal_pairs(convex_polygon)
    mep, ev1, ev2 = simple_mep(antipodal_evs)
    print(ev1)
    print(ev2)
    a, b, c, d = mep.a, mep.b, mep.c, mep.d()
    print(a, b, c, d)
