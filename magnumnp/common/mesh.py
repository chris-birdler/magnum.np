from magnumnp.common import logging

__all__ = ["Mesh"]

class Mesh(object):
    def __init__(self, n, dx, origin=(0,0,0)):
        self.n = n
        self.dx = dx
        self.origin = origin
        self.cell_volume = dx[0] * dx[1] * dx[2]
        logging.info_green("[Mesh] initialize mesh %dx%dx%d (size= %g x %g x %g)" % (self.n + self.dx))

    def __str__(self):
        return "%dx%dx%d_%gx%gx%g" % (self.n + self.dx)
