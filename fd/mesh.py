class Mesh(object):
    def __init__(self, n, dx):
        self.n = n
        self.dx = dx
        self.cell_volume = dx[0] * dx[1] * dx[2]

    def __str__(self):
        return "%dx%dx%d_%gx%gx%g" % (self.n + self.dx)
