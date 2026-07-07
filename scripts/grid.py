import numpy as np
from numba import jit

class Grid(object):

    def __init__(self, cell, mesh, grid):

        self._cell = cell
        self._mesh = mesh
        self._grid = grid

        nsm = mesh[0] * mesh[1] * mesh[2]
        vol = abs(np.dot(np.cross(cell[0],cell[1]), cell[2]))

        self._vol = vol
        self._dvol = vol / nsm

        dvector = np.copy(cell)
        dvector[0] = cell[0] / mesh[0]
        dvector[1] = cell[1] / mesh[1]
        dvector[2] = cell[2] / mesh[2]
        self._dvector = dvector


    def __add__(self, other):

        if isinstance(other, Grid):
            if ((self._cell == other._cell).all() and
                (self._mesh == other._mesh).all()):

                cell = self._cell
                mesh = self._mesh

                grid1 = self._grid
                grid2 = other._grid
                grid3 = grid1 + grid2

                return Grid(cell, mesh, grid3)

    def __sub__(self, other):

        if isinstance(other, Grid):
            if ((self._cell == other._cell).all() and
                (self._mesh == other._mesh).all()):

                cell = self._cell
                mesh = self._mesh
        
                grid1 = self._grid
                grid2 = other._grid
                grid3 = grid1 - grid2

                return Grid(cell, mesh, grid3)

    def __mul__(self, other):

        if isinstance(other, Grid):
            if ((self._cell == other._cell).all() and
                (self._mesh == other._mesh).all()):

                cell = self._cell
                mesh = self._mesh
                dvol = self._dvol

                grid1 = self._grid
                grid2 = other._grid
                integ = np.sum(grid1 * grid2 * dvol)

                return integ

    def sum(self):

        dvol = self._dvol
        grid = self._grid
        return np.sum(grid * dvol)


    @jit(nopython=True)
    def planeaverage(self, axis: int = 2):

        cell = self._cell
        grid = self._grid
        mesh = self._mesh
        size = np.shape(grid)

        x = np.zeros(mesh[axis])
        y = np.zeros(mesh[axis])
        l = np.linalg.norm(cell[axis])
        count = np.prod(size) / size[axis] # number of grid points along the plane

        # sum values along the given axis
        for iz in range(size[axis]):
            plane_sum = 0
            for isp in range(size[0]):
                for ix in range(size[(axis+1) % 3]):
                    for iy in range(size[(axis+2) % 3]):
                        if axis == 2:
                            plane_sum += grid[isp, ix, iy, iz]
                        elif axis == 1:
                            plane_sum += grid[isp, ix, iz, iy]
                        elif axis == 0:
                            plane_sum += grid[isp, iz, ix, iy]

            x[iz] = iz * (l / mesh[axis])  # coordinate along the axis
            y[iz] = plane_sum / count      # averaged value

        return x, y


    def line_cut(self, axis:int =2):

        cell = self._cell
        grid = self._grid
        mesh = self._mesh

        y = np.zeros(mesh[axis])
        l = np.linalg.norm(cell[axis])
        grid = np.sum(grid, axis=0)

        if axis == 2:
            y = grid[int(mesh[0]/2)-1,int(mesh[1]/2)-1,:]
        elif axis == 1:
            y = grid[int(mesh[0]/2)-1,:,int(mesh[2]/2)-1]
        elif axis == 0:
            y = grid[:,int(mesh[1]/2)-1,int(mesh[2]/2)-1]


        x = np.arange(mesh[axis]) * (l / mesh[axis])  # coordinate along the axis

        return x, y



    def plane_cut(self, axis:int =2):

        cell = self._cell
        grid = self._grid
        mesh = self._mesh

        grid = np.sum(grid, axis=0)

        if axis == 0:
            y = grid[int(mesh[0]/2)-1,:,:]
        elif axis == 1:
            y = grid[:,int(mesh[1]/2)-1,:]
        elif axis == 2:
            y = grid[:,:,int(mesh[2]/2)-1]

        return y




