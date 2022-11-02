import numpy as np


class ViewingCone(object):
    """
    A class representing viewing cones obtained by LINE21.
    """

    def __init__(self, name):

        if name.lower() == 'km11' or name.lower() == 'tofor':
            file_path = 'vc_data/TOFOR/TOFOR.vc'
        elif name.lower() == 'ab':
            file_path = 'vc_data/AB/AB.vc'
        else:
            file_path = name

        self.file = file_path

        VC = np.loadtxt(file_path)

        self.X   = VC[:,0]
        self.Y   = VC[:,1]
        self.Z   = VC[:,2]
        self.C   = VC[:,3]
        self.V   = VC[:,4]
        self.UX  = VC[:,5]
        self.UY  = VC[:,6]
        self.UZ  = VC[:,7]
        self.R   = VC[:,8]
        self.PHI = VC[:,9]
        if VC.shape[1]==11:
            self.OMEGA = VC[:,10]
        else:
            self.OMEGA = 4*np.pi * self.C / self.V

        self.Nvoxels = len(self.C)
        self.Vtot = np.sum(self.V)

        # Transform detector direction vectors to cylindrical coordinates
        self.UR   = self.UX*np.cos(self.PHI) + self.UY*np.sin(self.PHI)
        self.Uphi = -self.UX*np.sin(self.PHI) + self.UY*np.cos(self.PHI)

        # Collect cylindrical direction vectors
        self.U = np.vstack((self.UR, self.Uphi, self.UZ))

        # Reconstruct the grid used by LINE21
        self.Rvals = np.unique(self.R)
        self.Zvals = np.unique(self.Z)
        self.PHIvals = np.unique(self.PHI)
        self.NR = len(self.Rvals)
        self.NZ = len(self.Zvals)
        self.NP = self.NR * self.NZ    # number of poloidal grid points
        self.NPHI = len(self.PHIvals)   # number of toroidal grid points

        self.dR = self.Rvals[1] - self.Rvals[0]    # assumes uniform grid spacing
        self.dZ = self.Zvals[1] - self.Zvals[0]
        self.dPHI = self.PHIvals[1] - self.PHIvals[0]

        # Map each voxel to the corrresponding poloidal grid point
        self.IR =  np.round((self.R - self.Rvals[0]) / self.dR).astype(int)
        self.IZ =  np.round((self.Z - self.Zvals[0]) / self.dZ).astype(int)
        self.IP = self.IR + self.IZ*self.NR

        self.IPvals = np.unique(self.IP)

        # Map each poloidal grid point to the corresponding voxels
        IVOX = np.empty(self.NP, dtype='object')
        for i in range(self.NP):
            IVOX[i] = []

        for ivox, ip in enumerate(self.IP):
            IVOX[ip].append(ivox)

        self.IVOX = IVOX

        # Also count the number of voxels at each poloidal grid point
        NVOX = np.zeros(self.NP, dtype='int')

        for i in range(self.NP):
            NVOX[i] = len(IVOX[i])

        self.NVOX = NVOX

        # Generate poloidal projection of the viewing cone
        IRP = self.IPvals%self.NR
        IZP = self.IPvals//self.NR
        self.RP = self.Rvals[IRP]
        self.ZP = self.Zvals[IZP]

        CP  = np.zeros(len(self.IPvals))

        for i,ip in enumerate(self.IPvals):
            # Total weight of all voxels at the given poloidal location
            CP[i]  = np.sum(self.C[IVOX[ip]])

        self.CP = CP


    def __repr__(self):
        return 'Viewing cone: {}'.format(self.file)


    # Class methods
    def get_poloidal_index(self, r, z):
        """ Check which poloidal bin that we are in.
        Requires that the poloidal grid spacing is uniform."""
        r = np.atleast_1d(r)
        z = np.atleast_1d(z)

        ir = np.floor( ( r - self.Rvals[0] + self.dR/2.0) / self.dR ).astype(int)
        iz = np.floor( ( z - self.Zvals[0] + self.dZ/2.0) / self.dZ ).astype(int)

        ip = ir + iz*self.NR

        ip[ir<0]         = -1
        ip[ir>=self.NR] = -1
        ip[iz<0]         = -1
        ip[iz>=self.NZ] = -1

        outside = np.logical_not(np.in1d(ip, self.IP))
        ip[outside] = -1

        return ip


    def is_inside(self, r, z):
        """ Check whether the given points are inside the viewing cone. """

        ip = self.get_poloidal_index(r, z)

        is_inside = np.in1d(ip, self.IP)

        return is_inside


    def get_voxels(self, r, z):
        """ Return all voxels that include a given point in the poloidal plane. """

        ip = self.get_poloidal_index(r, z)

        return np.where(ip==-1, -1, self.IVOX[ip])


    def map_points_voxels(self, r, z):
        """
        Map given (r,z) points to each voxel that encloses them, under the
        assumption of toroidal symmetry. 
        
        Returns
        -------
        i_voxels: array of voxel indices. Each element correponds to 
            one point in the voxel with that index.

        i_points: array of point indices (same length as i_voxels). 
            For instance, i_point[j] tells us that that point is
            inside the voxel with index i_voxels[j].
            (A point with given r,z coordinates might be inside 
            multiple voxels, since toroidal symmetry is assumed).

        (DISABLED) n_points: array that, for each point-voxel pair, holds the 
            number of points in that voxel.
        """

        ip = self.get_poloidal_index(r, z)
        i_all_points = np.arange(len(ip))
        
        inside = np.in1d(ip, self.IP)
        ip = ip[inside]

        if len(ip) == 0:
            # No points inside any voxel
            return np.array([]), np.array([])

        i_points_inside = i_all_points[inside]

        ivox = self.IVOX[ip]
        nvox = self.NVOX[ip]
        
        #ip_inverse, n_counts = np.unique(ip, return_inverse=True, return_counts=True)[1:]
        #n_counts = n_counts[ip_inverse]

        i_voxels = np.concatenate(ivox)
        i_points = np.zeros_like(i_voxels)
        #n_points = np.zeros_like(i_voxels)

        i0 = 0
        for i,n in enumerate(nvox):
            i1 = i0 + n
            i_points[i0:i1] = i_points_inside[i] 
            #n_points[i0:i1] = n_counts[i]
            i0 = i1

        return i_voxels, i_points #, n_points
