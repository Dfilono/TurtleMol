import numpy as np
from scipy.spatial import cKDTree

class PendingKDManager:
    def __init__(self, radii, rebuildRate=500):
        self.radii = radii
        self.rebuildRate = int(rebuildRate)

        self.kdTree = None
        self._treeXYZ = np.empty((0,3), dtype=np.float64)
        self._treeR = np.empty((0,), dtype=np.float64)

        self._pendingXYZ = np.empty((0,3), dtype=np.float64)
        self._pendingR = np.empty((0,), dtype=np.float64)
        
        # Track max radius in comitted tree for KD queries
        self._treeRMax = 0.0

    def _atomsXYZR(self, mol):
        xyz = np.array([[a[1], a[2], a[3]] for a in mol], dtype=np.float64)
        r = np.array([self.radii[a[0].title()] for a in mol], dtype=np.float64)

        return xyz, r
    
    def overlapsMolecule(self, mol, tol):
        # If empty, can't overlap
        if self.kdTree is None and self._pendingXYZ.shape[0] == 0:
            return False
        
        tol = float(tol)

        for atom in mol:
            p = np.array([atom[1], atom[2], atom[3]], dtype=np.float64)
            pr = float(self.radii[atom[0]])

            # Check against comitted KDTree
            if self.kdTree is not None:
                query_r = pr + self._treeRMax + tol
                idxs = self.kdTree.query_ball_point(p, r=query_r)
                if idxs:
                    neighXYZ = self._treeXYZ[idxs]
                    neighR = self._treeR[idxs]
                    d = neighXYZ - p[None, :]
                    d2 = np.einsum("ij,ij->i", d, d)
                    cutoff = neighR + pr + tol
                    if np.any(d2 < cutoff * cutoff):
                        return True
            
            # Check against pending atoms
            if self._pendingXYZ.shape[0] != 0:
                d = self._pendingXYZ - p[None, :]
                d2 = np.einsum("ij,ij->i", d, d)
                cutoff = self._pendingR + pr + tol
                if np.any(d2 < cutoff * cutoff):
                    return True

        return False
    
    def addMolecule(self, mol):
        xyz, r = self._atomsXYZR(mol)
        if xyz.size == 0:
            return
        self._pendingXYZ = np.vstack([self._pendingXYZ, xyz])
        self._pendingR = np.concatenate([self._pendingR, r])

    def maybeRebuild(self, force=False):
        if (not force) and (self._pendingXYZ.shape[0] < self.rebuildRate):
            return
        
        # Merge pending into comitted arrays
        if self._pendingXYZ.shape[0] != 0:
            self._treeXYZ = np.vstack([self._treeXYZ, self._pendingXYZ])
            self._treeR = np.concatenate([self._treeR, self._pendingR])
            self._pendingXYZ = np.empty((0,3), dtype=np.float64)
            self._pendingR = np.empty((0,), dtype=np.float64)

        if self._treeXYZ.shape[0] == 0:
            self.kdTree = None
            self._treeRMax = 0.0
        else:
            self.kdTree = cKDTree(self._treeXYZ)
            self._treeRMax = float(np.max(self._treeR))