import os
import random
import re
import sys

import numpy as np
from sklearn.cluster import KMeans


def celllist(r, h, Ns, Ny=None, Nz=None):
    """Construct cell list in 3D

    This function takes the **real coordinates** of atoms `r` and the
    simulation box size `h`. Grouping atoms into Nx x Ny x Nz cells

    Parameters
    ----------
    r : float, dimension (nparticles, nd)
        *real* coordinate of atoms
    h : float, dimension (nd, nd)
        Periodic box size h = (c1|c2|c3)
    Ns : tuple, dimension (nd, )
        number of cells in x, y, z direction
    Ny : int
        if not None, represent number of cells in y direction, use with Nx = Ns
    Nz : int
        if not None, represent number of cells in z direction

    Returns
    -------
    cell : list, dimension (Nx, Ny, Nz)
        each element cell[i][j][k] is also a list recording all
        the indices of atoms within the cell[i][j][k].
        (0 <= i < Nx, 0 <= j < Ny, 0 <= k < Nz)
    cellid : int, dimension (nparticles, nd)
        for atom i:
        ix, iy, iz = (cellid[i, 0], cellid[i, 1], cellid[i, 2])
        atom i belongs to cell[ix][iy][iz]

    """
    if Ny is not None:
        if Nz is None:
            Ns = (Ns, Ny)
        else:
            Ns = (Ns, Ny, Nz)

    nparticle, nd = r.shape
    if nd != 3 or len(Ns) != 3:
        raise TypeError("celllist: only support 3d cell")

    # create empty cell list of size Nx x Ny x Nz
    cell = np.empty(Ns, dtype=object)
    for i, v in np.ndenumerate(cell):
        cell[i] = []

    # find reduced coordinates of all atoms
    s = np.dot(np.linalg.inv(h), r.T).T
    # fold reduced coordinates into [0, 1) as scaled coordinates
    s = s - np.floor(s)

    # create cell list and cell id list
    cellid = np.floor(s * np.array(Ns)[np.newaxis, :]).astype(np.int32)
    for i in range(nparticle):
        cell[tuple(cellid[i, :])].append(i)

    return cell, cellid


def verletlist(r, h, rv, atoms=None, near_neigh=None, vectorization=True):
    """Construct Verlet List (neighbor list) in 3D (vectorized)

    Uses celllist to achieve O(N)

    Parameters
    ----------
    r : float, dimension (nparticles, 3)
        *real* coordinate of atoms
    h : float, dimension (3, 3)
        Periodic box size h = (c1|c2|c3)
    rv : float
        Verlet cut-off radius

    Returns
    -------
    nn : int, dimension (nparticle, )
        nn[i] is the number of neighbors for atom i
    nindex : list, dimension (nparticle, nn)
        nindex[i][j] is the index of j-th neighbor of atom i,
        0 <= j < nn[i].

    """
    nparticles, nd = r.shape
    if nd != 3:
        raise TypeError("celllist: only support 3d cell")
    if atoms is not None:
        nparticles = atoms
    # first determine the size of the cell list
    c1 = h[:, 0]
    c2 = h[:, 1]
    c3 = h[:, 2]
    V = np.abs(np.linalg.det(h))
    hx = np.abs(V / np.linalg.norm(np.cross(c2, c3)))
    hy = np.abs(V / np.linalg.norm(np.cross(c3, c1)))
    hz = np.abs(V / np.linalg.norm(np.cross(c1, c2)))

    # Determine the number of cells in each direction
    Nx = np.floor(hx / rv).astype(np.int32)
    Ny = np.floor(hy / rv).astype(np.int32)
    Nz = np.floor(hz / rv).astype(np.int32)
    if Nx > 50:
        Nx = 12  # (Nx/20).astype(np.int)
    if Ny > 50:
        Ny = 12  # (Ny/20).astype(np.int)
    if Nz > 50:
        Nz = 12  # (Nz/20).astype(np.int)

    if Nx < 2 or Ny < 2 or Nz < 2:
        raise ValueError("Number of cells too small! Increase simulation box size.")

    # Inverse of the h matrix
    hinv = np.linalg.inv(h)
    cell, cellid = celllist(r, h, Nx, Ny, Nz)

    # initialize Verlet list
    nn = np.zeros(nparticles, dtype=int)
    nindex = [[] for i in range(nparticles)]
    if near_neigh is not None:
        global_nbr = []
    for i in range(nparticles):
        # position of atom i
        ri = r[i, :].reshape(1, 3)
        if near_neigh is not None:
            nbr_inds = []
            nbr_dist = []
        # find which cell (ix, iy, iz) that atom i belongs to
        ix, iy, iz = (cellid[i, 0], cellid[i, 1], cellid[i, 2])

        # go through all neighboring cells
        ixr = ix + 1
        iyr = iy + 1
        izr = iz + 1
        if Nx < 3:
            ixr = ix
        if Ny < 3:
            iyr = iy
        if Nz < 3:
            izr = iz

        for nx in range(ix - 1, ixr + 1):
            for ny in range(iy - 1, iyr + 1):
                for nz in range(iz - 1, izr + 1):
                    # apply periodic boundary condition on cell id nnx, nny, nnz
                    nnx, nny, nnz = (nx % Nx, ny % Ny, nz % Nz)

                    # extract atom id in this cell
                    ind = cell[nnx][nny][nnz].copy()
                    nc = len(ind)

                    # vectorized implementation
                    if vectorization:
                        if i in ind:
                            ind.remove(i)
                        rj = r[ind, :]
                        drij = pbc(rj - np.repeat(ri, len(ind), axis=0), h, hinv)

                        ind_nbrs = np.where(np.linalg.norm(drij, axis=1) < rv)[
                            0
                        ].tolist()
                        if near_neigh is None:
                            if len(ind_nbrs) > 0:
                                nn[i] += len(ind_nbrs)
                                nindex[i].extend([ind[j] for j in ind_nbrs])
                        else:
                            if len(ind_nbrs) > 0:
                                nbr_inds.extend([ind[j] for j in ind_nbrs])
                                nbr_dist.extend(
                                    np.linalg.norm(drij[ind_nbrs], axis=1).tolist()
                                )
                    else:
                        for k in range(nc):
                            j = ind[k]
                            # update nn[i] and nindex[i]
                            if i == j:
                                continue
                            else:
                                rj = r[j, :].reshape(1, 3)

                                # obtain the distance between atom i and atom j
                                drij = pbc(rj - ri, h, hinv)
                                if near_neigh is None:
                                    if np.linalg.norm(drij) < rv:
                                        nn[i] += 1
                                        nindex[i].append(j)
                                else:
                                    if np.linalg.norm(drij) < rv:
                                        nbr_dist.extend(
                                            np.linalg.norm(drij).reshape(
                                                -1,
                                            )
                                        )
                                        nbr_inds.append(j)
        if near_neigh is not None:
            if len(nbr_dist) >= near_neigh:
                nbr_dist = np.array(nbr_dist)
                nbr_inds = np.array(nbr_inds)
                nbr_inds = nbr_inds[np.argsort(nbr_dist)].tolist()
                to_add = []
                ct = 0
                for ll in range(len(nbr_inds)):
                    if (nbr_inds[ll] not in global_nbr) and ct < near_neigh:
                        to_add.append(nbr_inds[ll])
                        global_nbr.append(nbr_inds[ll])
                        ct += 1
                if ct == near_neigh:
                    nn[i] = near_neigh
                    nindex[i].append(to_add)

    return nn, nindex


def find_quartic_bonds(pos, h, rv, n_bead):
    """
    Find quartic bonds in a given position array.

    Parameters:
        pos (ndarray): The position array.
        h (float): The smoothing parameter.
        rv (float): The cutoff radius.
        n_bead (int): The number of beads in a chain.

    Returns:
        ndarray: An array of quartic bonds, where each bond is represented as a pair of atom indices.
    """
    n_atoms = pos.shape[0]
    nn, nindex = verletlist(pos, h, rv, near_neigh=3, vectorization=True)
    bonds = []
    tracks = []
    for i in range(n_atoms):
        if nn[i] > 0:
            chain_no = i // n_bead
            rej = np.linspace(
                chain_no * n_bead, (chain_no + 1) * n_bead - 1, n_bead, dtype=int
            )
            v_id = [l for l in nindex[i][0] if l not in rej]
            if len(v_id) > 0:
                if i + 1 not in tracks and v_id[0] + 1 not in tracks:
                    bonds.append([i + 1, v_id[0] + 1])
                    tracks.append(i + 1)
                    tracks.append(v_id[0] + 1)

    return np.array(bonds)


def uniform_bond_indices(bond_pos, n_bonds):
    """
    Given an array of bond positions and the number of desired bonds, this function uses K-means clustering to group the bond positions into clusters. It then selects the bond position closest to each cluster center as the representative bond for that cluster. The function returns an array of indices corresponding to the representative bonds.

    Parameters:
    - bond_pos (numpy.ndarray): An array of shape (n_bonds, n_dimensions) representing the positions of the bonds.
    - n_bonds (int): The desired number of representative bonds to select.

    Returns:
    - indices (numpy.ndarray): An array of shape (n_bonds,) containing the indices of the representative bonds.
    """
    kmeans = KMeans(n_clusters=n_bonds, init="k-means++").fit(bond_pos)
    centroids = kmeans.cluster_centers_
    indices = []
    for i in range(len(centroids)):
        distances = np.linalg.norm(bond_pos - centroids[i], axis=1)
        indices.append(np.argmin(distances))

    return indices


def excluded_volume_indices(scaled_pos, h, n_bonds):
    """
    Generate a list of indices for atoms that are excluded from bonded interactions.

    Parameters:
        scaled_pos (ndarray): The scaled positions of the atoms.
        h (ndarray): The lattice vectors.
        n_bonds (int): The number of bonds in the system.

    Returns:
        list: A list of indices for atoms that are excluded from bonded interactions.
    """
    f_pos = scaled_pos @ h
    np.random.shuffle(f_pos)
    r_cut = (np.linalg.det(h) * 3 / 4 / np.pi / n_bonds) ** (1 / 3) / 1.1
    nn_b, nindex_b = verletlist(f_pos, h, r_cut)
    selected_ind = []
    rejected_ind = []
    for i in range(f_pos.shape[0]):
        if i not in rejected_ind:
            selected_ind.append(i)
            rejected_ind.extend(nindex_b[i])

    return selected_ind


def acceptor_links(chain_1, chain_2, n_bead):
    """
    Generate the acceptor links for a given pair of chains.

    Parameters:
        chain_1 (int): The index of the first chain.
        chain_2 (int): The index of the second chain.
        n_bead (int): The number of beads in each chain.

    Returns:
        list: A list of acceptor links between the two chains.
    """
    acc = []
    acc1 = np.linspace(chain_1 * n_bead, (chain_1 + 1) * n_bead - 1, n_bead, dtype=int)
    acc2 = np.linspace(chain_2 * n_bead, (chain_2 + 1) * n_bead - 1, n_bead, dtype=int)
    acc.extend(acc1)
    acc.extend(acc2)
    return acc


def rejector_links(site_1, site_2, spacing):
    """
    Generates an array of sites around site_1 and site_2 with a given spacing.

    Parameters:
        site_1 (int): The first site value.
        site_2 (int): The second site value.
        spacing (int): The spacing value for generating sites around site_1 and site_2.

    Returns:
        list: An array of sites around site_1 and site_2 based on the spacing.
    """
    site = []
    s1 = np.linspace(site_1 - spacing, site_1 + spacing, 2 * spacing + 1, dtype=int)
    s2 = np.linspace(site_2 - spacing, site_2 + spacing, 2 * spacing + 1, dtype=int)
    site.extend(s1)
    site.extend(s2)
    return site


def uniform_spacing_indices(cand_bonds, n_bonds, n_bead):
    """
    Generate a list of indices for candidate bonds with uniform spacing.

    Parameters:
        cand_bonds (ndarray): The candidate bonds.
        n_bonds (int): The number of bonds.
        n_bead (int): The number of beads.

    Returns:
        list: A list of indices for candidate bonds with uniform spacing.
    """
    spacing = 500 * 500 // (4 * n_bonds)
    np.random.shuffle(cand_bonds)
    rej = []
    selected_ind = []
    for i in range(cand_bonds.shape[0]):
        if cand_bonds[i, 0] not in rej and cand_bonds[i, 1] not in rej:
            selected_ind.append(i)
            chain_1 = cand_bonds[i, 0] // n_bead
            chain_2 = cand_bonds[i, 1] // n_bead
            acc = acceptor_links(chain_1, chain_2, n_bead)
            rejs = rejector_links(cand_bonds[i, 0], cand_bonds[i, 1], spacing)
            rej.extend(np.intersect1d(acc, rejs))

    return selected_ind


def pbc(drij, h, hinv=None):
    """calculate distance vector between i and j

    Considering periodic boundary conditions (PBC)

    Parameters
    ----------
    drij : float, dimension (npairs, 3)
        distance vectors of atom pairs (Angstrom)
    h : float, dimension (3, 3)
        Periodic box size h = (c1|c2|c3)
    hinv : optional, float, dimension (3, 3)
        inverse matrix of h, if None, it will be calculated

    Returns
    -------
    drij : float, dimension (npairs, 3)
        modified distance vectors of atom pairs considering PBC (Angstrom)

    """
    # Check the input
    if len(drij.shape) == 1:  # Only one pair
        drij = drij.reshape(1, -1)
    if len(drij.shape) != 2:
        raise ValueError(
            "pbc: drij shape not correct, must be (npairs, nd), (nd = 2,3)"
        )
    npairs, nd = drij.shape
    if len(h.shape) != 2 or h.shape[0] != h.shape[1] or nd != h.shape[0]:
        raise ValueError("pbc: h matrix shape not consistent with drij")
    # Calculate inverse matrix of h
    if hinv is None:
        hinv = np.linalg.inv(h)

    dsij = np.dot(hinv, drij.T).T
    dsij = dsij - np.round(dsij)
    drij = np.dot(h, dsij.T).T

    return drij
