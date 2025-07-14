import numpy as np
import os
import h5py
import mdtraj as md
import tqdm
import ray
 
ray.init()
 
@ray.remote
def count_nS(traj_segs_path, x):
    counts = 0
    iter_dir = os.path.join(traj_segs_path, f'{x:06d}')
    if os.path.exists(iter_dir):
        counts = len([name for name in os.listdir(iter_dir)])
    return counts
 
@ray.remote
def process_segment(n_iter, iS, trajpath, refPDBfile, parentTraj, childTraj):
    try:
        coord0 = np.squeeze(md.load(f'{trajpath}/{parentTraj}', top=refPDBfile)._xyz)
    except OSError:
        coord0 = np.squeeze(md.load(refPDBfile, top=refPDBfile)._xyz)
 
    coord1 = np.squeeze(md.load(f'{trajpath}/{childTraj}', top=refPDBfile)._xyz)
 
    return iS, coord0, coord1
 
def main():
    traj_segs_path = '/projects/mdpoleto_data/share/4marcelo/data/westpa-drude/traj_segs'
    max_iter = 250
    refPDBfile = '/projects/mdpoleto_data/share/4marcelo/data/westpa-drude/bstates/01/basis.pdb'
    WEfolder = '/projects/mdpoleto_data/share/4marcelo/data/westpa-drude/'
    h5file = '/projects/mdpoleto_data/share/4marcelo/data/westpa-drude/west_with_coords.h5'
 
    print(f'Doing collectCoordinates on  WE file {h5file}')
    print('Preparing coordinates...')
 
    parentTraj = 'parent.pdb'
    childTraj = 'seg.pdb'
 
    top_temp = md.load(refPDBfile, top=refPDBfile)
    n_atoms = top_temp.topology.n_atoms
 
    with h5py.File(h5file, 'a') as f:
        for n_iter in tqdm.tqdm(range(1, max_iter + 1)):
            nS_future = count_nS.remote(traj_segs_path, n_iter)
            nS = ray.get(nS_future)
 
            coords = np.zeros((nS, 2, n_atoms, 3))
            dsetName = "/iterations/iter_%08d/auxdata/coord" % int(n_iter)
 
            coords_exist = False
            try:
                dset = f.create_dataset(dsetName, np.shape(coords))
            except (RuntimeError, ValueError):
                print('coords exist for iteration ' + str(n_iter) + ' NOT overwritten')
                coords_exist = True
                continue
 
            segment_futures = []
            for iS in range(nS):
                trajpath = WEfolder + "/traj_segs/%06d/%06d" % (n_iter, iS)
                segment_futures.append(process_segment.remote(n_iter, iS, trajpath, refPDBfile, parentTraj, childTraj))
 
            results = ray.get(segment_futures)
            for iS, coord0, coord1 in results:
                coords[iS, 0, :, :] = coord0
                coords[iS, 1, :, :] = coord1
 
            if not coords_exist:
                dset[:] = coords
 
            print(f"Wrote coords for iteration {n_iter}.")
 
if __name__ == "__main__":
    main()
