"""generate a solid file from an irregular shaped mask

"""
import logging
import os
import subprocess
from shutil import which
import numpy as np


def make_solid_file(clipped_mask, out_name, dx=1000, dz=1000):
    """Make a solid file and vtk file from a clipped mask, write to out_name, requires pfmask-to-sol from pftools

    Parameters
    ----------
    clipped_mask : ndarray
        the mask array cropped to the inner shape, with 1's in the mask, 0's outside
    out_name : str
        name for the output .vtk and .pfsol files
    dx : int, optional
        horizontal cell size (Default value = 1000)
    dz : int, optional
        vertical cell size (Default value = 1000)

    Returns
    -------
    batches : list
        list of batches (patches?) located
    """
    # TODO Add ability to buid features as patches on top layer of mask
    pf_mask_to_sol_path = find_mask_to_sol_exe()
    if pf_mask_to_sol_path is None:
        msg = 'Could not locate pfmask-to-pfsol utility needed to generate solid file (.pfsol)' \
              ' ensure PARFLOW_DIR environment variable is set'
        logging.exception(msg)
        raise Exception(msg)
    # TODO: Should not have to flip this mask. Find why it is upside down
    mask_mat = np.flip(clipped_mask, axis=1)
    #

    if len(mask_mat.shape) == 3:
        mask_mat = np.squeeze(mask_mat, axis=0)
    # create back borders
    # Back borders occur where mask[y+1]-mask[y] is negative
    # (i.e. the cell above is a zero and the cell is inside the mask, i.e. a 1)
    back_mat = np.zeros(mask_mat.shape)

    # create front borders
    # Front borders occur where mask[y-1]-mask[y] is negative
    # (i.e. the cell above is a zero and the cell is inside the mask, i.e. a 1)
    front_mat = np.zeros(mask_mat.shape)

    # create left borders
    # Left borders occur where mask[x-1]-mask[x] is negative
    left_mat = np.zeros(mask_mat.shape)

    # create right borders
    # Right borders occur where mask[x+1]-mask[x] is negative
    right_mat = np.zeros(mask_mat.shape)

    # deal with top and bottom patches
    # 3 = regular overland boundary
    # 4 = Lake
    # 5 = Sink
    # 6 = bottom
    # 8 = Stream
    # 9 = Reservoir

    top_mat = mask_mat * 3

    bottom_mat = mask_mat * 6

    # create header and write ascii files
    header = 'ncols ' + str(mask_mat.shape[1]) + '\n'
    header += 'nrows ' + str(mask_mat.shape[0]) + '\n'
    header += 'xllcorner 0.0\n'
    header += 'yllcorner 0.0\n'
    header += 'cellsize ' + str(dx) + '\n'
    header += 'NODATA_value 0.0\n'

    patches = {os.path.join(os.path.dirname(out_name), 'Back_Border.asc'): back_mat,
               os.path.join(os.path.dirname(out_name), 'Front_Border.asc'): front_mat,
               os.path.join(os.path.dirname(out_name), 'Right_Border.asc'): right_mat,
               os.path.join(os.path.dirname(out_name), 'Left_Border.asc'): left_mat,
               os.path.join(os.path.dirname(out_name), 'Bottom_Border.asc'): bottom_mat,
               os.path.join(os.path.dirname(out_name), 'Top_Border.asc'): top_mat}

    list_patches = list(patches.keys())
    for patch in patches:
        with open(patch, 'w') as fo:
            fo.write(header)
            np.savetxt(fo, patches[patch].reshape([-1, 1]), '%.3f', ' ')

    # Create solid domain file
    # This part assumes that you have installed pf-mask-utilities or ParFlow with PARFLOW_DIR set
    out_vtk = out_name + '.vtk'
    out_pfsol = out_name + '.pfsol'

    if os.path.isfile(out_vtk):
        os.remove(out_vtk)

    if os.path.isfile(out_pfsol):
        os.remove(out_pfsol)

    sub_cmd_str = [pf_mask_to_sol_path[0],
                   '--mask-back', list_patches[0],
                   '--mask-front', list_patches[1],
                   '--mask-right', list_patches[2],
                   '--mask-left ' + list_patches[3],
                   '--mask-bottom ' + list_patches[4],
                   '--mask-top ' + list_patches[5],
                   '--vtk', out_vtk,
                   '--pfsol', out_pfsol,
                   pf_mask_to_sol_path[1], str(dz)]
    logging.info(f'begin mask_to_sol subprocess, command executed: {sub_cmd_str}')
    create_sub = subprocess.run(sub_cmd_str, stdout=subprocess.PIPE)
    temp_list = create_sub.stdout.decode('utf-8').split('\n')
    batches = ''
    for line in temp_list:
        if 'Number of triangles in patch' in line:
            line = line.strip()
            batches += line.split()[-3] + ' '
    logging.info(f'identified batches in domain {batches}')
    return batches


def find_mask_to_sol_exe():
    """find the pf_mask_to_pfsol or mask-to-pfsol utility on the system

    Returns
    -------
    pf_mask_to_sol_path : tuple
        tuple of (path to executable, depth/z-bottom flag for argument) or None if no executable was found
    """
    pf_mask_to_sol_path = None
    possible_paths = {('mask-to-pfsol', '--depth'): [os.environ.get('PFMASKUTILS'), which('mask-to-pfsol')],
                      ('pfmask-to-pfsol', '--z-bottom'):
                          [f'{os.environ.get("PARFLOW_DIR")}/bin', which('pfmask-to-pfsol')]}
    for executable, possible_paths in possible_paths.items():
        executable_path = next((path for path in possible_paths if path is not None), None)
        if executable_path is not None:
            if os.path.isfile(os.path.join(executable_path, executable[0])):
                pf_mask_to_sol_path = (os.path.join(executable_path, executable[0]), executable[1])
                break
    if pf_mask_to_sol_path is None:
        logging.exception('Could not locate mask to solid utility!')
    else:
        logging.info(f'searching for mask_to_sol executable resulted in: {pf_mask_to_sol_path[0]}')
    return pf_mask_to_sol_path
