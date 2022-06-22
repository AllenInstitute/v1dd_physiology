import os
import numpy as np
import NeuroAnalysisTools.core.FileTools as ft
import NeuroAnalysisTools.core.ImageAnalysis as ia
import NeuroAnalysisTools.NwbTools as nt
import utils as utils

data_folder = r'\\allen\programs\mindscope\workgroups\surround\v1dd_in_vivo_new_segmentation\data\nwbs'

fns = ft.look_for_file_list(
    source=data_folder, 
    identifiers=[''], 
    file_type='nwb',
    is_full_path=False
    )

# fns = [fn for fn in fns if fn[9] not in '12345']
fns.sort()
print('\n'.join(fns))

for fi, fn in enumerate(fns):

    print(f'processing {fn}, {fi + 1}/{len(fns)} ...')

    vasmap_wf, vasmap_2p = utils.get_vasmap(fn[0:10])
    vasmap_wf = ia.array_nor(vasmap_wf)
    vasmap_2p = ia.array_nor(vasmap_2p)

    # rotate vasmaps recorded from deepscope to match standard view
    # these transformation are only for 3p
    vasmap_wf_r = ia.rigid_transform_cv2(vasmap_wf[:, ::-1], rotation=28).astype(np.float32)
    vasmap_wf_r = ia.array_nor(vasmap_wf_r)
    vasmap_2p_r = ia.rigid_transform_cv2(vasmap_2p, rotation=30).astype(np.float32)
    vasmap_2p_r = ia.array_nor(vasmap_2p_r)

    nwb_f = nt.RecordedFile(os.path.join(data_folder, fn))

    nwb_f.add_acquisition_image('vasmap_wf', vasmap_wf,
                                description='wide field surface vasculature map '
                                            'through cranial window original')
    nwb_f.add_acquisition_image('vasmap_wf_rotated', vasmap_wf_r,
                                description='wide field surface vasculature map '
                                            'through cranial window rotated to match standard view')
    nwb_f.add_acquisition_image('vasmap_2p', vasmap_2p,
                                description='2p surface vasculature map through '
                                            'cranial window original')
    nwb_f.add_acquisition_image('vasmap_2p_rotated', vasmap_2p_r,
                                description='2p surface vasculature map through '
                                            'cranial window rotated to match standard view')
    nwb_f.close()


