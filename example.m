Input_data = {'Data_TE1.nii.gz', 'Data_TE2.nii.gz', 'Data_TE3.nii.gz', 'Data_TE4.nii.gz'};
Output_prefix = 'XXX/output_prefix';
fn_mask = 'brain_mask.nii.gz';
fn_mask = 'TestMask.nii.gz';

fn_bvec = 'dwi.bvec';
fn_bval = 'dwi.bval';
TE = [71 101 131 161];% change it to the correct values


redim_pipe(Input_data, Output_prefix, fn_mask, fn_bvec, fn_bval, TE);
