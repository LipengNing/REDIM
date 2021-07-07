function redim_pipe(Input_data, Output_prefix, fn_mask, fn_bvec, fn_bval, TE,method,order)
% function redim_pipe(Input_data, Output_prefix, fn_mask, fn_bvec, fn_bval, TE,method,order)
% The main function to compute the joint relaxation-diffusion cumulant model fitting parameters
% Input_data: string array of the input files with mulitple TE
% Output_prefix: prefix of output files
% fn_mask: file name of the brain mask
% fn_bvec: file name of the b-vector file
% fn_bval: file name of the b-value file
% TE: TE values of the dwi files
% method: method of estimation algorithms:
%         'l' (linear model for logarithmic signals), 'wl' (weighted linear least
%         square), 'n' (nonlinear model for signals w/o logarithms), 'wm' (weighted
%          nonlinear: default)
% order: binary index about which moments of [r d r2 d2 rd r3 d3 r2d rd2] need to be estimated
%        default: [1 1 0 1 1 0 0 0 0]
% Lipeng Ning 07/07/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~exist('order','var'))
        order = [1 1 0 1 1 0 0 0 0];%[r d ~ d2 rd ~ ~ ~ ~]
end

if(~exist('method','var'))
        method = 'wn';%[r d r2 d2 rd r3 d3 r2d rd2]
end 



L = length(Input_data);
nii_mask = MRIread(fn_mask);
Mask = nii_mask.vol;
%%


[nx,ny,nz] = size(Mask);

bvec = load(fn_bvec);

if(size(bvec,1)==3)
    bvec = bvec';
end

bval = load(fn_bval);
bval = bval(:);
Nb = length(bval); 

if(size(bvec,1)==3)
    bvec = bvec';
end
bvec(bval==0,:) = 0;
id_bvec = zeros(Nb,1);
u = bvec(bval==min(bval(bval>0)),:);% the unique set of b vectors
Nu = size(u,1);
for k = 1:Nu
    c = abs(u(k,:)*bvec');
    id_bvec(c>0.99)=k;% label the gradient directions
end

S = zeros(nx,ny,nz,Nb,L); % all data
x = zeros(Nb,3,L);

for i = 1:L
    nii_data = MRIread(Input_data{i});
    data = nii_data.vol;
    S(:,:,:,:,i) = data;
    x(:,1,i) = TE(i)*ones(Nb,1);
    x(:,2,i) = bval;
    x(:,3,i) = id_bvec;
end

x = permute(x,[1 3 2]);
x = reshape(x,[Nb*L,3]);

Theta = REDIM_Data2Fit(S,x,Mask,method,order);
% 

save([Output_prefix '_theta.mat'], 'Theta','u'); %hthe fitting parameters

redim_theta2fig_basic(Theta,fn_mask,Output_prefix);
b_nonzero = unique(bval)*1e-3;
b_nonzero = b_nonzero(b_nonzero>0);
REDIM_Theta2DWI(Theta,fn_mask,Output_prefix,b_nonzero);
%%
bval_out = kron(b_nonzero(:),ones(Nu,1)); bval_out = [0;bval_out*1000];
dlmwrite([Output_prefix '_bvals.txt'],bval_out,'delimiter','\t');
bvec_out = repmat(u',1,length(b_nonzero)); bvec_out = [zeros(3,1) bvec_out];
dlmwrite([Output_prefix '_bvecs.txt'],bvec_out','delimiter','\t');
end