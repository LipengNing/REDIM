function REDIM_Theta2DWI(Theta,fn_mask,prefix,bval)
% function REDIM_Theta2DWI(Theta,Mask,pre_fix,order)
% convert REDIM parameters to relaxation-regressed nifty image
% Theta: estimation parameters
% fn_mask: file name of brain mask
% prefix: prefix of output file name
% bval: b-values of estimated data

mask_nii = MRIread(fn_mask);
mask = mask_nii.vol;
[nx,ny,nz,N] = size(Theta);

Nu = (N-4)/6; % the number of gradient directions

mask = mask(:);

Theta = reshape(Theta,nx*ny*nz,N);

nb = length(bval);

dwi = zeros(nx*ny*nz,nb,Nu);
S0 = zeros(nx*ny*nz,Nu);

parfor n = 1:nx*ny*nz
    if(mask(n))
        theta = Theta(n,:);
        if((~isnan(theta)) )
            theta_nb = UnPackParam(theta);
            for k = 1:Nu
                redim = vm_theta2param_basic(theta_nb(k,:));
                s = redim.s0*exp(-redim.d*bval+0.5*redim.d2*bval.^2);
                dwi(n,:,k) = s;
                S0(n,k) = redim.s0;
            end
        end
    end
end

S0 = mean(S0,2);
dwi = permute(dwi,[1 3 2]);
dwi = reshape(dwi,nx*ny*nz,Nu*nb);
dwi = cat(2,S0,dwi);
dwi = reshape(dwi,nx,ny,nz,[]);
mask_nii.vol = dwi;
fn_output = [prefix '_RelaxRegressed_dwi.nii.gz'];
MRIwrite(mask_nii,fn_output);


end

function theta_nb = UnPackParam(theta)
N = length(theta);
Nu = (N-4)/6;
theta_nb = zeros(Nu,10);

for i = 1:Nu
        theta_b = theta(5+[6*(i-1):(6*i-1)]); 
        % log s0, r, d, r^2, d^2, rd, r3, d3, r2d, rd2
        theta_nb(i,:) = [theta(1:2) theta_b(1) theta(3) theta_b(2:3) theta(4) theta_b(4:6)];
end
end