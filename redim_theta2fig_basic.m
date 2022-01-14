function redim_theta2fig_basic(Theta,Mask,pre_fix)
% function redim_theta2fig_basic(Theta,Mask,pre_fix)
% convert parameters to s0, r, r2, d,d2,rd


mask_nii = MRIread(Mask);
mask = mask_nii.vol;
[nx,ny,nz,N] = size(Theta);

Nu = (N-4)/6; % the number of gradient directions

mask = mask(:);

Theta = reshape(Theta,nx*ny*nz,N);

%measure = {'s0','r','d','r2','d2','rd'};
measure = {'r','rd'};

L = length(measure);
Measure = zeros(nx*ny*nz,L);

parfor n = 1:nx*ny*nz
    if(mask(n))
        theta = Theta(n,:);
        if((~isnan(theta)) )
            theta_nb = UnPackParam(theta);
            temp = zeros(Nu,L);
            for k = 1:Nu
                redim = vm_theta2param_basic(theta_nb(k,:));
                temp(k,1) = redim.r;
                temp(k,2) = redim.rd;
%                 
%                 for l = 1:L
%                     temp(k,l) = eval(['redim.' measure{l}]);
%                 end
            end
            
            for l = 1:L
                Measure(n,l) = mean(temp(:,l));% average over gradient directions
            end
            
        end
    end
end

for l = 1:L % export r,rd,
    m = reshape(Measure(:,l),nx,ny,nz);
    mask_nii.vol = m;
    fn_output = [pre_fix '_' measure{l} '.nii.gz'];
    MRIwrite(mask_nii,fn_output);
end

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
