function theta = redim_data2fit_1vox(signal, x, order)
% function theta_1x10 = dtd_skewness_all_data2fit(signal, x)
%
% Diffusion-Relaxation coupling estimation. 
% Input: signal is an 1xN (or Nx1) vector that is measured using different
% TE and b-values along different direction. x is a 3xN (or Nx3) matrix
% where the first row is the TE, the second row is the b-value and the
% third row is the gradient direction index
% order: the order of cumulant expansions
% Output: theta vector which consists of [log(S0) theta_mean_1x2
% theta_variance_1x3 and theta_3rd_moments_1x4]

if (nargin == 2)
    order = 3;
end

signal = signal(:);
% x(:,1):TE, x(:,2): b-value, x(:,3): directions
if(size(x,1)==3)
    x = x';
end
x = x(signal>0,:); % 
x(:,2) = x(:,2)*1e-3; % change the unit diff to ms/mu m^2

signal = signal(signal>0);

[X,A] = ParamToMatrix(x, order);% 

b = zeros(size(A,1),1);
B = A(3:2:end,:);
A = [-A; B];
b = [b; 3*ones(size(B,1),1)];% constrain the upper bound of diffusivity

options = optimoptions('lsqlin','display', 'off'); % run interior-point algorithm


theta = lsqlin(X,log(signal),A,b,[],[],[],[],[],options);
theta = theta';


end

function [X,A] = ParamToMatrix(x, order)
    if(nargin==1)
        order = 3;
    end
    
    b = unique(x(:,3)); % unique set of b-vectors
    nb = length(b);
    N = size(x,1);%the size of measurements
    if(order == 3)
        X = zeros(N,4+nb*6); 
        row_ind = 1;
        I = eye(4+nb*6);
        A = zeros(2+nb*2,4+nb*6);
        A(1:2,:) = I(2:3,:);% r, vr>0

        for i = 1:nb
            ind = find(x(:,3)==b(i));
            L = length(ind);
            %xi: t,b,t^2,b^2,tb,t^3,b^3,t^2b,tb^2
            Xi = cell2mat(cellfun(@vm_1x2_to_1x10,mat2cell(x(ind,1:2),ones(L,1),2), 'UniformOutput',false));
            
            X(row_ind:row_ind+L-1,:) =[Xi(:,[1 2 4 7]) zeros(L,6*(i-1)) Xi(:,[3 5 6 8 9 10]) zeros(L,6*(nb-i))];
            row_ind = row_ind+L;
            A((2*i+1):(2*i+2),:) = [zeros(2,4+6*(i-1)) eye(2,6) zeros(2,6*(nb-i))]; 
        end
    else
        if(order == 2)
            X = zeros(N,3+nb*3); 
            row_ind = 1;
            I = eye(3+nb*3);
            A = zeros(2+nb*2,3+nb*3);% 
            A(1:2,:) = I(2:3,:);% r, vr>0

            for i = 1:nb
                ind = find(x(:,3)==b(i));% parameters for the same direction
                L = length(ind);
                Xi = cell2mat(cellfun(@vm_1x2_to_1x6,mat2cell(x(ind,1:2),ones(L,1),2), 'UniformOutput',false));
                X(row_ind:row_ind+L-1,:) =[Xi(:,[1 2 4]) zeros(L,3*(i-1)) Xi(:,[3 5 6]) zeros(L,3*(nb-i))];
                row_ind = row_ind+L;
                A((2*i+1):(2*i+2),:) = [zeros(2,3+3*(i-1)) eye(2,3) zeros(2,3*(nb-i))]; %d, var(d) >0
            end
        end
    end
    

end