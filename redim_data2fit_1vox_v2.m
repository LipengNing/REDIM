function theta = redim_data2fit_1vox_v2(signal, x, order)
% function theta = redim_data2fit_1vox_v2(signal, x,order)
%
% Diffusion-Relaxation coupling estimation. 
% Input: signal is an 1xN (or Nx1) vector that is measured using different
% TE and b-values along different direction. x is a 3xN (or Nx3) matrix
% where the first row is the TE, the second row is the b-value and the
% third row is the gradient direction index
% order=[order_r, order_d]: the maximum orders for r and d in the cumulant expansions
% Output: theta vector which consists of [log(S0) theta_mean_1x2
% theta_variance_1x3 and theta_3rd_moments_1x4]
% update for _3echo: only the first order model is used to model relaxation
% LOR current only works when order = 2. Update is needed for order = 3

if (nargin == 2)
    order = [1 2]; % 1st order relaxation and 2nd order diffusion
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
%add constrains
b = zeros(size(A,1),1);%r, vr,d,vd >0
b = unique(x(:,3)); % unique set of b-vectors
nb = length(b);
nd = (size(A,2)-order(1)-1)/nb;% number of gradient-specific parameters
B = A(order(1)+2:nd:end,:); % diffusivity smaller than 3
A = [-A; B];
b = [b; 3*ones(size(B,1),1)];% constrain the upper bound of diffusivity

options = optimoptions('lsqlin','display', 'off'); % run interior-point algorithm

W = diag(signal.^2);
theta = lsqlin(W*X,W*log(signal),A,b,[],[],[],[],[],options);
theta = theta';
% reformat theta to a 1x10 vector
theta=ConvertFullParam(theta,x,order);
end

function [X,A] = ParamToMatrix(x, order)
    if(nargin==1)
        order = [1 2];
    end
    
    b = unique(x(:,3)); % unique set of b-vectors
    nb = length(b);
    N = size(x,1);%the size of measurements
    % build X and A matrices using 3rd order moments
    X = zeros(N,4+nb*6); 
    row_ind = 1;
    I = eye(4+nb*6);
    A = zeros(2+nb*2,4+nb*6);
    A(1:2,:) = I(2:3,:);% r, vr>0

    for i = 1:nb
        ind = find(x(:,3)==b(i));
        L = length(ind);
        %xi: 1 t,b,t^2,b^2,tb,t^3,b^3,t^2b,tb^2
        Xi = cell2mat(cellfun(@vm_1x2_to_1x10,mat2cell(x(ind,1:2),ones(L,1),2), 'UniformOutput',false));
        
        X(row_ind:row_ind+L-1,:) =[Xi(:,[1 2 4 7]) zeros(L,6*(i-1)) Xi(:,[3 5 6 8 9 10]) zeros(L,6*(nb-i))];
        row_ind = row_ind+L;
        A((2*i+1):(2*i+2),:) = [zeros(2,4+6*(i-1)) eye(2,6) zeros(2,6*(nb-i))]; %d vd>0
    end
    % remove elements in X and A with higher orders
    % order of t and d in 
    %Xi(:,[3 5 6 8 9 10])=b,b^2,tb,b^3,t^2b tb^2
    order_t=[0 0 1 0 2 1];
    order_d=[1 2 1 3 1 2];
    id=(order_t<=order(1))&(order_d<=order(2))&((order_t+order_d)<=max(order));
    %the columns of X and A to keep
    id_column = [ones(1,order(1)+1) zeros(1,3-order(1)) repmat(id,[1 nb])];
    X=X(:,id_column);
    A=A(:,id_column);
    
end


function theta_full=ConvertFullParam(theta,x,order)
b = unique(x(:,3)); % unique set of b-vectors
nb = length(b);
theta_full = zeros(1,4+nb*6);
order_t=[0 0 1 0 2 1];
order_d=[1 2 1 3 1 2];
id=(order_t<=order(1))&(order_d<=order(2))&((order_t+order_d)<=max(order));
id_column = [ones(1,order(1)+1) zeros(1,3-order(1)) repmat(id,[1 nb])];
theta_full(id_column) = theta;
end
