function Theta = REDIM_Data2Fit(S,x,mask,method,order)
% function Theta = REDIM_Data2Fit(S,x,mask,method,order)
% input: S: combined dMRI data with multiple TE
%        x: Nb*L x 3 matrix, x(:,1): TE values, x(:,2) b-value , x(:,3)
%        index of bvec
%       mask: brain mask
%       method: fitting method, (default 'wn', weighted nonlinear)
%       order: order of cumulant expansion used for fitting (default [1 1 0 1 1 0 0 0 0]);
% output: Theta: fitting paramters for REDIM model

if(~exist('order','var'))
        order = [1 1 0 1 1 0 0 0 0];%[r d r2 d2 rd r3 d3 r2d rd2]
end

if(~exist('method','var'))
        method = 'wn';%[r d r2 d2 rd r3 d3 r2d rd2]
end 


[nx,ny,nz,Nb,L] = size(S);
Nu = length(unique(x(:,3)))-1; % number of gradient directions
S = reshape(S,nx*ny*nz,Nb*L);
N_param = 4+6*Nu;

id_mask = find(mask>0);
S_mask = S(id_mask,:);
Theta_mask = zeros(length(id_mask),N_param);
% 
% construct measurement and constraint matrices
x(:,2) = x(:,2)/1000;% change unit to ms/um^2
[X,A,b] = ParamToMatrix(x, order);



parfor i = 1:length(id_mask)

        %theta = redim_data2fit_1vox_v3(S(i,:),x,order);
       %try
           theta = redim_data2fit_1vox(S_mask(i,:),X,A,b,method);
           theta = ConvertFullParam(theta,Nu,order);
           if(~isnan(theta))
               Theta_mask(i,:) = theta;
           end
       %catch
       %    fprintf('error index %d \n',i);
       %end

end

Theta = zeros(nx*ny*nz,N_param);
Theta(id_mask,:) = Theta_mask;
Theta = reshape(Theta,nx,ny,nz,N_param);
end



function theta = redim_data2fit_1vox(signal,X,A,b,method)
% function theta = redim_data2fit_1vox(signal, x,order)
%
% Diffusion-Relaxation coupling estimation. 
% Input: signal is an 1xN (or Nx1) vector that is measured using different
% TE and b-values along different direction. 
% X: measurement matrix X*theta = signal
% A,b: contraints A*theta<= b
% method: 'l' linear,'n' nonlinear,'wl' weighted linear,'wn' weighted
% nonlinear
% Output: estimation parameters

signal = signal(:);
% only use positive signals for fitting
X = X(signal>0,:);
signal = signal(signal>0); 

if(length(signal)>0)

    switch method
        case 'l'
            options = optimoptions('lsqlin','display', 'off','Algorithm','interior-point'); 
            theta = lsqlin(X,log(signal),A,b,[],[],[],[],[],options);
        case 'wl'
                options = optimoptions('lsqlin','display', 'off','Algorithm','interior-point'); 
                w = signal(:)/max(signal);
                theta = lsqlin(diag(sqrt(w))*X,sqrt(w).*log(signal),A,b,[],[],[],[],[],options);
        case 'n'
                % initialize use 'wl' based f
                options = optimoptions('lsqlin','display', 'off','Algorithm','interior-point'); 
                w = signal(:)/max(signal);
                theta0 = lsqlin(diag(sqrt(w))*X,sqrt(w).*log(signal),A,b,[],[],[],[],[],options);

                fun = @(x)sum((exp(X*x)-signal).^2);
                options2 = optimoptions('fmincon','display', 'off'); 
                if(~isnan(theta0))    
                    theta = fmincon(fun,theta0,A,b,[],[],[],[],[],options2);
                else
                    theta= theta0;
                end
         case 'wn'
                % initialize use 'wl' based f
                options = optimoptions('lsqlin','display', 'off','Algorithm','interior-point'); 
                w = signal(:)/max(signal);
                theta0 = lsqlin(diag(sqrt(w))*X,sqrt(w).*log(signal),A,b,[],[],[],[],[],options);
                fun = @(x)sum(w.*(exp(X*x)-signal).^2);
                options2 = optimoptions('fmincon','display', 'off'); 
                if(~isnan(theta0))    
                    theta = fmincon(fun,theta0,A,b,[],[],[],[],[],options2);
                else
                    theta= theta0;
                end
          otherwise
                error('incorrect method');
    end

    theta = theta';
    % reformat theta to a 1x10 vector
else
    theta = nan(1,size(X,2));
end


end



function [X,A,b] = ParamToMatrix(x, order)
% build measurement and constriant matrix
% X*theta = y;
% A*theta <= b;
    if(nargin==1)
        order = [1 1 0 1 1 0 0 0 0];%[r d r2 d2 rd r3 d3 r2d rd2]
    end

    bvec = unique(x(x(:,3)>0,3)); % unique set of b-vectors
    nb = length(bvec); % nonzero gradient
    N = size(x,1);%the size of measurements
    % build (measurement matrix) X and (constraint) A matrices using 3rd order moments
    X = zeros(N,4+nb*6); 

    I = eye(4+nb*6);
    A = zeros(2+nb*3,4+nb*6);
    A(1:2,:) = I(2:3,:);% r, vr>0
    % build measurement matrix
    for i = 1:N
        xi = vm_1x2_to_1x10(x(i,1:2));
        if(x(i,3)==0) % zero gradient
            X(i,:) = [xi([1 2 4 7]) zeros(1,nb*6)];
        else
            id_u = x(i,3);
            X(i,:) = [xi([1 2 4 7]) zeros(1,6*(id_u-1)) xi(1,[3 5 6 8 9 10]) zeros(1,6*(nb-id_u))];
        end
    end
    % build constraint matrix
    for i = 1:nb
        A((3*i):(3*i+1),:) = [zeros(2,4+6*(i-1)) eye(2,6) zeros(2,6*(nb-i))]; %d vd>0
        A((3*i+2),:) = [zeros(1,4+6*(i-1)) 1 -3 0 0 0 0 zeros(1,6*(nb-i))]; %d-bmax*vd>0 monotonical decay
    end
    % remove elements in X and A with higher orders
    % order of t and d in 
    %Xi(:,[3 5 6 8 9 10])=b,b^2,tb,b^3,t^2b tb^2
    order_r = order([1 3 6]);
    order_rd = order([2 4 5 7 8 9]);
    %the columns of X and A to keep
    id_column = [1 order_r repmat(order_rd,[1 nb])];
    X=X(:,id_column>0);
    A=A(:,id_column>0);
    A = A(sum(abs(A),2)>0,:);%remove zero rows
    
    b = zeros(size(A,1),1);%r, vr,d,vd >0
    B = A((sum(order_r)+1):3:end,:); % diffusivity smaller than 3
    A = [-A; B];
    b = [b; 3*ones(size(B,1),1)];% constrain the upper bound of diffusivity

end


function theta_full=ConvertFullParam(theta,nu,order)

theta_full = zeros(1,4+nu*6);
order_r = order([1 3 6]);
order_rd = order([2 4 5 7 8 9]);
%the columns of X and A to keep
id_column = [1 order_r repmat(order_rd,[1 nu])];
theta_full(id_column>0) = theta;
end
