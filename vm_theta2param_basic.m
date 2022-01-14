function redim = vm_theta2param_basic(theta)
% function drp = vm_theta_to_drp(theta_1x9)
% convert the cumulant expansion parameters to diffusion-relaxation
% properties
% mean relaxation and diffusivity

redim.s0 = exp(theta(1));
redim.r = theta(2); 
redim.d = theta(3);
%center moments
redim.r2 = theta(4);
redim.d2 = theta(5);
redim.rd = theta(6)/sqrt(2);

end


