% Understading the parameters
% The algorithm runs in two steps: (1) PSF (h) estimation and (2) final
% image (u) non-blind deconvolution
%
% Step (1) is iterative (with maxiter iterations) and alternates between 
% min_u and min_h. It uses only a central section of the input images (maxROIsize)
% and can be multiscale (MSlevels)
%
% min_u and min_h are iterative with maximum of maxiter_u and maxiter_h
% iterations or stop earlier if the relative change is less than ccreltol
%
% The optimization tasks min_u and min_h are:
% min_u gamma/2 | Hu - g |^2 + alpha_u*PHI(vx,vy) + beta_u/2*| Dx u - vx -
% ax|^2 + beta_u/2*|Dy u - vy - ay |^2
%
% min_h gamma/2 | Uh - g |^2 + alpha_h*R(w) +
% beta_h/2*| h - w - b |^2
%
% Note that gamma depends on noise: 10dB -> 1e1, 20dB -> 1e2, etc.
% However in step (1) I always set gamma=1e1 (or 1e2) and in step (2) (see below) 
% I use the gamma corresponding to the actual noise level.
%
% Step (2) uses estimated PSF from step (1) to perform deconvolution on 
% the whole image. It performs min_u with parameters gamma_nonblind and
% beta_u_nonblind instead of gamma and beta_u.

% PSF estimation parameters
% size of the image central section, where PSF estimation will be calculated 
maxROIsize = [1024 1024];
% multiscale levels {1,2,..} (if =1 do not use hierachical approach) 
MSlevels = 4;

% parameters
PAR.verbose = 2; %{0 = no messages,1 = text messages,2 = text and graphs}

% common parameters to both PSF and image estimation
% data term weight gamma
PAR.gamma = 1e2;
% which Lp norm to use
PAR.Lp = 0.3;

% PSFs estimation
PAR.beta_h = 1e4*PAR.gamma;
PAR.alpha_h = 1e1*PAR.gamma;
PAR.centering_threshold = 20/255; % threshold used for PSF centering between iterations, <=0 for no centering

% image estimation
PAR.beta_u = 1e0*PAR.gamma;
PAR.alpha_u = 1e-2*PAR.gamma; 

% non-blind image estimation (final u-step)
% gamma_nonblind will be used at the end by the final u-step on 
% the full image (optional)
PAR.gamma_nonblind = 2e1*PAR.gamma;
PAR.beta_u_nonblind = 1e-2*PAR.gamma_nonblind;
PAR.Lp_nonblind = 1;

% number of iterations
PAR.maxiter_u = 10;
PAR.maxiter_h = 10;
PAR.maxiter = 10;
PAR.ccreltol = 1e-3;