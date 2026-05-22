function [Spec,f]=acf2spec_wlsq(acf,acf_var,lag,n_freqs,alpha_Tikhonov)
% acf2spec_wlsq ACF-to-Power-spectra for irregular samplings
% This function accurately transforma an irregularly sampled ACF to
% the corresponding power-spectra using explcitit calculation of
% the DFT-matrix for irregularly sampled ACF to a power-spectra on
% a regular frequency-grid.
% 
% The FFT-like transform ("like" since irragularly sampled) is done
% with SVD and feather-light 0th-order Tikhonov regularization.
% 
% Calling:
%  [Spec,f]=acf2spec_wlsq(acf,acf_var,lag,n_freq,alpha_Tikhonov)
% Input:
%  acf
%  acf_var
%  lag
%  n_freq
%  alpha_Tikhonov
% 
%  Copyright © Gustavsson and Vierinen 20190909,
%  bjorn.gustavsson@uit.no, juha.vierinen@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

persistent U V invZ prev_lags % To facilitate speed-up
Tlsq_inv = 1;

% If no acf-variance is given use a constant variance
if isempty(acf_var)
  acf_var = ones(size(acf));
end

% Determine time-resolution
dt = median(diff(lag));

% If no specific number of frequencies are requested use the
% theoretical maximum - SVD and Tikhonov-regularization should
% automatically determine the intrinsic spectral resolution based
% on the available lags and lag-variance.
if nargin < 4 || isempty(n_freqs)
  if dt < 1e-6 % This means we're looking at plasma-line-wide spectra
    n_freqs = numel(acf); % So half the theoretical spectral
                          % resolution is reasonable.
  else
    n_freqs = numel(acf)*2; % for the ionline we need to fly closer
                            % to the sun
    Tlsq_inv = 1;
  end
end

% This is a good default value for Tikhonov, if modifying this be
% sensible and keep the value between 1e-1 and 100. An alpha of 1
% is a guarantee that the conversion does not amplify noise. See
% for example Tarantola.
if nargin < 5 || iempty(alpha_Tikhonov)
  alpha_Tikhonov = 1;
end

% Ensure that the number of frequencies are odd, a bit superflous
n_freqs = floor(n_freqs/2)*2-1;
% Neat frequency array, ensure that a theory matrix for an evenly
% sampled ACF would be orthogonal.
f = fftfreq(dt,n_freqs);

% Re-Normalize acf-estimates (assuming they are uncorrelated)
acf = acf./sqrt(acf_var);
% Create theory-matrix re-normalized with the measurement variances
M = [cos(2*pi*lag(:)*f);sin(2*pi*lag(:)*f)]./repmat(sqrt(acf_var(:)),2,numel(f));
% These two steps gives the correctly weighted
% least-square-estimate of the power-spectra. The theory-matrix
% expresses the discretized Fourier-transform-relation between the
% real-valued power-spectrum and the complex-conjugate-symmetric
% ACF, for an ACF at unevenly sampled lags.
%
% Since we cash the SVD-decomposition of M it might be preferable
% to use an unweigthed least-square estimate. This assumes that the
% acf_var does not vary signifficantly from one range to the next
% over the ranges where the cashing works (i.e. identical or
% near-identical lags)
try
  if Tlsq_inv
    % For the robust inverse use SVD and 0th-order Tikhonov inverse.
    if numel(lag)==numel(prev_lags) && max(abs(lag-prev_lags))<eps(max(lag))
      % Just use previous SVD-inverse
    else
      [U,S,V] = svd(M);
      LAMBDA = diag(S);
      % The Tikhonov-regularization is done by damping - that only
      % significantly affects singular-values smaller that
      % alpha_Tikhonov.
      invZ = diag(LAMBDA./(LAMBDA.^2 + alpha_Tikhonov));
      prev_lags = lag;
    end
    Spec = V(:,1:size(invZ,2))*invZ*U(:,1:size(invZ,1))'*[real(acf(:));imag(acf(:))];
  else
    Spec = M\[real(acf(:));imag(acf(:))];
  end
catch
  % Tooo many comments?
  Spec = nan(size(f));
end

end



function f = fftfreq(dt,n)
% FFTFREQ - fft-frequency-array
%   
% Calling:
%   f = fftfreq(dt,n)
% Input:
%   dt - sampling duration (s), double scalar
%   n  - number of frequencies,  scalar int
% Output:
%   f - frequency array (Hz), dobule array

if rem(n,2) == 0
  f = [0:n/2-1,-n/2:-1] / (dt*n); %   if n is even
else
  f = [0:(n-1)/2,-(n-1)/2:-1] / (dt*n); %   if n is odd
end

f = fftshift(f);

end

