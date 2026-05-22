function [Spec,Xfit]=acf2spec_GaussFit(acf,lag,f,amp0,std0,freq0)
% acf2spec_GaussFit ACF-to-Power-spectra by means of fitting a Gaussin ACF.
% 
% Calling:
%  [Spec,f]=acf2spec_GaussFit(acf,lag,f,freq0)
% Input:
%  acf
%  lag
%  f
%  freq0
% 
%  Copyright © Ilkka Virtanen 2026,
%  Ilkka.i.virtanen@oulu.fi
%  This is free software, licensed under GNU GPL version 2 or later

persistent Xprev

% separate real and imaginary parts
riacf = [real(acf(:));imag(acf(:))];

ntestfreqs = length(freq0);

Xfitn = {};
fval = [];
exitflag = [];
output = {};
X0n = {};

for itest = 1:ntestfreqs
    
    % Initial values for the fitted parameters.
    % X(1) = power
    % X(2) = standard deviation (width)
    % X(3) = Doppler shift (peak frequency)
    X0 = [amp0(itest),std0(itest),freq0(itest)];
    
    Xscale = [1 , 1e5 , 1e6];

    X0n{itest} = X0./Xscale;

    % fms_opts = optimset('fminsearch');
    % fms_opts.Display = 'off';%'off';%'final';
    % fms_opts.MaxFunEvals=1e4;
    % fms_opts.MaxIter=1e6;
    % fms_opts.TolFun=1e-10;
    % fms_opts.TolX=1e-10;
    
    
    %[Xfitn{itest},fval(itest),exitflag(itest),output{itest}] = fminsearch( @(X) costGauss(X,lag,riacf,Xscale) , X0n{itest} , fms_opts);


    %options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
    % Set custom stopping criteria
    options = optimoptions('lsqnonlin', ...
                           'Algorithm','levenberg-marquardt',...
                           'MaxIterations', 10000, ...      % Default is often 400
                           'StepTolerance', 1e-15, ...      % Tighten step requirement
                           'FunctionTolerance', 1e-10,...      % Stop if function stops improving
                           'optimalityTolerance',1e-12,...
                           'display','none');
    
    [Xfitn{itest},~,fval(itest),exitflag(itest),output{itest}] = lsqnonlin(@(X) costGauss(X,lag,riacf,Xscale), X0n{itest}, [], [], options);

    
%     plot(riacf)
%     hold on
%     X0fn = X0n{itest}
%     Xfn = Xfitn{itest}
% exitflag
%     plot(GaussModelAcf(Xfn(1)*Xscale(1),Xfn(2)*Xscale(2),Xfn(3)*Xscale(3),lag))
%     plot(GaussModelAcf(X0fn(1)*Xscale(1),X0fn(2)*Xscale(2),X0fn(3)*Xscale(3),lag))
%     hold off
%     drawnow
%     pause(1)

end
fval(exitflag~=1) = Inf;
[mincost,imincost] = min(fval);

if all(exitflag~=1)
    Xfit = [NaN,NaN,NaN];
else
    
    %imincost = 1;
    
    %Xfitn = X0n;
    
    Xfit = Xfitn{imincost}.*Xscale;
    X0 = X0n{imincost}.*Xscale;
    %Xfit = X0;
    %Spec = (10^Xfit(1))*exp( - ( Xfit(2)^2 * (2*pi*(f'-Xfit(3))).^2 ) / 2 );
end
% MATLAB code for the Power Spectral Density (PSD)
Spec = (10.^Xfit(1) / (Xfit(2) * sqrt(2*pi))) .* exp(-(f' - Xfit(3)).^2 / (2 * Xfit(2)^2));



% figure(19)
% plot(riacf)
% hold on
% plot(GaussModelAcf(Xfit(1),Xfit(2),Xfit(3),lag));
% hold off
% pause

end


function riacfGauss = GaussModelAcf(amp,sigma,doppler,lag)
%disp([amp,sigma,doppler])

% sigma is width of the spectrum, not the acf
acfGauss = (10^amp) * exp(-lag.^2 * (2*pi^2*sigma^2)).*exp(1i*2*pi*doppler*lag);
riacfGauss = [real(acfGauss(:));imag(acfGauss(:))];

end


function cost = costGauss(X,lag,riacf,Xscale)
cost = sum( ( riacf - GaussModelAcf(X(1)*Xscale(1),X(2)*Xscale(2),X(3)*Xscale(3),lag) ).^2 );
%disp(cost)
end
