function [Spec,Xfit]=spec_GaussFit(spec,f,amp0,std0,freq0)
% spec_GaussFit Fit Gaussin spectra to plasma line data
% 
% Calling:
%  [Spec,f]=spec_GaussFit(spec,f,amp0,std0,freq0)
% Input:
%  spec
%  f
%  amp0
%  std0
%  freq0
%  
%
%  Copyright © Ilkka Virtanen 2026,
%  Ilkka.i.virtanen@oulu.fi
%  This is free software, licensed under GNU GPL version 2 or later

%persistent Xprev

% Initial values for the fitted parameters.
% These are given as inputs and are guessed based on the raw spectra and the automatic peak frequency detection
X0 = [amp0,std0,freq0];
    
Xscale = [1 , 1e5 , 1e6];

X0n = X0./Xscale;

fms_opts = optimset('fminsearch');
fms_opts.Display = 'off';%'off';%'final';
fms_opts.MaxFunEvals=1e4;
fms_opts.MaxIter=1e6;
fms_opts.TolFun=1e-10;
fms_opts.TolX=1e-10;
    

[Xfitn,fval,exitflag,output] = fminsearch( @(X) costGauss(X,f(:),spec,Xscale) , X0n , fms_opts);


% % Set custom stopping criteria
% options = optimoptions('lsqnonlin', ...
%                        'Algorithm','levenberg-marquardt',...
%                        'MaxIterations', 10000, ...      % Default is often 400
%                        'StepTolerance', 1e-15, ...      % Tighten step requirement
%                        'FunctionTolerance', 1e-10,...      % Stop if function stops improving
%                        'optimalityTolerance',1e-12,...
%                        'display','none');

% [Xfitn,~,fval,exitflag,output] = lsqnonlin(@(X) costGauss(X,f,spec,Xscale), X0n, [], [], options);


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

%end

if exitflag~=1
    Xfit = [NaN,NaN,NaN];
else
    Xfit = Xfitn.*Xscale;
    X0 = X0n.*Xscale;
end
%Xfit = X0;

% calculate the fitted spetrum. Could be skipped to gain speed. 
Spec = 10.^Xfit(1) .* exp(-(f' - Xfit(3)).^2 / (2 * Xfit(2)^2));



% figure(19)
% plot(riacf)
% hold on
% plot(GaussModelAcf(Xfit(1),Xfit(2),Xfit(3),lag));
% hold off
% pause

end


function specGauss = GaussModelSpectrum(amp,sigma,doppler,f)
%disp([amp,sigma,doppler])
specGauss = (10^amp) * exp(-(f-doppler).^2/(2*sigma^2));

end


function cost = costGauss(X,f,spec,Xscale)
cost = sum( ( spec - GaussModelSpectrum(X(1)*Xscale(1),X(2)*Xscale(2),X(3)*Xscale(3),f) ).^2 );
%disp(cost)
end
