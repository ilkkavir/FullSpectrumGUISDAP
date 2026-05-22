% dirthe.m: DIRect THEory for Inocherent Scatter radar measurements
% GUISDAP v.1.60 96-05-27 Copyright Asko Huuskonen and Markku Lehtinen
% 
% Available also as a fast mex-version
% Input parameters:
% param: plasma parameters in scaled units
% p_coeffg: radar factor for each measurement 
% f_womega: spectral ambiguity function for each measurement
% kd2: (k D)^2, k is radar k-vector and D is Debye length for scale parameters
% p_om: frequency axis for f_womega
% pldfvv: plasma dispersion function interpolation table
% Output parameter:
% theo: theoretical values for the measurements
% See also: spec, transf
function theo=dirthe(param,p_coeffg,f_womega,kd2,p_om,pldfvv,p_m0,fb_womega)
    global local path_GUP
    global k_radar v_epsilon0 v_Boltzmann v_elemcharge p_N0 p_T0
    
    if libisloaded('libguisdap')
        
        if ~isreal(f_womega)
            f_womega=real(f_womega);
        end
        
        nom=prod(size(p_om));
        [womM,womN]=size(f_womega);
        [ns,aaN]=size(param);
        [coefM,coefN]=size(p_coeffg);
        nion=prod(size(p_m0));
        
        acfPr=libpointer('doublePtr',zeros(coefM+aaN,coefN));
        scr=libpointer('doublePtr',zeros(nion+2,3+4*nom));
        scr1=libpointer('doublePtr',zeros((nion+1)*5+nom+aaN,ns));
        womPr=libpointer('doublePtr',f_womega);
        
        calllib('libguisdap','DirtheCalc',ns,aaN,param,p_coeffg,womM,womPr,kd2,nom,p_om, ...
	        real(pldfvv),imag(pldfvv),acfPr,0,p_m0,nion,fb_womega,scr,scr1);
        
        theo=acfPr.value;
        
    else
        
        param=real(param); % hyi hyi
        
        %nonphys=physlim(param,p_m0);
        
        [nin0,tit0,mim0,psi,vi]=transf(param,p_m0);
        s=spec(nin0,tit0,mim0,psi,vi,kd2,p_om,pldfvv);
        %add clutter signal
        b=param(5+length(p_m0)); % Broadband noise not from pulse
        f0=find(p_om==0); %DC spike
        if ~isempty(f0)
            s(f0)=param(6+length(p_m0))+s(f0);
        end
        
        theo=[p_coeffg.*(f_womega*s)+fb_womega*b;col(param)];
        
    end

    % Modify the Ne value to match with the cold plasma approximation used in the plasma line prior. 
    % IV 2026
    npar = length(param);
    A = (3 * k_radar(1)^2 * v_epsilon0 * v_Boltzmann ) / (v_elemcharge^2);
    % add the Te contribution to get plasma line Ne under cold plasma assumption
    theo(end-npar+1) = param(1) + A * param(2)*param(3)*p_T0 / p_N0;

end

