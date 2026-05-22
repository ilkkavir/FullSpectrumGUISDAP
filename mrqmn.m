% mrqmn: Least squares fit routine (Marquart-Levenberg method)
% GUISDAP v.8.50 07-01-14 Copyright EISCAT, Huuskonen&Lehtinen
%
% input arguments:
% a        : parameters to fit
% ym       : measurements
% variance : data variance, either diagonal of full
% ftol     : stop execution, when step for all paramters less than ftol
% itmax    : maximum number of iterations
% kd2 p_coeffg, f_womega, p_om, pldfvv : parameters need by dirthe
% Output parameters
% aa    : The least-squares point
% chi2  : residual of the fit
% its   : number of iterations
% alpha : the error covariance matrix at the second last iteration point
%
% See also: dirthe
%  function [aa,chi2,its,alpha]=mrqmn(a,ym,variance,ftol,itmax,kd2,p_coeffg,f_womega,p_om,pldfvv,p_m0)
function [aa,chi2,its,alpha]=mrqmn(a,ym,variance,ftol,itmax,kd2,p_coeffg,f_womega,p_om,pldfvv,p_m0,physlim,fb_womega,fitprop)
    global local path_GUP
    global p_D0 p_N0 p_T0 p_om0 k_radar % IV 2026
    global v_Boltzmann v_epsilon0 v_elemcharge % IV 2026

    persistent printmessage
    
    diagonal=min(size(variance))==1;

    % force the matlab version to get the modified matlab dirthe-function
    if false
        %if diagonal & libisloaded('libguisdap')
        
        [ns,aaN]=size(a);
        [varianceM,varianceN]=size(variance);
        [coefM,coefN]=size(p_coeffg);
        [womM,womN]=size(f_womega);
        nom=prod(size(p_om));
        nion=prod(size(p_m0));
        ymPr=libpointer('doublePtr',ym);
        variancePr=libpointer('doublePtr',variance);
        p_coeffgPr=libpointer('doublePtr',p_coeffg);
        f_womegaPr=libpointer('doublePtr',f_womega);
        fb_womegaPr=libpointer('doublePtr',fb_womega);
        p_omPr=libpointer('doublePtr',p_om);
        
        aaOutPr=libpointer('doublePtr',zeros(ns,aaN));
        chi2Pr=libpointer('doublePtr',1);
        itsPr=libpointer('doublePtr',1);
        alphaPr=libpointer('doublePtr',zeros(aaN,aaN));
        
        calllib('libguisdap','MrqmndiagCalc',ns,aaN,a,ymPr,variancePr,varianceM,varianceN,ftol,itmax, ...
		p_coeffgPr,coefM,coefN,womM,f_womegaPr,kd2,nom,p_omPr,aaOutPr,chi2Pr,itsPr, ...
		alphaPr,real(pldfvv),imag(pldfvv),physlim,p_m0,nion,fb_womegaPr,fitprop);
        
        aa=aaOutPr.value;
        chi2=chi2Pr.value;
        its=itsPr.value;
        alpha=alphaPr.value;
        
    else

        if isempty(printmessage)
            printmessage = 0;
        end
        printmessage = printmessage + 1;
        if printmessage == 100
            printmessage = 0;
            disp('   ')
            disp('MATLAB mrqmn forced and dirthe modified for plasma line input -- IV 2026')
            disp('   ')
        end
        
        if diagonal
            vv=find(variance~=0);
            va=find(variance((end-length(a)+1):end)~=0);
            invsigma=1 ./variance(vv);
            apu=invsigma*ones(1,length(va));
        else
            diagvar=diag(variance);
            vv=find(diagvar~=0);
            va=find(diagvar((end-length(a)+1):end)~=0);
            % The apriori variance values for the parameters (especially for collision frequency)
            % are sometimes so low in scaled units, that a simple 
            % invsigma=inv(variance)
            % produces unnecessary warnings about the condition number of the matrix!
            % Therefore we treat the apriori part separately here.
            %invsigma=diag(1./diag(variance));
            %[M,N]=size(variance);
            %NP=length(a);
            %invsigma(1:M-NP,1:M-NP)=inv(variance(1:M-NP,1:M-NP));
            invsigma=inv(variance(vv,vv));
        end
        lambda=0.001;
        aa=a; validder=0;
        ya=dirthe(aa,p_coeffg,f_womega,kd2,p_om,pldfvv,p_m0,fb_womega);
        dyda=zeros(length(ya),length(va));
        if diagonal
            chi2=(ym(vv)-ya(vv))'*(invsigma.*(ym(vv)-ya(vv)));
        else
            chi2=(ym(vv)-ya(vv))'*invsigma*(ym(vv)-ya(vv));
        end
        its=0;
        % fit logs
        ag=find(fitprop==1); as=find(fitprop==2); %log!
        aa(ag)=log(aa(ag)); aa(as)=asinh(aa(as)/2); %log!
        physlim(:,ag)=log(physlim(:,ag)); physlim(:,as)=asinh(physlim(:,as)/2); %log!
        while its<itmax
            if ~validder
                its=its+1;
                for i=1:length(va)
                    aa2=aa; aa2(va(i))=aa(va(i))+0.0001;
                    aa2(ag)=exp(aa2(ag)); aa2(as)=2*sinh(aa2(as)); %log!
                    dyda(:,i)=(ya-dirthe(aa2,p_coeffg,f_womega,kd2,p_om,pldfvv,p_m0,fb_womega))/0.0001; %calculate derivatives;
                end
                if diagonal
                    alpha=dyda(vv,:)'*(apu.*dyda(vv,:));
                    beta=dyda(vv,:)'*(invsigma.*(ya(vv)-ym(vv)));
                else
                    alpha=dyda(vv,:)'*invsigma*dyda(vv,:);
                    beta=dyda(vv,:)'*invsigma*(ya(vv)-ym(vv));
                end
            end
            da=((alpha+lambda*eye(size(alpha)))\beta)';
            % disp(da);
            chi2old=chi2; yaold=ya; aaold=aa;
            aa(va)=aa(va)+da;
            d=find(aa(va)<physlim(1,va)); aa(va(d))=physlim(1,va(d));
            d=find(aa(va)>physlim(2,va)); aa(va(d))=physlim(2,va(d));
            % da=aa(va)-aaold(va);
            aa2=aa; aa2(ag)=exp(aa2(ag)); aa2(as)=2*sinh(aa2(as)); %log!
            ya=dirthe(aa2,p_coeffg,f_womega,kd2,p_om,pldfvv,p_m0,fb_womega);
            if diagonal
                chi2=(ym(vv)-ya(vv))'*(invsigma.*(ym(vv)-ya(vv)));
            else
                chi2=(ym(vv)-ya(vv))'*invsigma*(ym(vv)-ya(vv));
            end
            if chi2>=chi2old
                chi2=chi2old; aa=aaold; ya=yaold; lambda=lambda*10; validder=1;
                if max(abs(da))<ftol, OK=0; break, end
            else
                lambda=lambda/10; validder=0;
                if max(abs(da))<ftol, OK=0; break, end
            end
        end
        
        al=inv(alpha);
        alpha=zeros(length(aa)); alpha(va,va)=al;
        aa(ag)=exp(aa(ag)); aa(as)=2*sinh(aa(as)); %log
        aa2=ones(length(aa),1); aa2(ag)=aa(ag); aa2(as)=sqrt(4+aa(as).^2); %log
        alpha=alpha.*(aa2*aa2');
        
    end
    
end

