function plotResultFile(fname,nlim,tlim,vlim,olim,hlim,trmax)

    if nargin < 7
        trmax = 3;
    end
    
    if nargin < 6
        hlim = [];
    end
    if nargin < 5
        olim = [-.1 1.1];
    end
    if nargin < 4
        vlim = [-1 1]*200;
    end
    if nargin < 3
        tlim = [0 3000];
    end
    if nargin < 2
        nlim = [1e10,1e12];
    end
    if isempty(olim)
        olim = [-.1 1.1];
    end
    if isempty(vlim)
        vlim = [-1 1]*200;
    end
    if isempty(tlim)
        tlim = [0 3000];
    end
    if isempty(nlim)
        nlim = [1e10,1e12];
    end

    load(fname);
    %    figure;
    figpos = get(gcf,'position');
    set(gcf,'position',[figpos(1:2), 900 900])

    % r_param(r_status~=0,:) = NaN;
    % r_param(r_res(:,1)>10,:) = NaN;
    % r_param(r_param(:,3)>trmax,:) = NaN;


    if exist('r_param_filter')
        r_param = r_param_filter;
        r_error = r_error_filter;
    end
    splot = false;
    %    if exist('r_param_rcorr_smooth')
        %r_param_rsmooth = r_param_rcorr_smooth;
        %r_error_rsmooth = r_error_rcorr_smooth;
    if exist('r_apriori')
        r_param_rsmooth = r_apriori;
        r_error_rsmooth = r_apriorierror;
        splot = true;
    end
        
    
    subplot(2,2,1)
    errorbar(r_param(:,1),r_h,r_error(:,1),'horizontal','k.')
    if splot
        hold on
        errorbar(r_param_rsmooth(:,1),r_h,r_error_rsmooth(:,1),'horizontal','m-')
        hold off
    end
    xlim(nlim)
    if ~isempty(hlim)
        ylim(hlim)
    end
    xlabel('N_e [m^{-3}]')
    ylabel('Altitude [km]')
    set(gca,'xscale','log')

    subplot(2,2,2)
    errorbar(r_param(:,2).*r_param(:,3),r_h, sqrt((r_error(:,2).^2).*r_param(:,3)+(r_error(:,3).^2).*r_param(:,2)),'horizontal','r.')
    hold on
    errorbar(r_param(:,2),r_h,r_error(:,2),'horizontal','b.')
    if splot
        errorbar(r_param_rsmooth(:,2),r_h,r_error_rsmooth(:,2),'horizontal','m-')
        errorbar(r_param_rsmooth(:,2).*r_param_rsmooth(:,3),r_h,sqrt(r_error_rsmooth(:,2).^2 .* r_param(:,3) + r_error_rsmooth(:,3).^2.*r_param(:,2)),'horizontal','m-')
    end
    xlim(tlim)
    if ~isempty(hlim)
        ylim(hlim)
    end
    xlabel('Temperature [K]')
    ylabel('Altitude [km]')
    legend('T_e','T_i','location','southeast')
    %    legend('boxoff')
    hold off

    subplot(2,2,3)
    errorbar(r_param(:,5),r_h,r_error(:,5),'horizontal','k.')
    if splot
        hold on
        errorbar(r_param_rsmooth(:,5),r_h,r_error_rsmooth(:,5),'horizontal','m-')
        hold off
    end
    xlim(vlim)
    if ~isempty(hlim)
        ylim(hlim)
    end
    ylabel('Altitude [km]')
    xlabel('V_i [ms^{-1}]')
    
    subplot(2,2,4)
    errorbar(r_param(:,6),r_h,r_error(:,6),'horizontal','k.')
    if splot
        hold on
        errorbar(r_param_rsmooth(:,6),r_h,r_error_rsmooth(:,6),'horizontal','m-')
        hold off
    end
    xlim(olim)
    if ~isempty(hlim)
        ylim(hlim)
    end
    ylabel('Altitude [km]')
    xlabel('[O^+]')

    sgtitle(datestr(datetime(r_time(1,:))+(datetime(r_time(2,:))-datetime(r_time(1,:)))/2))

end
