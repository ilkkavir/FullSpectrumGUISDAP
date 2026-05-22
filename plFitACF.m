
% integrated plasma line ACF data (set display_spectra=-1 in guisdap)
acfdir0 = '2007-06-17_ipy_00_60@32p';
acfdir1 = '2007-06-17_ipy_01_60@32p';

% list the data files
ff0 = dir(fullfile(acfdir0,'*.mat'));
ff1 = dir(fullfile(acfdir1,'*.mat'));

nfile = length(ff0);

time_start = NaT;
time_end = NaT;
DF_fit_all = {};
peak_alt_all = {};
peak_amp_all = {};
peak_frq_fit_all = {};
plcount = 0;

plcols = [];
plffit = [];
for ifile=1:nfile
    load(fullfile(ff0(ifile).folder,ff0(ifile).name));
    dd1 = load(fullfile(ff1(ifile).folder,ff0(ifile).name));
    r_acf = conj(r_acf) + (dd1.r_acf);
    nh = size(r_acf,2);
    for ih=1:nh
        % ACF to spectra using GUISDAP acf2spec_wlsq
        [spec,f] = acf2spec_wlsq(r_acf(:,ih),[],r_lag(:,ih));
        nfreq = length(f);
        if ih==1
            r_spec = NaN(nfreq,nh);
        end
        r_spec(:,ih) = spec - median(spec);
    end
    for ih=1:nh
        r_spec(:,ih) = r_spec(:,ih) - median(r_spec(:,ih));
    end
    % the edge-detection algorithm written with help of gemini AI
    [cols,srow] = plfind_geminiAI(r_spec,0);
    plcols(:,ifile) = cols;
    
    subplot(1,2,1)
    rhd = r_h - median(diff(r_h))/2;
    pcolor(f,rhd,r_spec');shading flat;hold on
    plot(f(cols),r_h)
    colorbar
    hold off
    
    drawnow


    % I was initially fitting the ACF data, but it turns out to be too tricky. The fit below is for the spectra. 
    
    spec_final = r_spec.*NaN;
    % Xfit = [1,1e4,0];
    ffit = [];
    sfit = [];
    afit = [];
    for ih=1:nh
        %[specfit,Xfit] = acf2spec_GaussFit(r_acf(:,ih),r_lag(:,ih),f,[log10(mean(abs(r_acf(1:10,ih)))),Xfit(1)],[1e4 , Xfit(2)],[f(cols(ih)) , Xfit(3)]);
        %[specfit,Xfit] = acf2spec_GaussFit(r_acf(:,ih),r_lag(:,ih),f,[log10(mean(abs(r_acf(1:10,ih))))],[1e4],[f(cols(ih))]);
        %[specfit,Xfit] = acf2spec_GaussFit(r_acf(:,ih),r_lag(:,ih),f,f(cols(ih)),10);

        % fit Gaussian spectra at all altitudes
        [specfit,Xfit] = spec_GaussFit( r_spec(:,ih) , f , log10(max(r_spec(:,ih))) , 5e4 , f(cols(ih)) );

        spec_final(:,ih) = specfit;
        ffit(ih) = Xfit(3);
        sfit(ih) = Xfit(2);
        afit(ih) = Xfit(1);
        % if any(isnan(Xfit))
        %     Xfit = [1,1e4,0];
        % end
    end

    % if the fitted peak frequency is far from that given by the edge detection, there was probably no plasma line
    rminds = abs(ffit - f(cols)) >1e5;%> 3*sfit;% | sfit > 1e5 | sfit < 2e3;
    rminds = rminds | isnan(ffit) | isnan(sfit);
    spec_final(:,rminds) = NaN;
    ffit(rminds) = NaN;
    sfit(rminds) = NaN;
    afit(rminds) = NaN;

    % make a plasma line data struct that can be used by Tikemani's prior model
    if ~all(rminds)
        plcount = plcount + 1;
        time_start(plcount) = datetime(r_time(1,:));
        time_end(plcount) = datetime(r_time(2,:));
        DF_fit_all{plcount} = sfit(~rminds);
        peak_alt_all{plcount} = r_h(~rminds);
        peak_amp_all{plcount} = afit(~rminds);
        peak_frq_fit_all{plcount} = ffit(~rminds);
    end
    

    % draw the plasma line frequencies. 
    plffit(:,ifile) = ffit;
    subplot(1,2,2)
    pcolor(f,rhd,spec_final');shading flat;
    hold on
    plot(ffit,r_h);
    hold off
    colorbar
    drawnow
    %pause
end


outfile = ['Plasma_line_fits_',acfdir0,acfdir1,'.mat'];
save(outfile,'','DF_fit_all','peak_alt_all','peak_amp_all','peak_frq_fit_all','time_end','time_start');

figure
pcolor(1:nfile,rhd,plffit);shading flat;caxis(quantile(plffit(:),[.01,.99])); colorbar;colormap turbo
