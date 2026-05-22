function plotResultDir(dirname,nlim,tlim,vlim,olim,hlim,trmax)

    if nargin < 7
        trmax = [];
    end
    
    if nargin < 6
        hlim = [];
    end
    if nargin < 5
        olim = [];
    end
    if nargin < 4
        vlim = [];
    end
    if nargin < 3
        tlim = [];
    end
    if nargin < 2
        nlim = [];
    end
    if isempty(olim)
        olim = [];
    end
    if isempty(vlim)
        vlim = [];
    end
    if isempty(tlim)
        tlim = [];
    end
    if isempty(nlim)
        nlim = [];
    end


    ff = dir(fullfile(dirname,'*.mat'));

    outdir = [dirname,'_figures'];

    mkdir(outdir)
    
    figure('visible','off')
    
    for k=1:length(ff)                        
        plotResultFile(fullfile(ff(k).folder,ff(k).name))
        [~,b,~] = fileparts(ff(k).name);
        outfile = fullfile(outdir,b);
        print(outfile,'-dpng')
    end

    close(gcf)
end
