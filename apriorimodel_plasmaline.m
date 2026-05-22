function [apriori2,apriorierror2] = apriorimodel_plasmaline(apriori,apriorierror,heights,fit_altitude)
% [apriori2,apriorierror2] = apriorimodel_plasmaline(apriori,apriorierror,heights,fit_altitude)

    global d_time

    % Physical constants defined in guisdap
    global v_lightspeed v_epsilon0 v_elemcharge v_electronmass v_Boltzmann k_radar k_radar0 p_T0 p_N0

    % convert Ne in apriori into "plasma line Ne under cold plasma approximation"
    A = (3 * k_radar(1)^2 * v_epsilon0 * v_Boltzmann ) / (v_elemcharge^2);
    apriori(:,1) = apriori(:,1) + A*apriori(:,2).*apriori(:,3)*p_T0/p_N0;
    

    
    apriori2 = apriori;
    apriorierror2 = apriorierror;
    
    %% ------------------------------------------------------------------------
    % Load plasma line data
    % -------------------------------------------------------------------------
    %data_path = '/home/tibag_temp/Desktop/IONTRACE/GUISDAP_OUT/plasma_30Jan2022_beata_2_2_00sec_fit_altitude_2022-01-30_mod_new.mat';
    %%%3km gate updated
    %data_path='/home/tibag_temp/Desktop/IONTRACE/GUISDAP_OUT/Plasma_line_30Jan2022_beata_2_2_00sec_2022-01-30_3km_range.mat';
    %data_path='/home/ilkkavir/results/2026-Tikemani-plasmalines/plasmalinecode_20260505/Plasma_line_30Jan2022_beata_2_2_00sec_2022-01-30_3km_range.mat';
    data_path='/home/ilkkavir/results/2026-Tikemani-plasmalines/ipytest202605/Plasma_line_fits_2007-06-17_ipy_00_60@32p2007-06-17_ipy_01_60@32p.mat';
    
    if ~exist(data_path, 'file'), return; end
    data = load(data_path);
    
    %% ------------------------------------------------------------------------
    % Robust HH:MM Matching
    % -------------------------------------------------------------------------
    %
    % Changed this into a more flexible version -- IV
    %
    %
    gv_h = d_time(1,4);
    gv_m = d_time(1,5);
    %
    %mat_vec = datevec(data.time_start); 
    %mat_h = mat_vec(:,4);
    %mat_m = mat_vec(:,5);
    %
    %time_idx = find(mat_h == gv_h & mat_m == gv_m);
    
    ion_start = datetime(d_time(1,:));
    ion_end = datetime(d_time(2,:));
    
    plasma_centre =  data.time_start + .5*(data.time_end-data.time_start);
    
    time_idx = find(plasma_centre > ion_start & plasma_centre < ion_end);
    
    if isempty(time_idx)
        return
    end
    
    %% ------------------------------------------------------------------------
    % Physical constants
    % -------------------------------------------------------------------------
    %c = 299792458; eps0 = 8.854187e-12; e_charge = 1.602176e-19;
    %m_e = 9.109383e-31; k_B = 1.380649e-23;
    %f_radar = 933e6; k_radar = 2*pi*f_radar/c; k_plasma_line = 2*k_radar;
    %
    % use the GUISDAP globals, IV 202605
    
    %% ------------------------------------------------------------------------
    % Extract and Process All Peaks for this Time
    % -------------------------------------------------------------------------
    all_pl_alt = []; all_pl_frq = []; all_pl_df = [];
    for i = 1:length(time_idx)
        c_idx = time_idx(i);
        if iscell(data.peak_alt_all)
            alt_val = data.peak_alt_all{c_idx};
            frq_val = data.peak_frq_fit_all{c_idx};
            df_val  = data.DF_fit_all{c_idx};
        else
            alt_val = data.peak_alt_all(c_idx, :);
            frq_val = data.peak_frq_fit_all(c_idx, :);
            df_val  = data.DF_fit_all(c_idx, :);
        end
        frq_val = frq_val(:);
        df_val = df_val(:);
        valid = ~isnan(alt_val(:)) & ~isnan(frq_val(:)) & (frq_val(:) ~= 0);
        if any(valid)
            all_pl_alt = [all_pl_alt; alt_val(valid(:))];
            %all_pl_frq = [all_pl_frq; abs(frq_val(valid(:))-6.4e6)]; 
            all_pl_frq = [all_pl_frq; abs(frq_val(valid(:))+4.0e6)]; % ipy
            all_pl_df  = [all_pl_df; abs(df_val(valid(:)))];
        end
    end
    
    if isempty(all_pl_alt), return; end
    
    %% ------------------------------------------------------------------------
    % Map to Nearest GUISDAP Gate
    % -------------------------------------------------------------------------
    % Storage for averaging multiple peaks per gate
    gate_ne_accumulator = cell(length(heights), 1);
    gate_err_accumulator = cell(length(heights), 1);
    gate_alt_accumulator = cell(length(heights), 1);
    
    for k = 1:length(all_pl_alt)
        % 1. Find the nearest GUISDAP height index
        [~, nearest_gate_idx] = min(abs(heights - all_pl_alt(k)));
        
        % 2. Calculate Ne for this specific peak
        % Use Te from the nearest gate
        Te_prior = apriori(nearest_gate_idx, 2) * apriori(nearest_gate_idx, 3);
        
        fr = all_pl_frq(k);
        df = all_pl_df(k);

        % swiched to physical constants already defined in GUISDAP
        %    term1 = ( (fr * 2 * pi)^2 * m_e * eps0 ) / (e_charge^2);
        %    term2 = (3 * k_plasma_line^2 * eps0 * k_B * Te_prior) / (e_charge^2);
        term1 = ( (fr * 2 * pi)^2 * v_electronmass * v_epsilon0 ) / (v_elemcharge^2);
        term2 = (3 * k_radar(1)^2 * v_epsilon0 * v_Boltzmann * Te_prior) / (v_elemcharge^2);
        %ne_val = term1 - term2;
        ne_val = term1; % Use the cold plasma approximation. The corresponding changes must be done in mrqmn.m and dirthe.m
        %        ne_err = (8 * pi^2 * m_e * eps0 / e_charge^2) * fr * df; ??
        ne_err = (8 * pi^2 * v_electronmass * v_epsilon0 / v_elemcharge^2) * fr * df;
        
        % 3. Accumulate
        gate_ne_accumulator{nearest_gate_idx} = [gate_ne_accumulator{nearest_gate_idx}; ne_val];
        gate_err_accumulator{nearest_gate_idx} = [gate_err_accumulator{nearest_gate_idx}; ne_err];
        gate_alt_accumulator{nearest_gate_idx} = [gate_alt_accumulator{nearest_gate_idx}; all_pl_alt(k)];
    end
    
    %% ------------------------------------------------------------------------
    % Update Apriori with Averaged Results
    % -------------------------------------------------------------------------
    updated_count = 0;
    nplines = 0;
    for ipl = 1:length(heights)
        if ~isempty(gate_ne_accumulator{ipl})
            nplines = nplines + length(gate_ne_accumulator{ipl});
        end
    end
    if nplines > 4
        for ipl = 1:length(heights)
            if ~isempty(gate_ne_accumulator{ipl})
                if abs(mean(gate_alt_accumulator{ipl}) - heights(ipl)) < 3
                    %apriori2(ipl, 1) = mean(gate_ne_accumulator{ipl});
                    %apriorierror2(ipl, 1) = mean(gate_err_accumulator{ipl});
                    % variance-weighted average
                    apriorierror2(ipl, 1) = sqrt(1./sum(1./gate_err_accumulator{ipl}.^2));
                    apriori2(ipl, 1) = sum(gate_ne_accumulator{ipl}./(gate_err_accumulator{ipl}.^2)).*apriorierror2(ipl,1).^2;
                    apriorierror2(ipl, 6) = 0.5; 
                    disp([heights(ipl) apriori2(ipl,1) apriorierror2(ipl,1)])
                    updated_count = updated_count + 1;
                end
            end
        end
    end
    %% ------------------------------------------------------------------------
    % Save Diagnostics
    % -------------------------------------------------------------------------
    if updated_count > 0
        fprintf('HH:MM MATCH SUCCESS: %02d:%02d | Assigned %d peaks to %d gates.\n', ...
                gv_h, gv_m, length(all_pl_alt), updated_count);
        
        %save_folder = "/home/tibag_temp/Desktop/IONTRACE/GUISDAP_OUT/Plasma_line/check_apriori_00sec_mc120_new_only_fit_alt_new_new/";
        save_folder='/home/ilkkavir/results/2026-Tikemani-plasmalines/plasmalinecode_20260505/check_apriori_00sec_mc120_new_only_fit_alt_new_new/';
        if ~exist(save_folder, 'dir'), mkdir(save_folder); end
        
        guisdap_start = datenum(d_time(1,1:6));
        tstr = datestr(guisdap_start,'yyyymmdd_HHMMSS');
        fname = fullfile(save_folder, "plasmaline_debug_" + tstr + ".mat");
        
        % Save variables for verification
        pl_alt = all_pl_alt;
        pl_frq = all_pl_frq;
        save(fname, "pl_alt", "pl_frq", "heights", "apriori2", "apriorierror2");
    end
    
end
