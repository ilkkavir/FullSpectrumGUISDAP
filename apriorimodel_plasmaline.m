function [apriori2,apriorierror2] = apriorimodel_plasmaline(apriori,apriorierror,heights,fit_altitude)
% [apriori2,apriorierror2] = apriorimodel_plasmaline(apriori,apriorierror,heights,fit_altitude)
%
% Ne prior from plasma line data and O+ ion fraction fit.

    % Load the plasma line data
    load("/home/tibag_temp/Desktop/IONTRACE/GUISDAP_OUT/Plasma_line/plasma_30Jan2022_beata_2_2_time_2022-01-30.mat")
    global d_time 
    
    apriori2 = apriori;
    apriorierror2 = apriorierror;

    %% Physical Constants
    c = 299792458;
    eps0 = 8.854187e-12;
    e_charge = 1.602176e-19;
    m_e = 9.109383e-31;
    k_B = 1.380649e-23;
    f_radar = 933e6; 
    k_radar = 2 * pi * f_radar / c;
    k_plasma_line = 2 * k_radar; 

    % Height Gate
    gate_size = 3; % km
    half_gate = gate_size / 2;
    low_e = heights - half_gate;
    upp_e = heights + half_gate;

    % Check if plasma line variables exist in the loaded .mat
    if ~exist('peak_alt_fit_all', 'var'), return; end

    % Convert d_time to datenum for filtering
    guisdap_start = datenum(d_time(1,1:6));
    guisdap_end   = datenum(d_time(2,1:6));

    mat_start = datenum(time_start);
    mat_end   = datenum(time_end);
    
    % Filter plasma line by current integration time (60 second)
    time_idx = find(mat_start >= guisdap_start - 1/60 & mat_end <= guisdap_end + 1/60);

    if isempty(time_idx), return; end
    
    % Get plasma line data for this specific time interval
    pl_alt = peak_alt_fit_all(time_idx);
    pl_frq = peak_frq_fit_all(time_idx);
    pl_df  = DF_fit_all(time_idx);

    % use loop
    num_gates = size(apriori, 1);
    
    for ipl = 1:num_gates
        
        % Find all plasma line peaks that fall within the current ion line gate
        idx_gate = find(pl_alt >= low_e(ipl) & pl_alt < upp_e(ipl));
        
        if ~isempty(idx_gate)
            
            % Column 2: Ti, Column 3: Tr (Te/Ti)
            Te_prior = apriori(ipl, 2) * apriori(ipl, 3);
            
            % Plasma line frequency shift (Hz); 6MHz (frequency band=2) 
            fr = abs(pl_frq(idx_gate) - 6e6);  
            
            % Calculate Ne (m^-3) with Debye correction
            % formula: f_p^2 = f_r^2 - 3*v_th^2*k^2 / (2pi)^2
            % Ne = ( (fr*2*pi)^2 * m_e * eps0 / e^2 ) - (3 * k^2 * eps0 * k_B * Te / e^2)
            term1 = ( (fr .* 2 .* pi).^2 .* m_e .* eps0 ) ./ (e_charge^2);
            term2 = (3 * k_plasma_line^2 * eps0 * k_B * Te_prior) / (e_charge^2);
            ne_points = term1 - term2;
            
            % Ne Prior (Column 1)
            apriori2(ipl, 1) = mean(ne_points);
            
            % Ne Error (Column 1) in Ne unit; Derivative: dNe/dfr = (8 * pi^2 * m_e * eps0 / e^2) * fr
            df = pl_df(idx_gate);
            ne_err_points = (8 * pi^2 * m_e * eps0 / e_charge^2) .* fr .* df;
            apriorierror2(ipl, 1) = mean(ne_err_points);
            
            apriorierror2(ipl, 6) = 0.5;
        end
    end
end
