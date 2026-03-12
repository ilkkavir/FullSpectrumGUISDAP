function [apriori2,apriorierror2] = apriorimodel_plasmaline(apriori,apriorierror,heights,fit_altitude)
%
%
% [apriori2,apriorierror2] = apriorimodel_plasmaline(apriori,apriorierror,heights,fit_altitude)
%
% Ne prior from plasma line data and O+ ion fraction fit when plasma lines are available. The default GUISDAP prior otherwise. 
%
% set iono_model='plasmaline' in GUISDAP to activate this function. 
%
% INPUT:
%        apriori       initial prior parameters as returned ionomodel_bafim
%        apriorierror  prior errors calculated in ionoomdel. The matrix is not used in this version of bafim.
%        heights       centre points of height gates
%        fit_altitude  the fit_altitude matrix used internally in GUISDAP
%
% OUTPUT:
%        apriori2      modified prior values with Ne from plasma lines when available, and O+ ion fraction fit enabled in these cases. 
%        apriorierror2 standard deviations of the prior values. 
% 
%
%
% Ilkka Virtanen & Tikemani Bag, 2026-
%
    
    % global variables used in BAFIM. We will probably need only d_time in this case. 
    global r_param r_error r_res r_status d_time path_GUP result_path v_Boltzmann

    % copy the original prior model
    apriori2 = apriori;
    apriorierror2 = apriorierror;

    % The variables
    %
    % heights        an array of height gate center point altitudes (km)
    % apriori2       a length(heights) x N array of prior values. The first six columns are for Ne, Ti, Tr=Te/Ti, Coll, Vi, and [O+].
    % apriorierror2  an array of standard deviations for the prior values.
    % fit_altitude   included for compatibility with other prior models. Not needed here.
    % d_time         start and end time of the integration period. 


    % Then replace the Ne prior values and relax the O+ fraction prior standard deviation when plasma line data are available as explained below
    % 
    % 1. find all altitudes from which we have plasma line data during the integration period (from d_time(1,:) to d_time(2,:))
    % 2. At each gate (ipl):  If there is plasma line data from gate ipl, replace values in apriori2 and apriorierror2 as follows
    %
    %    - apriori2(ipl,1) = <Ne from plasma lines>
    %    - apriorierror2(ipl,1) = <width of the plasma line>  (in units of Ne)
    %    - apriorierror2(ipl,6) = 0.5

    % Unfortunately, edges of the height gates are not easily available. You will need to find  approximate gate edges from the center points in the heights array.
    %
    % Averaging will be needed if there are many plasma lines from one height gate. This may happen, because finer altitude resolution is used for the plasma line analysis.
    %
    % We can first assume that the plasma line frequency shift is equal to the local plasma frequency. 
    % To include the Debye legth correction to the plasma line frequency (See equation 1 of Fredriksen et al., 1989), you can use the IRI Te, which can be calculated as apriori(:,2).*apriori(:,3).
    %
    % Center frequencies of the plasma line channels for different EISCAT experiments and their versions
    %
    % Tromso UHF radar
    % beata 1.0  -5.0 MHz
    % beata 2.0  -3.6 MHz, -6.0 MHz, -8.4 MHz,
    % beata 2.2  -2.5 MHz, -4.9 MHz
    % beata 2.3  -3.85 MHz
    %
    % There are more, including many ESR experiments in particular, but will need to investigate a bit to find the frequencies. 
    % 
    % 
    % 
    %
    
    
end

