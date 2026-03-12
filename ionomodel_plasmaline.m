function [altitude,ne,te,ti,coll,cO,cM2,cH]=ionomodel_plasmaline(heights,modinfo)
%
% Model ionosphere for GUISDAP analysis with plasma line Ne. This function just calls ionomodel_iri
% to get the default values. The plasma line data are included in apriorimodel_plasmaline.m
% 
%
% Ilkka Virtanen, 2026- 
%

[altitude,ne,te,ti,coll,cO,cM2,cH] = ionomodel_iri(heights,modinfo);



end

