function setGlobalGPUParam(osf_in, wg_in, sw_in)

global osf % oversampling: 1.5 1.25
global wg % kernel width: 5 7
global sw % parallel sectors' width: 12 16

osf = osf_in; 
wg = wg_in; 
sw = sw_in; 
end

function [osf_in, wg_in, sw_in] = getGlobalGPUParam

global osf % oversampling: 1.5 1.25
global wg % kernel width: 5 7
global sw % parallel sectors' width: 12 16

osf_in = osf; 
wg_in = wg; 
sw_in = sw;

end



