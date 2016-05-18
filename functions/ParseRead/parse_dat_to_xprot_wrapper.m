% parse Siemens dat-file to xprotocol and mrprot structures
%  Steven Baete, NYU LMC CBI, Augustus 2013
%  based on the eja_rdMeas.zip-package of E. Auerbach, CMRR, 2013

% in:
%   name: filename
% out:
%   xprot: the xprot-data structured in several structures
%   mrprot: the xprot-data structured in several structures

function [xprot,mrprot] = parse_dat_to_xprot_wrapper(name)

fid = fopen(name);

rawarr = '';
found = false;
while (~found)
    line = fgets(fid);
    rawarr = [rawarr line];
    if strfind(line,'PMULearnPhase')
        found = true;
    end;
end;
ind = strfind(rawarr,'XProtocol');
for i = 1:length(ind)
    mrprot{i} = parse_dat_to_xprot(rawarr((ind(i)-10):end));
end
ind = strfind(rawarr,'### ASCCONV BEGIN');
for i = 1:length(ind)
    mrprot{i} = parse_mrprot(rawarr((ind(i)):end));
end
fclose(fid);