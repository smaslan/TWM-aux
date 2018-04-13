clear all;

mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);

try
    load('sine_ncorr.matsc','-v6','res','vr','p','s');
end

names = fieldnames(res{1});

rr = {};
for k = 1:numel(res)    
    rr{end+1} = struct();
    for n = 1:numel(names)
        data = getfield(res{k},names{n});
        if ~isscalar(data)
            res{k} = setfield(res{k},names{n},single(data));
        else
            rr{k} = setfield(rr{k},names{n},data);
        end
    end
    k    
end

save('sine_ncorr_min.matsc','-v7','rr','vr','p','s');

save('sine_ncorr_sgl.matsc','-v6','res','vr','p','s');






