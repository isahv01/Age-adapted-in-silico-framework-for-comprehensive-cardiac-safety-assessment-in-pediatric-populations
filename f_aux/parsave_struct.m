% this function fixed the problem of not being able to use 'save' inside a
% parfor loop. 
function parsave_struct(filename, S)
    save(filename, '-struct', 'S');
end

