function Forces= ReadingForces(path, fFaces, time)
    %%Reading the boundaries starting face and number of faces
    filename = path+string(time)+'\force';

    %%Using the boundaries starting face and number of faces to get faces
    opt = detectImportOptions(filename);
    opt = setvartype(opt, [1:3], 'double');
    opt = setvaropts(opt, 'Prefixes', {'('}, 'Suffixes', {')'});
    Forces = readmatrix(filename,opt);
    Forces(:, 4:end) = []; 
    Forces(length(fFaces)+1:end, :) = [];
end