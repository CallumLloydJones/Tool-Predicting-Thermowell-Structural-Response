function faces= ReadingFacesFluid(path)
    %%Reading the boundaries starting face and number of faces
    filename = append(path,'constant\polyMesh\boundary');
    filetext = fileread(filename);
    %Read from the word Thermowell to fixed bottom
    expr  = 'Thermowell';
    expr2 = 'Inlet';
        
    begining = regexp(filetext,expr);
    ending = regexp(filetext,expr2);
    text = filetext(begining:end);
    text = regexp(text,'\d*','Match'); % Obtain numbers

    StartFace = str2double(text(3)); % Assign numbers to correct label
    nFaces = str2double(text(2));

    %%Using the boundaries starting face and number of faces to get faces
    opt = detectImportOptions(append(path,'constant\polyMesh\faces'));
    opt = setvartype(opt, [2], 'char');
    opt = setvaropts(opt, 'Prefixes', {'('}, 'Suffixes', {')'});
    
    faces = readtable(append(path,'constant\polyMesh\faces'),opt);
    %Read faces from start face to end face
    faces = faces(StartFace+1:StartFace+nFaces,:);
    faces = faces.Var2;
    faces = cellfun(@(x) strsplit(x, ' '), faces, 'UniformOutput', false);
    faces = vertcat(faces{:});
    faces = str2double(faces);
end