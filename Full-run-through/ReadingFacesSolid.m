function faces= ReadingFacesSolid(path)
    %%Reading the boundaries starting face and number of faces
    filename = append(path,'constant\polyMesh\boundary');
    filetext = fileread(filename);
    %Read from the word Thermowell to fixed bottom
    expr  = 'Thermowell';
    expr2 = 'fixed_Bottom';
        
    begining = regexp(filetext,expr);
    ending = regexp(filetext,expr2);
    text = filetext(begining:ending);
    text = regexp(text,'\d*','Match'); % Obtain numbers

    StartFace = str2double(text(2)); % Assign numbers to correct label
    nFaces = str2double(text(1));

    %%Using the boundaries starting face and number of faces to get faces
    opt = detectImportOptions(append(path,'constant\polyMesh\faces'));
    opt = setvartype(opt, [1:3], 'double');
    opt = setvaropts(opt, 'Prefixes', {'3('}, 'Suffixes', {')'});
    faces = readmatrix(append(path,'constant\polyMesh\faces'),opt);
    %Read faces from start face to end face
    faces = faces(StartFace+1:StartFace+nFaces,:); 
end