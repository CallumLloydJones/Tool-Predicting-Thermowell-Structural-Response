function [fAreas,fCentres] = ReadingCentres(path,faces)
    %%Reading the boundaries starting face and number of faces
    filename = append(path,'constant\polyMesh\points');

    %%Using the boundaries starting face and number of faces to get faces
    opt = detectImportOptions(filename);
    opt = setvartype(opt, [1:3], 'double');
    opt = setvaropts(opt, 'Prefixes', {'('}, 'Suffixes', {')'});

    points = readmatrix(filename,opt);
    
    fCentres=zeros(length(faces),3); % Initialise F vector
    for i= 1:length(faces)
        %Calculate Area of faces of which traction is applied
        verts= points(faces(i,:)+1,:);
        fCentres(i,:) = mean(verts);
    end
       
    fAreas = zeros(length(faces),1); % Initialise F vector
    for i=1:length(faces)
        P1= points(faces(i,1)+1,:);
        P2= points(faces(i,2)+1,:);
        P3= points(faces(i,3)+1,:);
        fAreas(i)=0.5*norm(cross(P1-P2,P1-P3));
    end
    
end