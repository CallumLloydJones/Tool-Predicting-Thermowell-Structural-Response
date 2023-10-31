function F= ApplyF(force,nodes,faces)
    F=zeros(length(nodes),3); % Initialise F vector
    for T= 1:length(force)
        %Calculate Area of faces of which traction is applied
        P1= nodes(faces(T,1)+1,:);
        P2= nodes(faces(T,2)+1,:);
        P3= nodes(faces(T,3)+1,:);
        A=0.5*norm(cross(P1-P2,P1-P3));
        for i=faces(T,:) % Apply tractions to the nodes of each face
            F(i+1,:)= F(i+1,:) + force(T,:)/3;
        end
    end
    F= reshape(F, [length(nodes)*3,1]); %Reshape vector to 1 column
end