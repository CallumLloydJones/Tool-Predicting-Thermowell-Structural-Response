function centre = elCentre(nodes,elements)
    centre = zeros(length(elements),3);
    for e = 1:length(elements)
        a = elements(:,e);
        npos = [nodes(a(1),:); nodes(a(2),:); nodes(a(3),:); nodes(a(4),:);];
        centre(e,:) = mean(npos);
    end
end
