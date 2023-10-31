function [eqStress, p] = stress(nodes, elements, U, E, nu, t)
    G = (1+nu)/E;
    D= [1/E -nu/E -nu/E 0 0 0;
        -nu/E 1/E -nu/E 0 0 0;
        -nu/E -nu/E 1/E 0 0 0;
        0 0 0 G 0 0;
        0 0 0 0 G 0;
        0 0 0 0 0 G];
    
    derNxi = [1 0 0 -1;0 1 0 -1; 0 0 1 -1;];
    UU= reshape(U, [length(nodes), 3, length(t)]);
    eqStress = zeros(length(elements),length(t));
    p = zeros(length(elements),length(t));
    for ti=1:length(t)
        for e = 1:length(elements)
            a = elements(:,e);
            npos = [nodes(a(1),:); nodes(a(2),:); nodes(a(3),:); nodes(a(4),:);]';
            u = [UU(a(1),:,ti); UU(a(2),:,ti); UU(a(3),:,ti); UU(a(4),:,ti);];
    
            f = npos * derNxi';
            derNx = inv(f)'*derNxi;
            
            strain = 0.5 * ((derNx*u) + (derNx*u)');
            strainV = [strain(1,1) strain(2,2) strain(3,3) 2*strain(2,3) 2*strain(1,3) 2*strain(1,2)]';
            sigmaV = D\strainV;
            
            eqStress(e,ti) = sqrt(1.5*((sigmaV(1)-sigmaV(2))^2 + (sigmaV(2)-sigmaV(3))^2 + (sigmaV(3)-sigmaV(1))^2));
            p(e,ti) = (1/3)*sum(sigmaV(1:3));
        end
    end
end
