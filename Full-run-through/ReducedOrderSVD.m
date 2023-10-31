function U= ReducedOrderSVD(M, K, C, F, nModes, nFreq, rank, t)
    [Mode, ~] = eigs(K, M, nModes, 'smallestabs');                          % Calculate first nModes mode shapes
    M_phi= Mode'*M*Mode;
    K_phi= Mode'*K*Mode;
    C_phi= Mode'*C*Mode;
    m= diag(M_phi); k= diag(K_phi); c= diag(C_phi);

    %Steady state for is removed
    steady = mean(mean(F,2));
    F = F-steady;
    
    [U1,S1,V1] = svd(F);

    %Frequency sample length and stepping
    Fs = 1/(t(end)-t(1)); 
    L = length(t)-rem(length(t),2);                                          % L needs to be even 
    omega = 2*pi*Fs*(-L/2:L/2);                                               % Include negative frequnecies

    U= zeros(length(M), length(t));
    %To reconstruct F = U1S1V1' + U2S2V2' + U3S3V3' + ... + UmSmVm'
    for a = 1:rank
        %Lets truncate and take the most nFreq important USV 
        shape= U1(:,a);
        load= S1(a,a)*V1(:,a)';
        
        % Calculate T for every time
        F_nom= Mode'* shape; 

        % Fourier transform of the load and shifted
        fftL = fft(load); 
        xi= fftshift(fftL);
        [~,locs]= findpeaks(abs(xi),'Npeaks',nFreq,'SortStr','descend');
        
        for I= 1:nModes
            for i= locs % 1:number of peaks for this mode
                %Calculate X for that peak omega and mode
                X=(F_nom(I)*xi(i)/(k(I)+1j*c(I)*omega(i)-m(I)*omega(i)^2));      % Damping must be included   
                %Calculate U
                U=U + (Mode(:,I)*X*exp(1j*omega(i)*t)/length(t));               % Divide by L to act as a ifft
            end 
        end
    end
    
    %Steady state force will also be accounted for
    steady= steady*ones(size(F));
    [U2,S2,V2] = svd(steady);

    shape= U2(:,1);
    load= S2(1,1)*V2(:,1)';
    
    % Calculate T for every time
    F_nom= Mode'* shape; 

    % Fourier transform of the load and shifted
    fftL = fft(load); 
    xi= fftshift(fftL);
    [~,locs]= findpeaks(abs(xi),'Npeaks',1,'SortStr','descend');

    for I= 1:nModes
        %Calculate X for that peak omega and mode
        X=F_nom(I)*xi(locs)/k(I);      % Damping must be included   
        %Calculate U
        U=U + (Mode(:,I)*X*exp(1j*0*t)/length(t));               % Divide by L to act as a ifft
    end    
    
end