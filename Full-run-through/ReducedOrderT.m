function U= ReducedOrderT(M, K, C, F, nModes, nFreq, t)
    [Mode, ~] = eigs(K, M, nModes, 'smallestabs');                          % Calculate first nModes mode shapes
    M_phi= Mode'*M*Mode;
    K_phi= Mode'*K*Mode;
    C_phi= Mode'*C*Mode;
    m= diag(M_phi); k= diag(K_phi); c= diag(C_phi);
    
    % Calculate T for every time
    T= Mode'* F;    
    
    %Frequency sample length and stepping
    Fs = 1/(t(end)-t(1)); 
    L = length(t)-rem(length(t),2);                                         % L needs to be even 
    omega = 2*pi*Fs*(-L/2:L/2);                                             % Include negative frequencies
    
    U= zeros(length(M), length(t));    
    for I= 1:nModes
        %Fourier transform of T for each mode and shifted
        fftT= fft(T(I,:));
        xi= fftshift(fftT);                                                 % Shift FFT to incorperate peaks at negative frequencies
        
        %Find peak frequency indices
        [~,locs]= findpeaks(abs(xi),'Npeaks',nFreq,'SortStr','descend');

        for i= locs                                                        % 1:number of peaks for this mode
            %Calculate X for that peak omega and mode
            X=(xi(i)/(k(I)+1j*c(I)*omega(i)-m(I)*omega(i)^2));              % Damping must be included   
            %Calculate U
            U= U+(Mode(:,I)*X*exp(1j*omega(i)*t)/L);                        %Divide by L to act as a ifft
        end 
    end
end