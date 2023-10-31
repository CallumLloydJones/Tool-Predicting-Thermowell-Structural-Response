function St = plotliftdrag(rootdir, mesh, start, U, D, label)
    data = dlmread(append(rootdir,'postProcessing\forceCoeffs\0\coefficient.dat'), "\t", 13, 0);
    t = data(:,1);
    Cd = data(:,2);
    Cl = data(:,4);
    
    figure;
    plot(t(t>start), Cd(t>start));
    xlabel("time (s)");
    ylabel("Drag Coefficent");
    grid;
    legend([mesh]);
    saveas(gcf,append(rootdir+'\CD.png'))
    
    Cd= Cd(t > start );
    Cl= Cl(t > start );
    t = t(t > start);

    figure;
    plot(t(t>start), Cl(t>start));
    xlabel("time (s)");
    ylabel("Lift Coefficent");
    grid;
    saveas(gcf,append(rootdir+'\CL.png'))

    figure;
    N = length(Cl);
    fftCl = fft(Cl);
    dt = t(2) - t(1);
    f = linspace(0, 1 / dt, N);
    plot(f(1:floor(N/2)), abs(fftCl(1:floor(N/2))) / N);
    xlabel("Frequency (Hz)");
    ylabel("Lift Coefficent");
    xlim([0,200])
    legend([label]);
    grid;
    saveas(gcf,append(rootdir+'\CLfft.png'))
    
    
end
