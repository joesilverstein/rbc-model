function statdisp(stdpct,corrcon,corrlag)
disp('Standard Deviations')
disp(mean(stdpct,3))
disp('Contemporaneous Correlations')
disp(mean(corrcon,3))
x=size(corrlag);
lag=x(1);
lag=(lag-1)/2;
disp('Correlation(Y(t),Variable(t+j))')
disp([[-lag:lag]' mean(corrlag,3)])
disp('Sample Standard Deviations of the Above Statistics:')
disp(' ')
disp('Standard Deviations')
disp(std(stdpct,1,3))
disp('Contemporaneous Correlations')
disp(std(corrcon,1,3))
disp('Correlation(Y(t),Variable(t+j))')
disp([[-lag:lag]' std(corrlag,1,3)])
