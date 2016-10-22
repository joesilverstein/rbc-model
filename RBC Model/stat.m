function [stdpct,corrcon,corrlag]=stat(sim,hp,nlags)
%Returns summary statistics for a simulation.
%Inputs: SIM-columns correspond to variables and rows to observations.
%        HP-Hodrick-Prescott filter.
%        NLAGS-Number of leads and lags to use when computing CORRLAG
%               (generally 5).
%Outputs: STDPCT-Percent standard deviation.
%         CORRCON-Contemporaneous corrleation matrix.
%         CORRLAG-Correlations of NLAGS leads and lags of each variable with the 
%                  FIRST column of SIM (generally output).

%The calling program should invoke this function for each simulation and
%store the results.  The function HPF.M can be used to create the matrix HP.
%The function STATDISP.M can be used to compute averges
%and standard deviations of these statistics and display the results.

%Number of observations is the number of rows in SIM.
%Each column of SIM corresponds to a variable.

[nobs,nvars]=size(sim);
simd=hp*sim;
stdpct=std(simd,1)*100; %Computes percent standard deviations.

corrcon=cov(simd,1); %Computes contemporaneous correlation matrix.
vsimd = diag(corrcon);
corrcon = corrcon./sqrt(vsimd*vsimd');

%Compute lead-lag correlations with the first column of SIM.
for l=-nlags:0
    corlag=cov([simd(1-l:nobs,1) simd(1:nobs+l,:)],1);
    vsimd=diag(corlag);
	corlag=corlag./sqrt(vsimd*vsimd');
	corrlag(l+nlags+1,:)=corlag(1,2:nvars+1);
    corlag=cov([simd(1:nobs+l,1) simd(1-l:nobs,:)],1);
	vsimd=diag(corlag);
	corlag=corlag./sqrt(vsimd*vsimd');
	corrlag(-l+nlags+1,:)=corlag(1,2:nvars+1);
end
