function hpfilt=hpf(n,l)
%Calculates a matrix that can be used to compute cyclical component.
%     inputs:   n = number of observations in data set
%               l = lambda (see Prescott, "Theory Ahead of Measurement..")
%     output:  hpfilt = n x n matrix    
nobs=n;
a=zeros(nobs,nobs);
a(1,1)=1+l;
a(1,2)=-2*l; a(1,3)=l;
a(2,1)=-2*l; a(2,2)=1+5*l;
a(2,3)=-4*l; a(2,4)=l;
a(nobs-1,nobs-3)=l;
a(nobs-1,nobs-2)=-4*l;
a(nobs-1,nobs-1)=1+5*l;
a(nobs-1,nobs)=-2*l;
a(nobs,nobs-2)=l;
a(nobs,nobs-1)=-2*l; 
a(nobs,nobs)=1+l;
for i=3:nobs-2
     a(i,i)=1+6*l;
     a(i,i+1)=-4*l;
     a(i,i+2)=l;
     a(i,i-1)=-4*l;
     a(i,i-2)=l;
end
hpfilt=eye(nobs)-inv(a);
