% function [rho,df,rho_sig95] = correlate(xx,yy)
%  computes 
%      the (unlagged) correlation between xx and yy
%      the avg. no. of degrees of freedom in xx,yy
%      the 95% confidence level, rho_sig95
%
% error in D.F. fixed April 2009 (KAK)
%
function [rho,df,rho_sig95] = correlate(xx,yy)
xx=xx(:);yy=yy(:);
mx=length(xx);my=length(yy);
% compute correlation
xp=xx-mean(xx);
yp=yy-mean(yy);
rho=(yp'*xp)/sqrt((xp'*xp)*(yp'*yp));
% compute covariance
covx=real(xcov(xx,xx,'unbiased'));covx=covx/covx(mx); 
covy=real(xcov(yy,yy,'unbiased'));covy=covy/covy(my);
% integrate until zero crossings
% USE TRAPEZOID RULE INSTEAD OF RECTANGLE
taux=0.;
i=mx;
	while covx(i)>=0
	taux=taux+(covx(i)+covx(i+1))/2;
	i=i+1;
	end
tauy=0.;
i=my;
	while covy(i)>=0
	tauy=tauy+(covy(i)+covy(i+1))/2;
	i=i+1;
	end

% find degrees of freedom
df=min(mx,floor((mx/taux+my/tauy)/2));
%
% find 95% confidence correlation
%
b=2*1.96/sqrt(df-3);
rho_sig95=(exp(b)-1)/(1+exp(b));
