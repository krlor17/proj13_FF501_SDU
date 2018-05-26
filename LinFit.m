function [a,b,chisq,ts,tf] = LinFit(times, lnC, npoints, err,tmin,tmax)
	% This function takes the lattice time points $times
	% and the corresponding log(c(t)) values and returns
	% $a : the slope of the best linear fit on npoints
	% $b : the intersection of the best linear fit on npoints
	% $chisq : chi squared value of the best fit
	% $times are the timepoints (x-values) for the fit
	% $npoints is the length of the regression interval
	% $err are the errors for each point determined from resampling
	% $NIntervals is the number of intervals of consecutive $npoints
	% on half of the data points
	NIntervals=(tmax-tmin)-npoints+2;
	% $fits is the matrix that will contain all the fits
	fits=zeros(5,NIntervals);
	for i=1:NIntervals
		ts=tmin+i-1;		%earliest time in current fit interval
		tf=tmin+i+npoints-2;	%latest time in current fit interval
		%% linear fit, uses a Vandermonde matrix and least squares regression
		f = polyfit(times(ts:tf),lnC(ts:tf), 1);
       		fits(1,i) = f(1);	%storing the slope of the linear fit
		fits(2,i) = f(2);	%storing the intersection
		func=f(1)*times(1:64) + f(2)*linspace(1,1,64)';
		%% loop for calculating the chi square value of $func fitted to $lnC
		for k=ts:tf
			fits(3,i) = fits(3,i)+(1/(npoints-2))*((func(k)-lnC(k))^2)/((err(k))^2);
		end	
		fits(4,i)=ts;
		fits(5,i)=tf;
	end
	%we find the smallest chi squared value  and its index
	[val,index] = min(fits(3,:));
	%the function returns the best fit for npoints
	a=fits(1,index);
	b=fits(2,index);
	chisq=fits(3,index);
	ts=fits(4,index);
	tf=fits(5,index);
end
