function [err lnerr replicas lnreplicas] = JKR(DATA,tpoints)
%this function takes a dataset DATA of 2 columns
%column 1 is the lattice time points
%column 2 is average c(t) values
% $tpoints is the number of lattice times in each run
% Returns:
% $err the std. deviation of the non-log Jackknife replicas
% $lnerr the std. deviation of the log Jackknife replicas
% $replicas the non-log Jackknife replicas
% $lnreplicas the log Jackknife replicas
	%% we consider only half the data, as we are not interested
	%% in the parity partner
	t=round(tpoints/2);			 
	N=round(length(DATA(:,2))/tpoints);	%Number of replicas
	replicas=zeros(t,N);
	lnreplicas=zeros(t,N);
	
	%making the non-log replicas
	C=ReStruct(DATA(:,2),tpoints,1);
	MEAN=mean(C')';
	for i=1:N
		replicas(:,i)=(N*MEAN-C(:,i))/(N-1);
	end
	err=zeros(t,1);
	SUM=err;
	for i=1:N
		SUM=SUM+(MEAN-replicas(:,i)).^2;
	end

	for i=1:t
		err(i)=sqrt(SUM(i))*sqrt((N-1)/N);
	end

	%making the log replicas
	C=ReStruct(DATA(:,2),tpoints,0);
	MEAN=mean(C')';
	for i=1:N
		lnreplicas(:,i)=(N*MEAN-C(:,i))/(N-1);
	end
	lnerr=zeros(t,1);
	SUM=lnerr;
	for i=1:N
		SUM=SUM+(MEAN-lnreplicas(:,i)).^2;
	end

	for i=1:t
		lnerr(i)=sqrt(SUM(i))*sqrt((N-1)/N);
	end
end
