%%% Data analysis FF501%%%
function [] = DataJK(datafile,name)
% $datafile is the name of the file the data is imported from
% $name is the name of the particle e.g. 'pion'

%%% Initialization of output files %%%
%Latex table containing data
fileID = fopen(['results/' name '/' name '_table.tex'], 'w');
fprintf(fileID, '\\begin{table} \n \\centering \n \\caption{fits %s} \n \\label{tab:fits_%s} \n',name,name);
fprintf(fileID, '\\begin{tabular}{|c| >{\\centering}m{100 pt} | c| c|}\\hline \n');
fprintf(fileID, '\\# points & mass $ \\si{ \\frac{MeV }{ c^2} } $ & $ \\chi^2 $ & time interval \\\\ \\hline \n');
fclose(fileID);

%Latx table containing the final quotable reults
fileID = fopen(['results/' name '/' name '_final_table.tex'], 'w');
fprintf(fileID, '\\begin{table} \n \\centering \n \\caption{%s} \n \\label{tab:%s} \n',name,name);
fprintf(fileID, '\\begin{tabular}{|c|c|c|c|c|}\\hline \n');
fprintf(fileID, 'avg. mass $\\hat{m} \\si{ \\frac{MeV }{ c^2} } $ & error avg. mass  $\\si{ \\frac{MeV }{ c^2} }& median mass $\\hat{m} \\si{ \\frac{MeV }{ c^2} } $ & error median $\\si{ \\frac{MeV }{ c^2} } & Systematic error $ \\si{ \\frac{MeV }{ c^2} } $ \\\\ \\hline \n');
fclose(fileID);

%Latex table for effective mass estimate
fileID = fopen(['results/' name '/' name '_effmass_table.tex'], 'w');
fprintf(fileID, '\\begin{table} \n \\centering \n \\caption{} \n \\label{tab:%s effm} \n',name);
fprintf(fileID, '\\begin{tabular}{|c| >{\\centering}m{100 pt} | c|}\\hline \n');
fprintf(fileID, ' effecive mass $ \\si{ \\frac{MeV }{ c^2} } $ & error $ \\si{ \\frac{MeV }{ c^2} } $  & time interval \\\\ \\hline \n');
fclose(fileID);

%Data file containing results in lattice units
fileID=fopen(['results/' name '/' name '_results_jackknife.dat'],'w');
fprintf(fileID,' %-12s  \t %-12s \t %-12s \t %-12s \t %-12s \t %-12s \t %-12s \r\n','#points     ','mass        ','error mass','coeff       ','chi sq     ','ts','tf');
fclose(fileID);

%Data file containing results in physical units (OBS. not coefficient A)
fileID=fopen(['results/' name '/' name '_results_jackknife_units.dat'],'w');
fprintf(fileID,'%s \r\n', 'Unit of the mass is MeV/c^2');
fprintf(fileID,'%-12s \t %-12s \t %-12s \t %-12s \t %-12s \t %-12s \t %-12s \r\n','#points     ','mass        ','error mass','coeff       ','chi sq      ','ts','tf');
fclose(fileID);

%%% Parameters %%%
tpoints=128;		%Number of points in time interval
startn=5;		%Lowest  # points in the regression intervals
tmin=10;		%Earliest time to consider
tmax=40;		%Latest time to consider in fit

%%% DATA %%%
Imp = importdata(datafile,' ',1);
times= Imp.data(1:tpoints,1);		 %the first time series
c=ReStruct(Imp.data(:,2),tpoints,0);	 %matrix for log(c(t)) values
cj=ReStruct(Imp.data(:,2),tpoints,1)';	 %matrix for c(t) values
MEAN=mean(cj)';
lnMEAN=mean(c')';
%Jackknife resampling
[errM lnerrM replicas lnreplicas]=JKR(Imp.data(:,:),tpoints);

masses=zeros(1,(tmax-tmin-startn+2));
errors=masses;
chis=masses;

%%% Evaluaton loop %%%
for npoints=startn:(tmax-tmin+1)

%% DATA Analysis %%
[b A chisq ts tf]=LinFit(times,lnMEAN, npoints,lnerrM,tmin,tmax);	%Lin. fit of ln(c(t))
% estimating the error
N=length(replicas);
SUMm=0;
jackm=zeros(N,1);
for i=1:N
[jackm(i) garb1 garb2 garb3 garb4]=LinFit(times, lnreplicas(:,i),npoints,lnerrM,tmin,tmax);
SUMm=SUMm+(jackm(i)-b)^2;
end
err=sqrt((N-1)/N)*sqrt(SUMm);

%% Results %%
m=abs(b);		%the mass is the negative slope
hca=5628.7;		%conversion factor to physical units

masses(1,(npoints-startn +1))=m;
errors(1,(npoints-startn +1))=err;
chis(1,(npoints-startn +1))=chisq;

%file containing results in lattice units
fileID=fopen(['results/' name '/' name '_results_jackknife.dat'],'a');
fprintf(fileID,'%12d \t %12d \t %12d \t %12d \t %12d \t %12d \t %12d \r\n',npoints,m,err, A, chisq,ts,tf);
fclose(fileID);

%file containing results in physical units
fileID=fopen(['results/' name '/' name '_results_jackknife_units.dat'],'a');
fprintf(fileID,'%12d \t %12d \t %12d \t %12d \t %12d \t %12d \t %12d \r\n',npoints,m*hca,err*hca,A, chisq, ts, tf);
fclose(fileID);

%latex table of results in physical units
fileID = fopen(['results/' name '/' name '_table.tex'], 'a');
fprintf(fileID, '$ %8.2i $ & $ %8.2f  \\pm %8.2f $ ', npoints, vpa(m*hca,5),vpa(err*hca,4));
fprintf(fileID, ' & $  %s $ & $ %i - %i $ \\\\ \\hline',  latex(vpa(chisq,3)),ts,tf);
fprintf(fileID, '\n');
fclose(fileID);

end

%%% End of evaluation loop %%%

% ending latex table
fileID = fopen(['results/' name '/' name '_table.tex'], 'a');
fprintf(fileID, '\\end{tabular}\n');
fprintf(fileID, '\\end{table}\n');
fclose(fileID);


%%% Systematic error and final results %%%

%% using the average mass with jackknife for error propagation
mhat=mean(masses);
Nm=length(masses);
mrep=zeros(1,Nm);
for i=1:Nm
	mrep(i)=(Nm*mhat-masses(i))/(Nm-1);
end
summ=0;
for i=1:Nm
	summ=summ+(mhat-mrep(i))^2
end
errhat=sqrt(summ);

%%using the 'median' of the sorted masses
sortmasses=sort(masses);
mederr=0;
medm=masses(round(length(masses))); 	%central value
for i=1:length(masses)
	if medm == masses(i)
		mederr=errors(i);
	end
end

%%systematic error from Range
[mmax,indexmax]=max(masses);
[mmin,indexmin]=min(masses);
syserr=(mmax + errors(indexmax) - mmin + errors(indexmin))/2;

% Ouput final table
fileID = fopen(['results/' name '/' name '_final_table.tex'], 'a');
fprintf(fileID, '$ %8.2f $ & $\\pm  %8.2f$ & $ %8.2f $ & $\\pm  %8.2f$ &  $\\pm %8.2f$ ', vpa(mhat*hca,5), vpa(errhat*hca,5), vpa(medm*hca,5), vpa(mederr*hca,5), vpa(syserr*hca,5));
fprintf(fileID, '\n');
fprintf(fileID, '\\end{tabular}\n');
fprintf(fileID, '\\end{table}\n');
fclose(fileID);


%%% Effective mass %%%
N=length(replicas);
t=round(tpoints/2);
effM=zeros(1,t-1);
for i=1:(t-1)
	effM(i)=log(MEAN(i)/MEAN(i+1));
end

effMj=zeros(t-1,N);

for j=1:N
	for i=1:(t-1)
		effMj(i,j)=log(replicas(i,j)/replicas(i+1,j));
	end	
end
%%Error for effective mass
errMj=zeros(t-1,1);
SUM=errMj;
MJ=mean(effMj')';
for i=1:N
	SUM=SUM+(MJ-effMj(:,i)).^2;
end

for i=1:t-1
	errMj(i)=sqrt(SUM(i))*sqrt((N-1)/N);
end

%output effective mass table
fileID = fopen(['results/' name '/' name '_effmass_table.tex'], 'a');
for i=1:t-1
	fprintf(fileID, '$ %8.2f $ & $ %8.2f$  & $  %i-%i $ \\\\ \\hline \n ', effM(i)*hca, errMj(i)*hca, times(i),times(i+1));
end
fprintf(fileID, '\\end{tabular}\n');
fprintf(fileID, '\\end{table}\n');
fclose(fileID);



%%% Plots %%%

%custom color
lightblue=1/255*[120,120,120];

%% plot of mean data with error
plot1=figure(1)
h=errorbar(times(1:64),MEAN,errM);
h.Color='k'
h.CapSize=3;
h.Marker='o';
h.MarkerSize=2;
axis([0 64 (min(MEAN)-2*min(errM(1:64))) (max(MEAN)+2*max(errM(1:64)))])
set(get(h,'Parent'),'Yscale','log')
xlabel('Lattice time','FontName','MathJax_typewriter','FontSize',14)
ylabel('Correlation function C(t)','FontName','MathJax_typewriter','FontSize',14)
plot1.PaperUnits='centimeters';
plot1.PaperPosition=[0,0,12,10];
print(plot1,['results/' name '/' name '_mean_err'],'-dpng')

%% plot of effective mass curve
plot2=figure(2)
h=errorbar((times(1:63)+1),effM(1:63)*hca,errMj(1:63)*hca,'Color','k');
h.CapSize=3;
h.Marker='o';
h.MarkerSize=2;
xlabel('Lattice time','FontName','MathJax_typewriter','FontSize',14)
ylabel('Effective mass [MeV/c^2]','FontName','MathJax_typewriter','FontSize',14)
xlim([0 64])
plot2.PaperUnits='centimeters';
plot2.PaperPosition=[0,0,12,10];
print(plot2,['results/' name '/' name '_effmass'],'-dpng')

%% plot of effective mass with avg. estimate
plot3=figure(3)
hold on
h=errorbar((times(1:63)+0.5),effM(1:63)*hca,errMj(1:63)*hca,'Color','k');
h.CapSize=3;
h.Marker='o';
h.MarkerSize=2;
xlabel('Lattice time','FontName','MathJax_typewriter','FontSize',14)
ylabel('Effective mass [MeV/c^2]','FontName','MathJax_typewriter','FontSize',14)
xlim([0 64])
l1=line([0 64],[mhat*hca mhat*hca],'Color','blue','LineStyle','--')
l2=line([0 64],[(mhat-errhat-syserr)*hca (mhat-errhat-syserr)*hca],'Color',lightblue,'LineStyle','--')
line([0 64],[(mhat+errhat+syserr)*hca (mhat+errhat+syserr)*hca],'Color',lightblue,'LineStyle','--')
hold off
legend([h l1 l2],'Effective mass','avg. estimate','error')
plot3.PaperUnits='centimeters';
plot3.PaperPosition=[0,0,12,10];
print(plot3,['results/' name '/' name '_effmass_estimateline'],'-dpng')

%% plot of effective mass with median estimate
plot4=figure(4)
hold on
h=errorbar((times(1:63)+0.5),effM(1:63)*hca,errMj(1:63)*hca,'Color','k');
h.CapSize=3;
h.Marker='o';
h.MarkerSize=2;
xlabel('Lattice time','FontName','MathJax_typewriter','FontSize',14)
ylabel('Effective mass [MeV/c^2]','FontName','MathJax_typewriter','FontSize',14)
xlim([0 64])
l1=line([0 64],[medm*hca medm*hca],'Color','blue','LineStyle','--')
l2=line([0 64],[(medm-mederr-syserr)*hca (medm-mederr-syserr)*hca],'Color',lightblue,'LineStyle','--')
line([0 64],[(medm+mederr+syserr)*hca (medm+mederr+syserr)*hca],'Color',lightblue,'LineStyle','--')
hold off
legend([h l1 l2],'Effective mass','median estimate','error')
plot4.PaperUnits='centimeters';
plot4.PaperPosition=[0,0,12,10];
print(plot4,['results/' name '/' name '_effmass_median'],'-dpng')

%% plot of Signal / Noise ratio
% Signal is the mean correlator data
% Noise is the error estimated by jackknife resampling
plot5=figure(5)
mval=MEAN;
tval=linspace(0,63,64);
R=mval./errM;
e=plot(tval,R')
e.Color='black';
xlabel('Lattice time','FontName','MathJax_typewriter','FontSize',14)
ylabel('Signal / Noise','FontName','MathJax_typewriter','FontSize',14)
xlim([0 64])
plot5.PaperUnits='centimeters';
plot5.PaperPosition=[0,0,12,10];
print(plot5,['results/' name '/' name '_SNratio'],'-dpng')

%end of function
end
