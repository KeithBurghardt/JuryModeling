(* ::Package:: *)

function [y, t] = Influencecdf_paramspace(file,dt,Tboundaries)
% convert simulations into a continuous distribution
% Version 1.0 (2015)
% Copyright (C) 2015 Keith Burghardt (University of Maryland, College Park)
% Distributed under GNU GPL v3 .0
% http://www.gnu.org/copyleft/gpl.html
% comes with ABSOLUTELY NO WARRANTY
% ---------------------------------------------------------------
% ---------------------------------------------------------------
% ---------------------------------------------------------------
directory = '/Users/networklab/Desktop/SimData/';
file = strcat(directory,file);
Tboundaries = reshape(Tboundaries,numel(Tboundaries),1);
t = Tboundaries;
[BinomialRange,MVMRange,TauRange,TauHungRange,AlphaRange,AlphaHungRange,N,NumVals] = GetRanges(file)
%{
% string positions where "number of jurors" will be listed
pos1=strfind(file,'N=');
pos2=strfind(file,';Binomial_p='); 
pos3=strfind(file,';MVM_p=');
pos4=strfind(file,';Alpha=');

pos5=strfind(file,';Alpha_hung=');
pos6=strfind(file,';Tau=');
pos7=strfind(file,';Tau_hung=');
pos8=strfind(file,';NumRuns=');
pos9=strfind(file,';NumTrials=');

N = str2num(file(pos1+2:pos2-1));
disp('N = ');
disp(N);
NumRuns = str2num(file(pos8+9:pos9-1));
NumTrials = str2num(file(pos9+11:pos9+12));
NumVals = NumRuns * NumTrials;

% below order is the order in which the data is written
BinomialRange=str2num(strrep(file(pos2+12:pos3-1),'_',':'));
MVMRange=str2num(strrep(file(pos3+7:pos4-1),'_',':'));
TauRange=str2num(strrep(file(pos6+5:pos7-1),'_',':'));
TauHungRange=str2num(strrep(file(pos7+10:pos8-1),'_',':'));
AlphaRange=str2num(strrep(file(pos4+7:pos5-1),'_',':'));
AlphaHungRange=str2num(strrep(file(pos5+12:pos6-1),'_',':'));
%}
fprintf('Parsing File Data...\n');
fid = fopen(file)
for i=1:30
    tline = fgetl(fid);
end

y = zeros(numel(BinomialRange),numel(MVMRange),numel(TauRange),numel(TauHungRange),numel(AlphaRange),numel(AlphaHungRange),N+1,numel(Tboundaries));
size(y)

% NumGuiltyVotes = (0:N).*NumVotes;
% NumInnocentVotes = (N:-1:0).*NumVotes;
count=0;
for b = 1:numel(BinomialRange)
    for mvm = 1:numel(MVMRange)
        for tau = 1:numel(TauRange)
            for tau_hung = 1:numel(TauHungRange)
                for a = 1:numel(AlphaRange)
                    for a_hung = 1:numel(AlphaHungRange)
                        % isolate times by votes
			data = textscan(fid,'% f',2*NumVals);% fscanf(fid,'% f');
			x = [data{1}(1:2:end) data{1}(2:2:end)];

                        % x=data(count*NumVals+1:(count+1)*NumVals,:);
                        % numel (x(:,1))/NumVals
                        for v = 0:N
			    % 1 timestep = 3*N seconds. divide by 60*60 = 3600 to convert to hours
                            xv = x(x(:,2)==v,1)*dt*N./(60*60);% normalize data to be within hours
                            FractVotes = numel (xv)./NumVals;

			    % if timesteps are larger than bins, we smooth by the timestep, dt/(60 60)
			    if dt/(60*60) > (Tboundaries(2)-Tboundaries(1))
				y(b,mvm,tau,tau_hung,a,a_hung,v+1,:) = FractVotes.*ksdensity(xv,t,'support','positive','function','cdf','bandwidth',dt/(60*60));
			    else
                                % t
                                h2 = histc(xv, Tboundaries);
	    	 	        % 5 minute steps: smooth out CDF/histogram
                                h2(end) = [];
                                % Compute distance using KS statistic
                                temp = reshape(cumsum(h2(end:-1:1)),1,numel(h2));
                            
                                if numel(xv) > 0
                                    cx = [1-temp (end:-1:1)./numel(xv) 1];
                                else
                                    cx = zeros(1,numel(Tboundaries));
                                end
                                y(b,mvm,tau,tau_hung,a,a_hung,v+1,:) = FractVotes.*cx;% ksdensity(max([exprnd(i_guilty,j,NumTrials);exprnd(i_innocent,N-j,NumTrials)]),t,'support','positive','function','cdf');
                             end
			end
                        count = count + 1;
                    end
                end
            end
        end
    end
end
fclose(fid);


end
