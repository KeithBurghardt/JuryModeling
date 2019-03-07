function [BinomialRange,MVMRange,TauRange,TauHungRange,AlphaRange,AlphaHungRange,N,NumVals] = GetRanges(file)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%     Find ranges from file     %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %string positions where "number of jurors" will be listed
    disp(file);
    pos1=strfind(file,'N=');
    pos2=strfind(file,';Binomial_p='); 
    pos3=strfind(file,';MVM_p=');
    pos4=strfind(file,';Alpha=');
    pos5=strfind(file,';Alpha_hung=');
    pos6=strfind(file,';Tau=');
    pos7=strfind(file,';Tau_hung=');
    if numel(pos7) == 0
        pos7=strfind(file,';HungRatio=');
    end

    pos8=strfind(file,';NumRuns=');
    pos9=strfind(file,';NumTrials=');
    disp([pos1 pos2 pos3 pos4 pos5 pos6 pos7 pos8 pos9]);

    N = str2num(file(pos1+2:pos2-1));
    NumRuns = str2num(file(pos8+9:pos9-1));
    NumTrials = str2num(file(pos9+11:pos9+12));
    NumVals = NumRuns * NumTrials;

    %below order is the order in which the data is written
    BinomialRange=str2num(strrep(file(pos2+12:pos3-1),'_',':'));
    MVMRange=str2num(strrep(file(pos3+7:pos4-1),'_',':'));
    TauRange=str2num(strrep(file(pos6+5:pos7-1),'_',':'));
    TauHungRange=str2num(strrep(file(pos7+10:pos8-1),'_',':'));
    if numel(strfind(file,';HungRatio=')) > 0
	TauHungRange=str2num(strrep(file(pos7+11:pos8-1),'_',':'));
    end
    AlphaRange=str2num(strrep(file(pos4+7:pos5-1),'_',':'));
    AlphaHungRange=str2num(strrep(file(pos5+12:pos6-1),'_',':'));
end
