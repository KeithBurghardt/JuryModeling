% make alpha vs Ttrial figure

%Data:
clear;
%load('Influence_Time_Q0=0.3.mat');

x=1.8.^(0.5:8.5);
%WA6Est=all_param_est(1:9,1,5);
%NE12Est=all_param_est(10:end,1,5);
%WA6Errors=reshape(all_param_quantiles(1:9,5,2:3),9,2);
%NE12Errors=reshape(all_param_quantiles(10:end,5,2:3),9,2);


load('Influence_FBinom_FullModel_Q0=0.3.mat');

CA8Est=all_param_est(2+9+1:2+9+9,1,5);
CA12Est=all_param_est(2+9+9+1:end,1,5);
CA8Errors=reshape(all_param_quantiles(2+9+1:2+9+9,5,2:3),9,2);
CA12Errors=reshape(all_param_quantiles(2+9+9+1:end,5,2:3),9,2);

%Taking parameter resolution into account when calvulating errors
%CA8Errors(4,:) = [0.013 0.019];
%CA8Errors(5,:) = [0.008 0.013];
%New changes made
CA8Est(4)=0.0154;
CA8Est(5)=0.0104;
CA8Errors(4,:)=[0.014 0.017];
CA8Errors(5,:)=[0.010 0.012];

CA12Est(4)=0.0146;
CA12Est(5)=0.011;
CA12Est(6)=0.0079;
CA12Est(7)=0.0068;
CA12Errors(4,:)=[0.014 0.016];
CA12Errors(5,:)=[0.010 0.013];
CA12Errors(6,:)=[0.007 0.009];
CA12Errors(7,:)=[0.005 0.007];
CA12Errors(8,:)=[0.0035 0.0045];

%set up figure

figure;
plot(1,2);
set(gca,'FontSize',20); 
set(gca,'xscale','log','yscale','log');
ax=gca;
axis([3 100 3*10^-3 10^-1.5]);
box on;
xlabel('T_{trial} (Hrs)','FontSize',25);
ylabel('$$\hat{\alpha}$$','Interpreter','Latex','FontSize',25);
ax.TickLength = [0.02, 0.04];
ax.YTick=[0.003 0.01 0.03];
ax.YTickLabel={'3\times10^{-3}','1\times10^{-2}','3\times10^{-2}'};
hold on;

%best fit line
plot(x,0.04*x.^(-1/2),'-','LineWidth',6,'Color',[0.5 0.5 0.5]);


%CA8
for i=4:6;
    gray=0.95-0.2*(i-4);
    plot(x(i),CA8Est(i),'o','MarkerFaceColor',[gray gray gray],'Color',[0 0.7 0],'LineWidth',2,'MarkerSize',20);
    CA8plot=errorbar(x(i),CA8Est(i),CA8Errors(i,1)-CA8Est(i),CA8Errors(i,2)-CA8Est(i),'o','Color',[0 0.7 0],'LineWidth',2,'MarkerSize',20);
end;

%CA12
for i=4:8;
    gray=0.95-0.2*(i-4);
    plot(x(i),CA12Est(i),'^','MarkerFaceColor',[gray gray gray],'Color',[0 0.7 0.7],'LineWidth',2,'MarkerSize',20);
    CA12plot=errorbar(x(i),CA12Est(i),CA12Errors(i,1)-CA12Est(i),CA12Errors(i,2)-CA12Est(i),'^','Color',[0 0.7 0.7],'LineWidth',2,'MarkerSize',20);

end;


%WA6
%for i=3:5;
%    if i==3;
%        gray = 1;
%    else
%        gray=0.95-0.2*(i-4);
%    end
%    plot(x(i),WA6Est(i),'pk','MarkerFaceColor',[gray gray gray],'Color',[0 0 0],'LineWidth',2,'MarkerSize',20);
%    WAplot=errorbar(x(i),WA6Est(i),WA6Errors(i,1)-WA6Est(i),WA6Errors(i,2)-WA6Est(i),'pk','LineWidth',2,'MarkerSize',20);
%end;
%NE12
%for i=4:6;
%    if i==3;
%        gray = 1;
%    else
%        gray=0.95-0.2*(i-4);
%    end
%    NEplot=errorbar(x(i),NE12Est(i),NE12Errors(i,1)-NE12Est(i),NE12Errors(i,2)-NE12Est(i),'x','Color',[0.85 0.5 0],'LineWidth',3,'MarkerSize',20);
%end;
