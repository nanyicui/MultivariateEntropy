function GaussianEstimate_200(XMatrix1, YMatrix1, EMatrix1, X1, Y1, YMatrix2, EMatrix2, YMatrix3, EMatrix3, YMatrix4, EMatrix4)
%CREATEFIGURE(XMATRIX1,YMATRIX1,EMATRIX1,X1,Y1,YMATRIX2,EMATRIX2,YMATRIX3,EMATRIX3,YMATRIX4,EMATRIX4)
%  XMATRIX1:  errorbar x matrix
%  YMATRIX1:  errorbar y matrix
%  EMATRIX1:  errorbar e matrix
%  X1:  vector of x data
%  Y1:  vector of y data
%  YMATRIX2:  errorbar y matrix
%  EMATRIX2:  errorbar e matrix
%  YMATRIX3:  errorbar y matrix
%  EMATRIX3:  errorbar e matrix
%  YMATRIX4:  errorbar y matrix
%  EMATRIX4:  errorbar e matrix

%  Auto-generated by MATLAB on 28-Feb-2012 16:24:10

% Create figure
figure1 = figure;

% Create subplot
subplot1 = subplot(1,2,2,'Parent',figure1);
% Uncomment the following line to preserve the X-limits of the axes
 xlim(subplot1,[4 13]);
% Uncomment the following line to preserve the Y-limits of the axes
 ylim(subplot1,[1.6 2.8]);
hold(subplot1,'all');

% Create multiple error bars using matrix input to errorbar
errorbar1 = errorbar(XMatrix1,YMatrix1,EMatrix1,'MarkerSize',8,...
    'MarkerFaceColor',[0 0 1],...
    'LineStyle','none');
set(errorbar1(1),'Marker','^','DisplayName','binless');
set(errorbar1(2),'Marker','o','DisplayName','binned');

% Create plot
plot(X1,Y1,'Parent',subplot1,'LineWidth',2,'LineStyle',':',...
    'DisplayName','correct value');

% Create xlabel
xlabel('log_{2}(Sample Size)');

% Create ylabel
ylabel('H_{diff} (bits)');

% Create title
%title('Number of bins = Sample Size / 16');

% % Create subplot
% subplot2 = subplot(2,2,3,'Parent',figure1);
% % Uncomment the following line to preserve the X-limits of the axes
%  xlim(subplot2,[4 13]);
% hold(subplot2,'all');
% 
% % Create multiple error bars using matrix input to errorbar
% errorbar2 = errorbar(XMatrix1,YMatrix2,EMatrix2,'MarkerSize',8,...
%     'MarkerFaceColor',[0 0 1],...
%     'LineStyle','none');
% set(errorbar2(1),'Marker','^','DisplayName','binless');
% set(errorbar2(2),'Marker','o','DisplayName','binned');
% 
% % Create plot
% plot(X1,Y1,'Parent',subplot2,'LineWidth',2,'LineStyle',':',...
%     'DisplayName','correct value');
% 
% % Create xlabel
% xlabel('log_2(Sample Size)');
% 
% % Create ylabel
% ylabel('H_{diff} /bits');
% 
% % Create title
% title('Number of bins = Sample Size / 8');

% Create subplot
subplot2 = subplot(1,2,1,'Parent',figure1);
% Uncomment the following line to preserve the X-limits of the axes
 xlim(subplot2,[4 13]);
hold(subplot2,'all');

% Create multiple error bars using matrix input to errorbar
errorbar2 = errorbar(XMatrix1,YMatrix3,EMatrix3,'MarkerSize',8,...
    'MarkerFaceColor',[0 0 1],...
    'LineStyle','none');
set(errorbar2(1),'Marker','^','DisplayName','binless');
set(errorbar2(2),'Marker','o','DisplayName','binned');

% Create plot
plot(X1,Y1,'Parent',subplot2,'LineWidth',2,'LineStyle',':',...
    'DisplayName','correct value');

% Create xlabel
xlabel('log_{2}(Sample Size)');

% Create ylabel
ylabel('H_{diff} (bits)');

% Create title
%title('Number of bins = Sample Size / 4');

% % Create subplot
% subplot4 = subplot(2,2,1,'Parent',figure1);
% % Uncomment the following line to preserve the X-limits of the axes
%  xlim(subplot4,[4 13]);
% hold(subplot4,'all');
% 
% % Create multiple error bars using matrix input to errorbar
% errorbar4 = errorbar(XMatrix1,YMatrix4,EMatrix4,'MarkerSize',8,...
%     'MarkerFaceColor',[0 0 1],...
%     'LineStyle','none');
% set(errorbar4(1),'Marker','^','DisplayName','binless');
% set(errorbar4(2),'Marker','o','DisplayName','binned');
% 
% % Create plot
% plot(X1,Y1,'Parent',subplot4,'MarkerFaceColor',[0 0 1],'MarkerSize',8,...
%     'LineWidth',2,...
%     'LineStyle',':',...
%     'DisplayName','correct value');
% 
% % Create xlabel
% xlabel('log_2(Sample Size)');
% 
% % Create ylabel
% ylabel('H_{diff} /bits');
% 
% % Create title
% title('Number of bins = Sample Size / 2');

% Create legend
legend(subplot2,'show');

% Create legend
legend(subplot1,'show');

% % Create legend
% legend(subplot2,'show');
% 
% % Create legend
% legend(subplot4,'show');
