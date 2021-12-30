%% Data extraction
clearvars
close all


alpha =  .04; % learning rate
beta =   15; % 1/Boltzmann temperature,0: random, large:rational
go_bias = .1;


Nmax = 0.8e3; % number of trials
Nsessions = 1;
Policy='softmax';%'greedy';%'softmax'


Nodor = 4;

correct_output = nan(Nsessions,Nmax,Nodor);

QInit = [.2.*ones(Nodor,1) , 0.*ones(Nodor,1)];


[file, path] = uigetfile('*.*','MultiSelect','on');
if ischar(file)
    file = cellstr(file); %if only one file was loaded, converts to cell
end

file = sort(file); %sorts files in alphabetical order

movavg = 10; %number of elements to average when plotting


for k = length(file):-1:1
    bdata = load([path file{k}]);
    
    if size(bdata(:,2),1) <300 || range(bdata(:,2)) < 3
        data(k) = [];
        continue
    else
        data(k).file = file{k}; %saving filename as string
        data(k).animal = str2num(file{k}(2:3)); %animal number
        data(k).date = data(k).file(9:12);
        
        
        data(k).odorID = bdata(:,2); data(k).odorID(bdata(:,3)==0 & ...
            bdata(:,4)==0) =[]; % odor number in sequence; removes trials that
        % were neither correct or incorrect (incomplete trials)
        data(k).odorID = padarray(data(k).odorID,[1500-length(data(k).odorID) 0],nan,'post');
        %padding array to allow for concatenation later
        %         data(k).trial = (1:size(data(k).odorID,1))'; % matrix of trial numbers
        %         data(k).trial = padarray(data(k).trial,[1500-length(data(k).trial) 0],nan,'post');
        data(k).allcorrect = bdata(:,3); data(k).allcorrect(bdata(:,3)==0 & ...
            bdata(:,4)==0) = []; %full session correct or not; removes incomplete trials
        data(k).allcorrect = padarray(data(k).allcorrect,[1500-length(data(k).allcorrect) 0],nan,'post');
        data(k).rxn = bdata(:,9) - bdata(:,8); data(k).rxn(bdata(:,3)==0 & ...
            bdata(:,4)==0) = []; %calculates reaction speed for full session
        data(k).rxn = padarray(data(k).rxn,[1500-length(data(k).rxn) 0],nan,'post');
        data(k).trialstart = bdata(:,7);
        data(k).trialstart(bdata(:,3)==0 & ...
            bdata(:,4)==0) =[];
        data(k).trialstart = padarray(data(k).trialstart,[1500-length(data(k).trialstart) 0],nan,'post');
        
        % Preparing data for model fit
        
        state = data(k).odorID;
        
        state(data(k).odorID == 2) = 3; state(data(k).odorID == 3) = 2; %switches
        %odors 2 & 3 to match model identities
        
        correct = data(k).allcorrect;
        
        action = zeros(size(state)); % 1 = mouse went to water port; 0 = it did not
        action(state==1 & correct==1) = 1;
        action(state==2 & correct==1) = 1;
        action(state==3 & correct==0) = 1;
        action(state==4 & correct==0) = 1;
        
        reward = zeros(size(state));
        reward(correct ==1 & (state == 1 | state ==2)) = 1; % when mouse receives water
        
        opts = optimset('display','iter', 'PlotFcns',{@optimplotx,...
            @optimplotfval},'MaxFunEvals',1e4,'TolX',1e-6); % debugging version
        opts_ = optimset('display','off','MaxFunEvals',1e3,'TolX',1e-6); % simple version
        
        Nrepeats = 1; %number of fittings, initial the parameters at random start point
        para_names = {'\alpha','\beta','go bias'};
        para_truth = [alpha, beta,go_bias];
        Nparas = length(para_truth);
        Para_opt = nan(Nparas,Nrepeats,Nsessions);
        
        
        for j = 1:4
            data(k).odor(j).rxn = ...
                data(k).rxn(data(k).odorID == j); % rxn times of odor j
            
            data(k).odor(j).correct = ...
                data(k).allcorrect(data(k).odorID ==j); %correct trials of odor j
            
            data(k).odor(j).trialstart = data(k).trialstart(data(k).odorID == j);
            
            
            data(k).odor(j).tlearned = [];
            data(k).odor(j).performance = mean(data(k).odor(j).correct,'omitnan');
            data(k).odor(j).finalperf = mean(data(k).odor(j).correct(end-19:end),'omitnan');
            data(k).odor(j).slidmean = data(k).odor(j).correct(~isnan([data(k).odor(j).correct]));
            data(k).odor(j).slidmean = movmean(data(k).odor(j).slidmean,20);
            data(k).odor(j).stdev_rxn = std(data(k).odor(j).rxn,'omitnan');
            data(k).odor(j).rxn_avg = mean(data(k).odor(j).rxn,'omitnan');
            
            goodorlast20 = mean(data(k).odor(1).slidmean(end-19:end) >= 0.9) ==1;
            
            if not(goodorlast20)
                truth90 = false(1,length(data(k).odor(1).slidmean));
                truth90(data(k).odor(1).slidmean>=.85) = true;
                last90 = find(not(truth90),1,'first')-1;
                if last90 == 0
                    last90 =1;
                end
                last90_time = data(k).odor(1).trialstart(last90);
                data(k).odor(j).slidmean = data(k).odor(j).slidmean(1:find(data(k).odor(j).trialstart<last90_time,1,'last'));
                data(k).odor(j).correct = data(k).odor(j).correct(1:find(data(k).odor(j).trialstart<last90_time,1,'last'));
                data(k).odor(j).rxn = data(k).odor(j).rxn(1:find(data(k).odor(j).trialstart<last90_time,1,'last')); % removes end of session if it...
                %dropped below 90 (program error or loss of motivation)
                
                
            end
            
            
            for i = 1:size(data(k).odor(j).correct,1)-2 %finding when 3 consecutive trials are first done
                if data(k).odor(j).correct(i) == 1 && ...
                        data(k).odor(j).correct(i+1) ==1 && ...
                        data(k).odor(j).correct(i+2) ==1
                    data(k).odor(j).tlearned = i; % number of trials before 3
                    % consecutive trials have been done correctly
                    break %ends for-loop once 3 consecutive correct trials have been found
                end
                if isempty(data(k).odor(j).tlearned)
                    data(k).odor(j).tlearned = inf; % if odor was never learned, sets to infinity
                end
            end
            if (j == 2 | j == 4) & data(k).odor(j).tlearned ~= inf
                finalperf90 = movmean(data(k).odor(j).correct,10,'omitnan') > 0.9*data(k).odor(j).finalperf;
                data(k).odor(j).tlearned_to_90 = find(finalperf90(data(k).odor(j).tlearned:end),1,'first') ...
                    + (data(k).odor(j).tlearned-1) - data(k).odor(j).tlearned;
            elseif (j == 2 | j == 4) & data(k).odor(j).tlearned == inf
                data(k).odor(j).tlearned_to_90 = inf;
            end
            if length(data(k).odor(j).correct)<20 || mean(data(k).odor(j).correct(end-19:end),'omitnan') < 0.5
                data(k).odor(j).tlearned = inf;
            end
            
            
            data(k).odor(j).slidmean = padarray(data(k).odor(j).slidmean,...
                [300-length(data(k).odor(j).slidmean) 0],nan,'post');
            data(k).odor(j).correct = padarray(data(k).odor(j).correct,...
                [300-length(data(k).odor(j).correct) 0],nan,'post');
            data(k).odor(j).rxn = padarray(data(k).odor(j).rxn,...
                [300-length(data(k).odor(j).rxn) 0],nan,'post');
            data(k).odor(j).trialstart = padarray(data(k).odor(j).trialstart,...
                [300-length(data(k).odor(j).trialstart) 0],nan,'post');
            
            cross60 = find([data(k).odor(j).slidmean] >= 0.6);
            if isempty(cross60)
                cross60 = inf;
                data(k).odor(j).cross60 = cross60(1);
            end
            data(k).odor(j).cross60 = cross60(1);
            
            
        end
        
        if not(goodorlast20)
            truth90 = false(1,length(data(k).odor(1).slidmean));
            truth90(data(k).odor(1).slidmean>=.85) = true;
            last90 = find(not(truth90),1,'first')-1;
            if last90 == 0
                last90 =1;
            end
            last90_time = data(k).odor(1).trialstart(last90);
            data(k).odorID(data(k).trialstart>last90_time) = [];
            data(k).allcorrect(data(k).trialstart>last90_time) = [];
            data(k).rxn(data(k).trialstart>last90_time) = [];
            data(k).trialstart(data(k).trialstart>last90_time) = [];
            data(k).odorID = padarray(data(k).odorID,[1500-length(data(k).odorID) 0],nan,'post');
            data(k).allcorrect = padarray(data(k).allcorrect,[1500-length(data(k).allcorrect) 0],nan,'post');
            data(k).rxn = padarray(data(k).rxn,[1500-length(data(k).rxn) 0],nan,'post');
            data(k).trialstart = padarray(data(k).trialstart,[1500-length(data(k).trialstart) 0],nan,'post');
        end
        
        data(k).nogodelay = data(k).odor(4).tlearned - data(k).odor(2).tlearned;
        data(k).perfdiff = data(k).odor(2).performance - data(k).odor(4).performance;
        data(k).cross60diff = data(k).odor(4).cross60 - data(k).odor(2).cross60;
        
        
        if isempty(data(k).nogodelay) || isnan(data(k).nogodelay) || data(k).odor(1).performance < 0.9 ...
                || data(k).odor(3).performance < 0.9 || length(data(k).odor(4).correct(~isnan(data(k).odor(4).correct))) < ...
                0.6*length(data(k).odor(1).correct(~isnan(data(k).odor(1).correct)))
            data(k) = [];
        end
    end
end

% adding stim

[stim,date] = xlsread(['stimdate' num2str(data(1).animal) '.xlsx']);

for i = 1:length(stim)
    if sum(strcmp({data.date},date(i)))>0
        data(strcmp({data.date}, date(i))).stim = stim(i);
    end
end
data_backup = data;
%% Plot moving average curves

color = {'k-','r-','c--','m--'}; %color/line properties for plotting

for k = 1:length(data)
    figure(k)
    for j = 1:4
        hold on
        plot(movmean(data(k).odor(j).correct,movavg),color{j},'LineWidth',1.5)
        plotbrowser on
    end
    
    legend('Odor 1 (Go)','Odor 2 (No-Go)','Odor 3 (Go w/ Stim)',...
        'Odor 4 (No-Go w/ Stim')
    title(['File: ' data(k).file],'Interpreter','none')
    xlabel('Trial Number')
    ylabel(sprintf('Moving Average Correct (%.f)',movavg))
    text(0.05,0.87,{'\bf Time Until',' Learned:\rm', ...
        sprintf(' Odor 1: %.f \n Odor 2: %.f \n Odor 3: %.f \n Odor 4: %.f',...
        data(k).odor(1).tlearned,data(k).odor(2).tlearned,...
        data(k).odor(3).tlearned,data(k).odor(4).tlearned),...
        sprintf(' No-go Delay: %.f \n',data(k).nogodelay)},...
        'FontSize',14,'Color','b','Interpreter','tex','VerticalAlignment','top','Units','normalized')
    text(0.16,0.87,{'\bf Optimal',' Parameters: \rm',...
        sprintf(' alpha: %.3f \n beta: %.3f \n go bias: %.3f',...
        Para_opt(:,1,1))},'FontSize',14,'Color','b','Interpreter','tex','VerticalAlignment','top',...
        'Units','normalized')
    ylim([0,1.1])
    grid on
    plotbrowser on
    
    %     saveas(figure(2*k),['Curve_' file{k}(1:end-4) '.png'])
    
end
%% Aligned binary curves
% extracts curves into 3D matrix: session x trial x odor (row, col, page)
% curves are aligned at trial learned (3 trial criteria), column 300 is
% first trial correct; or last trial if never learned

data(isnan([data.nogodelay])) = [];
data = data([data.stim] ==1);

odorID = [data.odorID]';
correct = [data.allcorrect]';

curve = nan(length(data),600,4);
%curve = correct(:,odorID == 2);

for j = 1:length(data)
    for i = 1:4
        if data(j).odor(i).tlearned == inf
            %%tlearned = size(data(j).odor(i).correct(~isnan(data(j).odor(i).correct)),1);
        else
            tlearned = data(j).odor(i).tlearned;
            curve(j,300-tlearned+1:300,i) = data(j).odor(i).correct(1:tlearned);
            curve(j,301:300+(300-tlearned),i) = data(j).odor(i).correct(tlearned+1:end);
        end
        
    end
end

figure('Position',[0 0 500 800])
for i = 1:4
    [row,col] = find(~isnan(curve(:,:,i)));
    
    subplot(4,1,i)
    plot(movmean(curve(:,min(col):end,i)',10,1,'omitnan'))
    %plot(movmean(curve(:,280:end,i)',10,1,'omitnan'))
end
%% Q-Learning: fitting across sessions; separated by stimulation
% each session has it's own cost fucntion. fminsearch uses a costfunction
% that is the sum of all other cost functions
animal = 'Striatum';

data(isnan([data.nogodelay])) = [];
data = data([data.stim] ==1);
states = cell(1,length(data));
actions = cell(1,length(data));
rewards = cell(1,length(data));
correct = cell(2,length(data));

Nsessions = 1;

for k = 1:length(data)
    for l = 1:2
        
        if abs(data(k).odor(l+2).tlearned) ~= inf && abs(data(k).odor(l).tlearned) ~= inf
            
            states{l,k} = data(k).odorID;
            
            %correct = data(k).allcorrect;
            
            states{l,k}(data(k).odorID == 2) = 3; states{l,k}(data(k).odorID == 3) = 2; %switches
            %odors 2 & 3 to match model identities
            correct{l,k} = data(k).allcorrect(states{l,k} == l | states{l,k} == l+2);
            states{l,k}(states{l,k} ~= l & states{l,k} ~= l+2) = [];
            
            actions{l,k} = zeros(size(states{l,k})); % 1 = mouse went to water port; 0 = it did not
            actions{l,k}(states{l,k}==l & correct{l,k}==1) = 1;
            %actions{k}(states{k}==2 & correct==1) = 1;
            actions{l,k}(states{l,k}==l+2 & correct{l,k}==0) = 1;
            %actions{k}(states{k}==4 & correct==0) = 1;
            
            rewards{l,k} = zeros(size(states{l,k}));
            rewards{l,k}(correct{l,k} ==1 & (states{l,k} == 1 | states{l,k} ==2)) = 1; % when mouse receives water
            
            states{l,k}(isnan(states{l,k})) = [];
            actions{l,k} = actions{l,k}(1:length(states{l,k}));
            rewards{l,k} = rewards{l,k}(1:length(states{l,k}));
            
        end
    end
end
close all;clc;
opts = optimset('display','iter', 'PlotFcns',{@optimplotx,...
    @optimplotfval},'MaxFunEvals',1e4,'TolX',1e-6); % debugging version
opts_ = optimset('display','off','MaxFunEvals',1e3,'TolX',1e-6); % simple version

Nrepeats = 10; %number of fittings, initial the parameters at random start point
para_names = {'\alpha','\beta','go bias'};
para_truth = [alpha, beta,go_bias];
Nparas = length(para_truth);
Para_opt = nan(Nparas,Nrepeats,Nsessions);
Para_init = [];

for kkk = 1:2
    for jjj=1:Nsessions
        for iii=1:Nrepeats
            Para_init(:,iii,jjj,kkk) = [0.5*rand;15*rand+15;0.4*rand];%randn(size(para_truth));%abs(para_truth.*(1 + .5*randn(size(para_truth))));
            [Para_opt(:,iii,jjj,kkk),NLL] = fminsearch(@(Para)SUM_GNG_loglikeli_action(Para,states(kkk,:),...
                actions(kkk,:),rewards(kkk,:), Policy,QInit),Para_init(:,iii,jjj,kkk),opts_);
            
        end
    end
    
end


titles = {'No-Stim Odors','Stim Odors'};
figure('Name','Stim Vs UnStim','Position',[0 500 1500 400]);
colors = {'#0072BD','#D95319'};
colorbox = {'b','r'};
ytextpos = [0.95 0.9];
t = tiledlayout(1,3);
t.TileSpacing = 'compact';
t.Padding = 'compact';

for ppp=1:Nparas
    ylims = {[0,0.15],[0,25],[-0.05,0.25]};
    nexttile
    for lll = 1:2
        hold on
        %         subplot(1,Nparas,ppp); hold on
        %bar(0,para_truth(ppp),'FaceColor','k','EdgeColor','none','FaceAlpha',.3 );
        para = reshape(squeeze(Para_opt(ppp,:,:,lll)),[],1);
        plot(squeeze(Para_init(ppp,:,:,lll)),para,'o','Color',colors{lll},'LineWidth',4,'MarkerSize',10);
        paramean(lll,ppp) = nanmean(para);
        parastd(lll,ppp) = nanstd(para);
        errorbar(nanmean(squeeze(Para_init(ppp,:,:,lll))),nanmean(para),nanstd(para),'o','LineWidth',5,'MarkerFaceColor','k','CapSize',0,'Color',colors{lll},'MarkerSize',10)
        ax = gca;
        ax.FontSize = 35;
        ax.FontName = 'Times New Roman';
        ylabel(para_names{ppp})
        xticks('auto')
        xticklabels('auto')
        ylim(ylims{ppp})
        %         xlabel('Initialized Value')
        %         text(0.1,ytextpos(lll),['\textbf{Average $' para_names{ppp} ' =$' num2str(nanmean(para),3) '}'],'Units','normalized','FontSize',14,...
        %             'Interpreter','latex','Color',colors{lll})
        %         text(0.1,ytextpos(lll)-0.1,['\textbf{$ \sigma $ = ' num2str(nanstd(para),3) '}'],'Units','normalized','FontSize',14,...
        %             'Interpreter','latex','Color',colors{lll})
        %         title(['\textbf{ ' sprintf('D%.f %s Optimized Parameter: $ %s $',data(1).animal,titles{lll},para_names{ppp}) '}'],'FontSize',16,'Interpreter','latex')
        title([para_names{ppp}])
        
    end
end
xlabel(t,'Initialized Value','FontSize',30,'FontName','Times New Roman')
% subplot(1,3,2)
lgd = legend({'Unstim model fit','Unstim mean+/-sd','Stim model fit','Stim mean+/-sd'},'Location','southoutside','FontSize',30);
lgd.Layout.Tile = 'east';
Parameter = {' ','$\alpha$ unstim','$\alpha$ stim','$\beta$ unstim','$\beta$ stim','$gobias$ unstim','$gobias$ stim'};
Mean = round([paramean(:,1)',paramean(:,2)',paramean(:,3)'],3,'significant');
STD = round([parastd(:,1)',parastd(:,2)',parastd(:,3)'],3,'significant');


t = num2cell([Mean;STD]);
t(:,2:7) = t;
t(:,1) = {'Mean';'$\sigma$'};
t = array2table(t,'VariableNames',Parameter);
table2latex(t,'test')
%% Q-Learning simulation

%NOTE: Run aligned binary curve section as well as fitting across sessions
%section before running this!!!

close all;
% Parameters
% clc;  clear all;
animal = 'Striatum No-Go';

rxntimeplot = true; %switch between true/false to plot rxn time

alpha_1 =  nanmean(Para_opt(1,:,:,1),'all'); % learning rate of unstim
alpha_2 =  nanmean(Para_opt(1,:,:,2),'all'); % learning rate of stim

beta_1 =  nanmean(Para_opt(2,:,:,1),'all'); % 1/Boltzmann temperature,0: random, large:rational
beta_2 =  nanmean(Para_opt(2,:,:,2),'all');

go_bias_1 =  nanmean(Para_opt(3,:,:,1),'all');
go_bias_2 =  nanmean(Para_opt(3,:,:,2),'all');

epsilon = 0.2;  % cost of Go, relative to drop of water reward
gamma =  .2; % laser stim,
epsilon_noise = 0.00; % noise in updating Q

zeta =1; %0.5; %for punishing learning rate on stimulated trials

%rng('default');

sigm = @(x) exp(x)./(1 + exp(x)); % sigmoid function
dsigm = @(x) exp(x)./(1 + exp(x)).^2; % derivative of sigmoid

Nmax = 110*4; % number of trials
Nsessions = 8;
Nmean = 10; % running mean for calculating %correct
Policy='softmax';%'greedy';%'softmax'
Trial_Type = {'GO';'GO-stim';'NOGO';'NOGO-stim'};
Odor = [[1,0,0,0];...% GO NOGO GO-stim NOGO-stim
    [0,1,0,0];...
    [0,0,1,0];...
    [0,0,0,1]];

Correct = [ [1];...% correct action
    [1];...
    [0];...
    [0]];
Reward_ = [[1-epsilon,0];...% Go
    [1-epsilon,0];...% GO-stim
    [ -epsilon,0];...      % NOGO
    [ -epsilon,0]];  % NOGO-stim

Nodor = size(Odor,1);
CPD = cumsum(1/Nodor.* ones(Nodor,1));
maxRT = nan(Nodor,Nsessions);

State = nan(Nmax,Nsessions); % Odor state
Action = nan(Nmax,Nsessions); %
Reward = nan(Nmax,Nsessions); % recieved eward
RT = nan(Nmax,Nsessions);

correct_output = nan(Nsessions,Nmax,Nodor);

QInit = [.2.*ones(Nodor,1) , 0.*ones(Nodor,1)];%+ epsilon_noise .*randn(Nodor,2) ; Q-initialization?
%QInit(4,:) = [0.8 0];

if rxntimeplot
    
    figure('Position',[0 0 600 800])
    
else
    
    figure('Position',[0 0 500 800])
    
end
for iii=1:Nsessions
    % initialization
    [~,~,Od] = histcounts(rand(1,Nmax),CPD);
    Od = Od +1; %Odor identity, random odor in each trial
    
    switch Policy
        case 'softmax'
            w_ = QInit;
            w  = QInit;
            Q = nan([size(w_),Nmax]);
            
    end
    action = nan(1,Nmax);
    reward = nan(1,Nmax);
    % learning
    for i=1:Nmax
        %w_ = w_ + epsilon_noise .*randn(size(w_)); % introduce noise
        switch Policy
            case 'softmax'
                odor = Odor( Od(i) , :);
                x_go = odor*w_(:,1)  ; x_nogo = odor*w_(:,2);
                if Od(i) ==2 || Od(i) == 4
                    alpha = alpha_2;
                    beta = beta_2;
                    go_bias = go_bias_2;
                else
                    alpha = alpha_1;
                    beta = beta_1;
                    go_bias = go_bias_1;
                end
                prob_go = sigm((x_go - x_nogo + go_bias).*beta); %
                RT(i,iii) = 1./abs((x_go - x_nogo));
                a = rand(1)< prob_go;
                action(i) = a;
                if a ==1; indx_a = 1; else ;indx_a = 2; end
                if a
                    sign_a = 1;  q = (x_go);  dq = dsigm(x_go);
                else
                    sign_a = -1; q = (x_nogo); dq = dsigm(x_nogo);
                end
                
                r = Reward_(Od(i),indx_a) ;
                reward(i) = r;
                
                
                dw = alpha.*(r-q);
                I = epsilon_noise .*(randn(size(w_)));% introduce noise
                I(Od(i),indx_a) = I(Od(i),indx_a) + dw;
                w = w_ + I;
                Q(:,:,i) = w; %  Q value , Nodor x Naction x Ntrial
                w_ = w;
                
        end
    end
    State(:,iii) = Od;
    Action(:,iii) = action;% 0: no go or 1:go
    Reward(:,iii) = reward;
    
    
    
    % plotting
    for j=1:Nodor
        jj = j;
        clear indx correct_action correct
        indx = (Od==j);
        correct_action = Correct(j);
        correct = double(action(indx) == correct_action);
        correct_output(iii,1:size(correct,2),j) = correct;
        
    end
    
    
end
clf
color = {'b','r'};
for j = 3:4
    
    
    if rxntimeplot
        nplotcol = 2; %number of plot columns
        fun1 = @(j) 2*j-1;
        fun2 = @(j) 2*j;
    else
        nplotcol =1;
        fun1 = @(j) j;
        fun2 = @(j) j;
    end
    hold on
    
    subplot(2,1,1)
    y = mean(movmean(correct_output(:,:,j)',10,1,'omitnan'),2,'omitnan');
    %         y = movmean(mean(correct_output(:,:,j),1,'omitnan'),10,2,'omitnan');
    err = std(movmean(correct_output(:,:,j)',10,1,'omitnan'),0,2,'omitnan');
    h = shadedErrorBar([],y,err,3,'lineProps',{color{j-2},'markerfacecolor',color{j-2}});
    %         text((Nmax/Nodor/2),.3,{['$\alpha = $', num2str(alpha_1,3)];['$\beta = $', num2str(beta_1,3)];...
    %             ['Go bias = ', num2str(go_bias_1,3)]},'Color',color{1},'FontSize',15,'Interpreter','latex');
    %         text((Nmax/Nodor/2+50),.3,{['$\alpha = $', num2str(alpha_2,3)];['$\beta = $', num2str(beta_2,3)];...
    %             ['Go bias = ', num2str(go_bias_2,3)]},'Color',color{2},'FontSize',15,'Interpreter','latex');
    
    ax = gca;
    ax.FontSize = 25;
    ax.FontName = 'Times New Roman';
    title('Q-Learning Model');
    xlabel('Trials');
    ylabel('Performance')
    ylim([-.1 1.1])
    hold on
    subplot(2,1,2)
    if j == 2
        jj = 3;
    elseif j ==3
        jj = 2;
    else
        jj = j;
    end
    hold on
    if jj == 1 || jj == 3
        y = mean(movmean(curve(:,300:end,jj)',10,1,'omitnan'),2,'omitnan');
        err = std(movmean(curve(:,300:end,jj)',10,1,'omitnan'),0,2,'omitnan');
        h = shadedErrorBar([],y,err,3,'lineProps',{color{j-2},'markerfacecolor',color{j-2}});
    else
        probcurve = mean(curve,1,'omitnan');
        y = movmean(probcurve(:,280:280+140,jj),20,1,'omitnan');
        
        err = std(movmean(movmean(curve(:,280:280+140,jj)',10,1,'omitnan'),10,1,'omitnan'),0,2,'omitnan');
        h = shadedErrorBar([],movmean(y,10),err,3,'lineProps',{color{j-2},'markerfacecolor',color{j-2}});
    end
    ax = gca;
    ax.FontSize = 25;
    ax.FontName = 'Times New Roman';
    
    ylim([-.1 1.1])
    title('Behavioral Data');
    xlabel('Trials');
    ylabel('Performance');
    
end
subplot(2,1,1)
legend({'Unstimulated','Stimulated'},'FontSize',20,'Location','southeast')
subplot(2,1,2)
legend({'Unstimulated','Stimulated'},'FontSize',20,'Location','southeast')
%% Model LL parameter test

close all
animal = 'Striatum';
alpha_step = 0.0014;
beta_step = 0.20;
gobias_step = 0.007;
nsamples = 500;
one7 = NLL*(-1) - log(7);
figure('Position',[0 0 400 1000])
colors = {'#0072BD','#D95319'};
nLL = nan(nsamples,3,2);
paraname = {'\alpha','\beta','gobias'};
for k = [2 4]
    for i = 1:3
        switch i
            case 1
                alpha_samples = nanmean(Para_opt(1,:,1,k*0.5)).*ones(1,nsamples) - ...
                    (0.5*nsamples.*ones(1,nsamples) - (1:nsamples))*alpha_step;
                beta_samples = nanmean(Para_opt(2,:,1,k*0.5)).*ones(1,nsamples);
                gobias_samples = nanmean(Para_opt(3,:,1,k*0.5)).*ones(1,nsamples);
            case 2
                beta_samples = nanmean(Para_opt(2,:,1,k*0.5)).*ones(1,nsamples) - ...
                    (0.5*nsamples.*ones(1,nsamples) - (1:nsamples))*beta_step;
                alpha_samples = nanmean(Para_opt(1,:,1,k*0.5)).*ones(1,nsamples);
                gobias_samples = nanmean(Para_opt(3,:,1,k*0.5)).*ones(1,nsamples);
            case 3
                gobias_samples = nanmean(Para_opt(3,:,1,k*0.5)).*ones(1,nsamples) - ...
                    (0.5*nsamples.*ones(1,nsamples) - (1:nsamples))*gobias_step;
                alpha_samples = nanmean(Para_opt(1,:,1,k*0.5)).*ones(1,nsamples);
                beta_samples = nanmean(Para_opt(2,:,1,k*0.5)).*ones(1,nsamples);
        end
        params = [alpha_samples;beta_samples;gobias_samples];
        for j = 1:nsamples
            
            nLL(j,i,k*0.5) = SUM_GNG_loglikeli_action(params(:,j),states(kkk,:),actions(kkk,:),rewards(kkk,:), Policy,QInit);
        end
        LL = -1.*nLL;
        subplot(3,1,i)
        hold on
        ffig{i,k*0.5} = plot(params(i,:),-1.*nLL(:,i,k*0.5),'--','Color',colors{k*0.5},'LineWidth',5);
        ax = gca;
        ax.FontSize = 20;
        ax.FontName = 'Times New Roman';
        xlim([params(i,1) params(i,end)])
        ylim([2*one7 0])
        xlabel(paraname{i})
        ylabel('LL')
        title(['Log-likelihood vs ' paraname{i}])
        hold on
        
        [~,ciindex] = mink(abs(LL(:,i,k*0.5) - one7*ones(size(LL(:,i,k*0.5),1),1)),2);
        if LL(ciindex(1),i,k*0.5)>0.95*one7 || LL(ciindex(1),i,k*0.5)<1.1*one7
            ciindex(1) = [];
        elseif LL(ciindex(2),i,k*0.5)>0.95*one7 || LL(ciindex(2),i,k*0.5)<1.1*one7
            ciindex(2) = [];
        end
        hold on
        cifig{i} = plot(params(i,ciindex),LL(ciindex,i,k*0.5),'-|','Color',colors{k*0.5},'MarkerSize',15,'LineWidth',3.5);
        

    end
end

subplot(3,1,3)
legend([ffig{3,1},ffig{3,2},cifig{3}],{'Stim','Unstim','95% CI'},'Location','southeast')
%% Trial Learned and Overall Performance

close all
pv = [];
animal = 'Prelimbic';

data(isnan([data.nogodelay])) = [];
data = data([data.stim] ==1);
fontsize = 40;

x = [zeros(length(data),1) ones(length(data),1)]';
tlearned = nan(2,length(data));
performance = nan(2,length(data));
finalperf = nan(2,length(data));
tlearned_to_90 = nan(2,length(data));
stdev_rxn = nan(2,length(data));
rxn_avg = nan(2,length(data));

for i = 1:length(data)
    for j = 1:2
    tlearned(j,i) = data(i).odor(j*2).tlearned;
    performance(j,i) = data(i).odor(j*2).performance;
    finalperf(j,i) = data(i).odor(j*2).finalperf;
    tlearned_to_90(j,i) = data(i).odor(j*2).tlearned_to_90;
    stdev_rxn(j,i) = data(i).odor(j*2).stdev_rxn;
    rxn_avg(j,i) = data(i).odor(j*2).rxn_avg;
    end
end

x1 = tlearned(1,:); x1(x1==inf)=nan;
y = tlearned(2,:);y(y==inf)=nan;
[h,p] = ttest2(x1,y,'Tail','left');
pv(1,1) = p;


figure('Name',num2str(data(1).animal),'Position',[250,250,1000,800])
clf
% subplot(1,4,1)
subplot(1,2,1)

plot(x,tlearned,'.--','MarkerSize',50,'LineWidth',5)
ax = gca;
ax.FontSize = fontsize;
ax.FontName = 'Times New Roman';
xlim([-.75 1.75])
xticks([0 1])
xticklabels({'Unstimulated','Stimulated'})
ylim('auto')
% title(sprintf('D%.f \n Time till learned: distribution',data(1).animal))
title('Trial learned')
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',14)
ylabel('Trial Learned')
% xlabel('No-Stim (0) vs Stim (1)','FontSize',14)
hold on
% text(0.05,0.95,[' p = ' num2str(p,3)],'Units','normalized','FontSize',18,'Color','black','Interpreter','latex')
tlearned = tlearned';
f = plot([0 1],[mean(tlearned(:,1).*[tlearned(:,1)~=inf],'omitnan'),mean(tlearned(:,2).*[tlearned(:,2)~=inf],'omitnan')],...
    '+','Color','black','MarkerSize',20,'LineWidth',5);
legend(f,'Mean')

x1 = performance(1,:); x1(x1==inf)=nan;
y = performance(2,:);y(y==inf)=nan;
[h,p] = ttest2(x1,y,'Tail','right');
pv(1,2) = p;

% subplot(1,4,2)
subplot(1,2,2)
plot(x,performance,'.--','MarkerSize',50,'LineWidth',5)
ax = gca;
ax.FontSize = fontsize;
ax.FontName = 'Times New Roman';
xlim([-0.75 1.75])
xticks([0 1])
xticklabels({'Unstimulated','Stimulated'})
ylim('auto')
% title(sprintf('D%.f \n Overall Performance: distribution',data(1).animal))
% title([ animal ': Overall Performance'])
title('Overall Performance')
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',14)
ylabel('Overall Performance')
% xlabel('No-Stim (0) vs Stim (1)','FontSize',14)
% text(0.05,0.95,[' p = ' num2str(p,3)],'Units','normalized','FontSize',18,'Color','black','Interpreter','latex')
performance = performance';
hold on
f= plot([0 1],[mean(performance(:,1).*[performance(:,1)~=inf],'omitnan'),mean(performance(:,2).*[performance(:,2)~=inf],'omitnan')],...
    '+','Color','black','MarkerSize',20,'LineWidth',5);
legend(f,'Mean')
% text([0.1 1.1],[mean(performance(:,1).*[performance(:,1)~=inf],'omitnan'),mean(performance(:,2).*[performance(:,2)~=inf],'omitnan')],...
%     {['$ \bar{x} = $' num2str(mean(performance(:,1).*[performance(:,1)~=inf],'omitnan'),3)],['$\bar{x} = $' num2str(mean(performance(:,2).*[performance(:,2)~=inf],'omitnan'),3)]},...
%     'VerticalAlignment','middle','HorizontalAlignment','left','Interpreter','latex','FontSize',16)

metrics = {' ','Trial Learned (unstim)','Trial Learned (stim)','Overall Performance (unstim)','Overall Performance (stim)'};
Mean = round([mean(tlearned(:,1).*[tlearned(:,1)~=inf],'omitnan'),mean(tlearned(:,2).*[tlearned(:,2)~=inf],'omitnan'),...
    mean(performance(:,1).*[performance(:,1)~=inf],'omitnan'),mean(performance(:,2).*[performance(:,2)~=inf],'omitnan')],3,'significant');

pv = round([pv,0,0],3,'significant');

t = num2cell([Mean;pv]);
t(:,2:5) = t;
t(:,1) = {'Mean';'p-value'};
t = array2table(t,'VariableNames',metrics);
table2latex(t,'test')
%% Sigmoid fitting CI: binary avg -> movavg

% requires you to run "align binary curves" section first
animal = 'Striatum';
close all
clear xydata
clear x1
clear cimean
clear xmean
clear jacobian
clear residual
clear x
clear ci
clear xmean
xydata{1} = [];
xydata{2} = [];
x1 = -299:300;
x = nan(1,3);
cimean = {};
xmean = {};
jacobian = {};
residual = {};
% titles = {'No-Go Unstimulated','No-Go Stimulated'};
plotcurve = figure('Position',[0 0 1000 1000]);
% plothist = figure('Position',[0 0 1000 800]);
% plotdistr = figure('Position',[0 0 1000 600]);
probcurve = mean(curve,1,'omitnan');
colors = {'#0072BD','#D95319'};
linewidth = 8;
for j = [2 4]
    

        current_curve = movmean(probcurve(1,:,j),10);
        xydata{j*0.5} = [xydata{j*0.5} [x1(~isnan(current_curve)); current_curve(~isnan(current_curve))]];

    fun = @(x,x1) x(3)./(1+exp(-x(1)*(x1-x(2))));
    if j ==2
        lb = [0.09,-10,0];
    else
        lb = [0,-10,0];
    end
    for i = 1:100
        
        x0 = randi(20,1,3);
        
        [x(i,:,j*0.5),resnorm,residual{j*0.5},exitflag,output,lambda,jacobian{j*0.5}] = lsqcurvefit(fun,x0,xydata{j*0.5}(1,:),xydata{j*0.5}(2,:),[0 -100 0],[10,200,1],...
            optimoptions('lsqcurvefit','Display','off'));
        
        ci(:,:,i,j*0.5) = nlparci(x(i,:,j*0.5),residual{j*0.5},'jacobian',jacobian{j*0.5});
    end
    
    cimean{j*0.5} = mean(ci(:,:,:,j*0.5),3,'omitnan');
    xmean{j*0.5} = mean(x(:,:,j*0.5),1);
    figure(plotcurve.Number)
    %     subplot(2,1,j*0.5)
    hold on
    plot(xydata{j*0.5}(1,:)',movmean(xydata{j*0.5}(2,:),10)','--','Color',colors{j*0.5},'LineWidth',linewidth*(2/3))
    hold on
    plot(x1,fun(xmean{j*0.5}(:,:),x1),'LineWidth',linewidth,'Color',colors{j*0.5})
    ax = gca;
    ax.FontSize = 40;
    ax.FontName = 'Times New Roman';
    xlim([-100 100])
    xticks(-100:20:100)
    ylabel('Performance')
    xlabel('Trial (relative to 3-Trial criteria)')
    title(['Sigmoid Fitting'])
%     text(0.05,0.85,sprintf('Alpha = %.4f CI (95): [%.3f, %.3f] \nTheta = %.3f CI (95): [%.3f, %.3f] \nEndPerf = %.3f CI (95): [%.3f, %.3f]',...
%         xmean{j*0.5}(1),cimean{j*0.5}(1,1),cimean{j*0.5}(1,2),xmean{j*0.5}(2),...
%         cimean{j*0.5}(2,1),cimean{j*0.5}(2,2),xmean{j*0.5}(3),cimean{j*0.5}(3,1),cimean{j*0.5}(3,2)),'Units',...
%         'normalized','FontSize',16,'Interpreter','latex')
%     text(0.05,0.65,'Where $$ f(x) = \frac{EndPerf}{1+e^{\alpha (x-\theta)}}$$','Units','normalized','Interpreter','latex',...
%         'FontSize',16)
   
    legend({'Unstim Data','Unstim Sigmoid Fitting','Stim Data','Stim Sigmoid Fitting'},'Location','southeast')
    ax.XTickLabel{ismember(ax.XTickLabel,'0')} = '3-Trial Point';
    figure(gcf)

end
%% Sigmoid fitting: movmean rmse vs movmean window
% requires you to run "align binary curves" section first
animal = 'Striatum';
close all
clear xydata
clear x1
clear cimean
clear xmean
clear jacobian
clear residual
clear x
clear ci
clear xmean
xydata{1} = [];
xydata{2} = [];
x1 = -299:300;
x = nan(1,3);
cimean = {};
xmean = {};
jacobian = {};
residual = {};
titles = {'No-Go Unstimulated','No-Go Stimulated'};
errfig = figure('Position',[0 0 500 300]);
paramfig = figure('Position',[0 0 500 900]);
probcurve = mean(curve,1,'omitnan');
maxwin = 50;
nrepeat = 10;
residualsum = nan(nrepeat,maxwin,2);
colors = {'blue','red'};
params = {'\alpha','\theta','EndPerf'};
marker = '.';
markersize = 15;
ylims = {[0,20],[-20,30],[0.5,1.1]};
for j = [2 4]
    if j ==2 && data(1).animal == 62
        lb = [0.09,-10,0];
    else
        lb = [0,-10,0];
    end
    for k = 1:maxwin
        xydata{j*0.5} = [];
        nonavg_curve = probcurve(1,:,j);
        current_curve = movmean(probcurve(1,:,j),k);
        xydata{j*0.5} = [x1(~isnan(current_curve)); current_curve(~isnan(current_curve))];
        
        fun = @(x,x1) x(3)./(1+exp(-x(1)*(x1-x(2))));
        for i = 1:nrepeat
            
            x0 = randi(20,1,3);
%             x0 = randn(1,3);
            [x(i,:,k,j*0.5),resnorm,residual,exitflag,output,lambda,jacobian{j*0.5}] = lsqcurvefit(fun,x0,xydata{j*0.5}(1,:),xydata{j*0.5}(2,:),lb,[10,200,1],...
                optimoptions('lsqcurvefit','Display','off'));
            residualsum(i,k,j*0.5) = sqrt(sum((fun(x(i,:,k,j*0.5),xydata{j*0.5}(1,:))-nonavg_curve(~isnan(current_curve))).^2)/size(xydata{j*0.5},2));
            
        end
        
    end
end

wind = 1:maxwin;
wind = wind.*ones(nrepeat,maxwin,2);
figure(errfig.Number)
unstim_plot = plot(wind(:,:,1),residualsum(:,:,1),marker,'MarkerSize',markersize,'Color',colors{1});
hold on
stim_plot = plot(wind(:,:,1),residualsum(:,:,2),marker,'MarkerSize',markersize,'Color',colors{2});
set(gca,'YScale','log')
grid on
ylim([0.2,0.4])
xlabel('Moving Mean Window')
ylabel('Error')
title('RMSE vs. moving mean window')
legend([unstim_plot(1) stim_plot(1)],{'Unstimulated','Stimulated'})
ax = gca;
ax.FontSize = 20;
ax.FontName = 'Times New Roman';
figure(paramfig.Number)
for jj = 1:size(x,2)
    
    subplot(size(x,2),1,jj)
    unstim_plot = plot(wind(:,:,1),squeeze(x(:,jj,:,1)),marker,'MarkerSize',markersize,'Color',colors{1});
    hold on
    stim_plot = plot(wind(:,:,1),squeeze(x(:,jj,:,2)),marker,'MarkerSize',markersize,'Color',colors{2});
    ylabel(params{jj})
    ylim(ylims{jj})
    xlabel('Moving mean Window')
    title([params{jj} ' vs. moving mean window'])
    if jj == 1
        set(gca,'yscale','log')
    end
    legend([unstim_plot(1) stim_plot(1)],{'Unstimulated','Stimulated'})
    ax = gca;
    ax.FontSize = 20;
    ax.FontName = 'Times New Roman';
%     set(gca,'yscale','log')
    grid on
    
end
%% Movmean different window sizes
close all
animal = 'Prelimbic';
plotcurve = figure('Position',[0 0 550 800]);
probcurve = mean(curve,1,'omitnan');
Legend = {};
windows = [1 5 10 20 40];
titles = {'No-Go Unstimulated','No-Go Stimulated'};

for i = windows
    for j = [2 4]
        if i ==1
            subplot(2,1,j*0.5)
            hold on
            plot(movmean(probcurve(:,:,j),i)','o','LineWidth',1.5,'Color',0.25.*ones(1,3))
            Legend{windows == i} = ['No Moving Mean'];

        else
            subplot(2,1,j*0.5)
            hold on
            plot(movmean(probcurve(:,:,j),i)','-','LineWidth',4.5)
            Legend{windows == i} = ['Window = ' num2str(i)];

        end
        
        if i == windows(end)
            ax = gca;
            ax.FontSize = 25;
            ax.FontName = 'Times New Roman';
            xlabel('Trial')
            ylabel('Performance')
            title(['Moving Mean, ' titles{j*0.5}],'FontName','Times New Roman')

        end
    end
end
legend(Legend,'Location','northwest','FontSize',20)