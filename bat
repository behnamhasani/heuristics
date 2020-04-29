

%% Insert Data

 data=InsertData();%%% initialize cities , calculate Euclidean distance

%%  Parameters Setting

nvar=data.N;       % Number of Variables (cities)

lb=0*ones(1,nvar);    % Lower Bound of variables
ub=1*ones(1,nvar);    % Upper Bound of variables

npop=20;        % Population size
maxiter=200;   % Number of generations
A=0.5;          % Loudness  (constant or decreasing)
r=0.5;         % Pulse rate (constant or decreasing)
Qmin=0;        % Frequency minimum
Qmax=0.1;      % Frequency maximum

Alpha=0.99;    % Reduction factor A
Lambda=0.01;   % Reduction factor r


data.lb=lb;
data.ub=ub;
data.nvar=nvar;
%% Initialize the population/solutions

tic
emp.v=[]; %%%%%% velocity
emp.x=[]; %%%%%% position of bats
emp.fit=[]; %%%%% fitness

pop=repmat(emp,npop,1);%%create population of bats/


for i=1:npop
    pop(i).x = unifrnd(lb, ub);%%%%% initialize position of bat(i)
    pop(i).v = 0; %%%% at first all bat have 0 velocity
    pop(i) = fitness(pop(i),data); %%% calculate fitness of bat(i)
end


% Find the initial gpop(X* or best) solution
[~,ind]=min([pop.fit]);%%%% finding index of bat with lowest fitness
gpop=pop(ind);


%% Main Loop
BEST=zeros(maxiter,1);
MEAN=zeros(maxiter,1);
A=A*ones(npop,1);
r=r*ones(npop,1);
r0=r;

for iter=1:maxiter
    
    % Loop over all bats/solutions
    
    for i=1:npop
        
        Q=unifrnd(Qmin,Qmax); % Random Frequency
        
        newpop.v=pop(i).v+(pop(i).x-gpop.x)*Q; % Update Velocity
        
        newpop.x=pop(i).x+newpop.v;            % Update Position
        
        
        
        % Pulse rate
        if rand>r(i)
            newpop.x=gpop.x+0.001*randn(1,nvar);  % The factor 0.001 limits the step sizes of random walks
        end
        
        
        newpop= fitness(newpop,data);
        
        
        % Update if the solution improves, or not too loud
        if (newpop.fit<=pop(i).fit) && (rand<A(i))
            pop(i)=newpop; 
            r(i)=r0(i)*(1-exp(-Lambda*iter));
            A(i)=A(i)*Alpha;
        end
        
        % Update the  gpop solution/X* (best solution)
        if (newpop.fit<=gpop.fit)
            gpop=newpop;
        end
  
    end
    
    

    
    
    
    BEST(iter)=gpop.fit;
    MEAN(iter)=mean([pop.fit]);

    disp(['Iter ' num2str(iter) ' BEST = ' num2str(BEST(iter))]);  %%display iteration and X*  
    
   % Plot Best Solution
  PlotBestSol(gpop,data,iter)    
end

%% Results

disp(' ')
x=gpop.x;
[~,x]=sort(x);
x=[x x(1)];
disp([ ' Best x = '  num2str(x)])
disp([ ' BEST fitness = '  num2str(gpop.fit)]);
disp([ ' Time = '  num2str(toc)]);

figure
semilogy(BEST,'r');
hold on
semilogy(MEAN,'b');
xlabel('Iteration');
ylabel('Fitness');
legend('BEST','MEAN');
title('Bat algorithm ')


