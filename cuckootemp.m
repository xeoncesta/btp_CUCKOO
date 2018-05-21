%OPTIMISATION OF ANALOG CIRCUITS USING CUCKOO SEARCH ALGORITHM

%DEVELOPED BY :
%AMITESH
%ABHISHEK
%GARV
%KARAN
%NSIT,DWARAKA SEC-3

% =============================================================== %
% Notes:                                                          %
% Different implementations may lead to slightly different        %
% behavour and/or results, but there is nothing wrong with it,    %
% as this is the nature of random walks and all metaheuristics.   %
% -----------------------------------------------------------------

function [bestnest,fmin]=cuckoo_search(n)
if nargin<1
% Number of nests (or different solutions)
n=5;
end

% Discovery rate of alien eggs/solutions
pa=0.25;
f=0*ones(1,1);
%% Change this if you want to get better results
% Tolerance
Tol=1.0e-8;
%% Simple bounds of the search domain
% Lower bounds
%nd=15; 
%Lb=-5*ones(1,nd); 
% Upper bounds
%Ub=5*ones(1,nd);
Lb=[244.74 135.43 133.16 4.89 640]; %  --for our ckt
Ub=[550.54 143.43 340.16 489.51 658];  % --for our ckt
%z=0*ones(n,1);
arr=0*ones(250,1);
temp=[1 1 1];
low_lim = [55 40 45 ];	%0.002 10 10^7 10 ];
up_lim = [80 50 60 ];	%0.0005 20^6 10000];
% Random initial solutions
for i=1:n
nest(i,:)=Lb+(Ub-Lb).*rand(size(Lb));
end


% Get the current best
fitness=55*ones(n,1);
[fmin,bestnest,nest,fitness]=get_best_nest(nest,nest,fitness);
%temp
N_iter=0;
%% Starting iterations
while ((fmin>Tol) && (N_iter<100))

    % Generate new solutions (but keep the current best)
     new_nest=get_cuckoos(nest,bestnest,Lb,Ub);   
     [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness);
    % Update the counter
      N_iter=N_iter+n; 
    % Discovery and randomization
      new_nest=empty_nests(nest,Lb,Ub,pa) ;
    
    % Evaluate this set of solutions
      [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness);
    % Update the counter again
      N_iter=N_iter+n;
    % Find the best objective so far  
    
	%if fnew>fmin,
        %  bestnest=best;
        %end
%	disp(strcat('gain',temp(1,1)));
 	
width_outputfile = fopen('./pm.txt','r');
parameters(3)=fscanf(width_outputfile,'%f\n');
fclose(width_outputfile);
clear width_outputfile; 
parameters(3) = parameters(3) +180;
temp(1,3)=parameters(3);

			if(fnew<low_lim(1))
                fmin = fnew*100;
            elseif (fnew>low_lim(1) && fnew<up_lim(1))
                fmin = fnew*10;
            elseif (fnew>up_lim(1) && temp(1,3)>40)
                fmin = fnew;
				bestnest=best;
             end
            
	
	
end %% End of iterations

%% Post-optimization processing
%% Display all the nests
disp(strcat('Total number of iterations=',num2str(N_iter)));
%disp(strcat('pm=',num2str(temp(1,3))));
fmin
bestnest
%grid on

 %plot(bestnest);
 title('WIDTH OF MOSFETS');
 xlabel('No of mos');
 ylabel('Width in um');

% surf(Lb,Ub,nest)
% xlabel('Lower lim');
% ylabel('Upper lim');
% zlabel('Nest for iteration');


%% --------------- All subfunctions are list below ------------------
%% Get cuckoos by ramdom walk
function nest=get_cuckoos(nest,best,Lb,Ub)
% Levy flights
n=size(nest,1);
% Levy exponent and coefficient
% For details, see equation (2.21), Page 16 (chapter 2) of the book
% X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

for j=1:n
    s=nest(j,:);
    % This is a simple way of implementing Levy flights
    % For standard random walks, use step=1;
    %% Levy flights by Mantegna's algorithm
    u=randn(size(s))*sigma;
    v=randn(size(s));
    step=u./abs(v).^(1/beta);
  
    % In the next equation, the difference factor (s-best) means that 
    % when the solution is the best solution, it remains unchanged.     
    stepsize=0.01*step.*(s-best);
    % Here the factor 0.01 comes from the fact that L/100 should the typical
    % step size of walks/flights where L is the typical lenghtscale; 
    % otherwise, Levy flights may become too aggresive/efficient, 
    % which makes new solutions (even) jump out side of the design domain 
    % (and thus wasting evaluations).
    % Now the actual random walks or flights
    s=s+stepsize.*randn(size(s));
   % Apply simple bounds/limits
   nest(j,:)=simplebounds(s,Lb,Ub);
end

%% Find the current best nest
function [fmin,best,nest,fitness]=get_best_nest(nest,newnest,fitness)
% Evaluating all new solutions
for j=1:size(nest,1)
    fnew=fobj(newnest(j,:));
    if fnew>=fitness(j)
       fitness(j)=fnew;
       nest(j,:)=newnest(j,:);
    end
% plot(fitness,'r')
% hold on
end
% Find the current best
[fmin,K]=max(fitness) ;
best=nest(K,:); 
plot(nest,'r');
hold on;


%plot(best);

%% Replace some nests by constructing new solutions/nests
function new_nest=empty_nests(nest,Lb,Ub,pa)
% A fraction of worse nests are discovered with a probability pa
n=size(nest,1);
% Discovered or not -- a status vector
K=rand(size(nest))>pa;

% In the real world, if a cuckoo's egg is very similar to a host's eggs, then 
% this cuckoo's egg is less likely to be discovered, thus the fitness should 
% be related to the difference in solutions.  Therefore, it is a good idea 
% to do a random walk in a biased way with some random step sizes.  
%% New solution by biased/selective random walks
stepsize=rand*(nest(randperm(n),:)-nest(randperm(n),:));
new_nest=nest+stepsize.*K;
for j=1:size(new_nest,1)
    s=new_nest(j,:);
  new_nest(j,:)=simplebounds(s,Lb,Ub);  
end

% Application of simple constraints
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);
  
  % Apply the upper bounds 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move 
  s=ns_tmp;

%% You can replace the following by your own functions
% A d-dimensional objective function
%function z=fobj(u)
%% d-dimensional sphere function sum_j=1^d (u_j-1)^2. 
%  with a minimum at (1,1, ...., 1); 
%z=sum((u-1).^2);
function z=fobj(u)
disp('we are here in objective_func(u)');
%paste coeffs to cir file later
%first 6 entry are width and next 6 are length
		%If fails then try optimising the width only and fix length..
%w_m1m2 = u(1);
%w_m3 = u(2);
%w_m4 = u(3);
%w_m5 = u(4);
%w_m6m7 = u(5);
%w_m8 = u(6);
%l_m1m2 = u(7);
%l_m3 = u(8);
%l_m4 = u(9);
%l_m5 = u(10);
%l_m6m7 = u(11);
%l_m8 = u(12);
w_m1 = u(1);

w_m2 = u(1);

w_m3 = u(2);

w_m4 = u(2);

w_m5 = u(3)+120;

w_m6 = u(4);

w_m7 = u(5);

w_m8 = u(3);


f1 = fopen('width.txt','w');		%Replace by our width
% fprintf(f1,'.param w_m1m2=%fu \r\n',w_m1m2);
% fprintf(f1,'.param w_m3=%fu \r\n',w_m3);
% fprintf(f1,'.param w_m4=%fu \r\n',w_m4);
% fprintf(f1,'.param w_m5=%fu \r\n',w_m5);
% fprintf(f1,'.param w_m6m7=%fu \r\n',w_m6m7);
% fprintf(f1,'.param w_m8=%fu \r\n',w_m8);
% fprintf(f1,'.param l_m1m2=%fu \r\n',l_m1m2);
% fprintf(f1,'.param l_m3=%fu \r\n',l_m3);
% fprintf(f1,'.param l_m4=%fu \r\n',l_m4);
% fprintf(f1,'.param l_m5=%fu \r\n',l_m5);
% fprintf(f1,'.param l_m6m7=%fu \r\n',l_m6m7);
% fprintf(f1,'.param l_m8=%fu \r\n',l_m8);
% fclose(f1);
fprintf(f1,'.param w_m1=%fu \r\n',w_m1);

fprintf(f1,'.param w_m2=%fu \r\n',w_m2);

fprintf(f1,'.param w_m3=%fu \r\n',w_m3);

fprintf(f1,'.param w_m4=%fu \r\n',w_m4);

fprintf(f1,'.param w_m5=%fu \r\n',w_m5);

fprintf(f1,'.param w_m6=%fu \r\n',w_m6);

fprintf(f1,'.param w_m7=%fu \r\n',w_m7);

fprintf(f1,'.param w_m8=%fu \r\n',w_m8);

fclose(f1);
clear f1;

%Execute ELDO 
eldo_run=1;
while eldo_run==1
        eldo_run = system('eldo ./temp.cir');		%our cmos_ota files
end

disp(' ELDO Ran !!!')

% Create outputs.txt file
    disp('Create the outputs files !!')

 	%run our scripts
    system('./extractdata');

%extraction of required data
% vmax
disp('Extracting from outputs.txt !!!')
width_outputfile = fopen('./vmax.txt','r');
parameters(1)=fscanf(width_outputfile,'%f\n');
mult_vmax=fscanf(width_outputfile,'%s\n');
fclose(width_outputfile);
clear width_outputfile; 

if(mult_vmax == 'm')
    mult_vmax = 0.001;
elseif(mult_vmax == 'K')
    mult_vmax = 1000;
elseif(isempty(mult_vmax))
    mult_vmax = 1;
elseif(mult_vmax == 'MEG')
    mult_vmax = 1000000;
end
%f3db
width_outputfile = fopen('./f3db.txt','r');
parameters(2)=8;
parameters(2)=fscanf(width_outputfile,'%f\n');
if(isempty(parameters(2)))
    parameters(2)=8;
end
mult_f3db=fscanf(width_outputfile,'%s\n');
fclose(width_outputfile);
clear width_outputfile; 

if(mult_f3db == 'm')
    mult_f3db = 0.001;
elseif(mult_f3db == 'K')
    mult_f3db = 1000;
elseif(isempty(mult_f3db))
    mult_f3db = 1;
elseif(mult_f3db == 'MEG')
    mult_f3db = 1000000;
end

%pm
width_outputfile = fopen('./pm.txt','r');
parameters(3)=fscanf(width_outputfile,'%f\n');
% if(isempty(parameters(3)))
%     parameters(3)=47;
% end
fclose(width_outputfile);
clear width_outputfile; 
parameters(3) = parameters(3) +180;


%power
% width_outputfile = fopen('./pow.txt','r');
% parameters(3)=fscanf(width_outputfile,'%f\n');		%changed index
% mult_pow=fscanf(width_outputfile,'%s\n');
% fclose(width_outputfile);
% clear width_outputfile; 
% if(mult_pow == 'MWatt')
%     mult_pow = 0.001;
% elseif(mult_pow == 'KWatt')
%     mult_pow = 1000;
% elseif(mult_pow == 'UWatt')
%     mult_pow = 0.000001;   
% elseif(isempty(mult_pow))
%     mult_pow= 1;
% elseif(mult_pow == 'MEGWatt')
%     mult_pow = 1000000;
% end
%sr
%width_outputfile = fopen('./sr.txt','r');
%parameters(5)=fscanf(width_outputfile,'%f\n');
%mult_sr=fscanf(width_outputfile,'%s\n');
%fclose(width_outputfile);
%clear width_outputfile; 

%if(mult_sr == 'M')
%    mult_sr = 0.001;
%elseif(mult_sr == 'K')
%    mult_sr = 1000;
%elseif(mult_sr == 'U')
%    mult_sr = 0.000001;   
%elseif(isempty(mult_sr))
%    mult_sr= 1;
%elseif(mult_sr == 'MEG')
%    mult_sr = 1000000;
%end


%area
%parameter(6) = w_m1m2*l_m1m2*2 + w_m3*l_m3 + w_m4*l_m4 + w_m5*l_m5 + w_m6m7*l_m6m7 + w_m8*l_m8;
%IC


parameters(1) = parameters(1)*mult_vmax;
parameters(2) = parameters(2)*mult_f3db;

%parameters(4) = parameters(4)*mult_pow;
%parameters(5) = parameters(5)*mult_sr;
z = parameters(1); %*parameters(2);
x = parameters(2);
y= parameters(3);
disp('/n Curr obj ended');
disp('/n setting temp');
temp(1,1)=z;


%temp(1,3)=y;

%temp(1,2)=x;

%z
%x

%z
%%%%% ============ end ====================================


