% ======================================================== % 
% Files of the Matlab programs included in the book:       %
% Xin-She Yang, Nature-Inspired Metaheuristic Algorithms,  %
% Second Edition, Luniver Press, (2010).   www.luniver.com %
% ======================================================== %    

% -------------------------------------------------------- %
% Firefly Algorithm for constrained optimization using     %
% for the design of a spring (benchmark)                   % 
% by Xin-She Yang (Cambridge University) Copyright @2009   %
% -------------------------------------------------------- %
%fin ans  13APR 307.6102  399.4075  339.8255  690.6627  715.1937
%   15APR
% .param w_m1=342.913762u 
% .param w_m2=342.913762u 
% .param w_m3=489.179523u 
% .param w_m4=489.179523u 
% .param w_m5=241.866358u 
% .param w_m6=412.882248u 
% .param w_m7=804.686580u 
% .param w_m8=241.866358u


function fa_ndim
% parameters [n N_iteration alpha betamin gamma]
para=[20 5 0.5 0.2 1];

% Simple bounds/limits for d-dimensional problems
d=5;
%Lb=[209.7 1.212 4.194 38.39 13.28 ];
Lb=[244.74 135.43 133.16 4.89 640];
%Ub=[419.5804 606 419.4 800 1000.67 ];
Ub=[550.54 145.43 340.16 489.51 658];


% Initial random guess
u0=Lb+(Ub-Lb).*rand(1,d);



[u,fval,NumEval]=ffa_mincon(u0,Lb,Ub,para);

% Display results
bestsolution=u
bestojb=fval
total_number_of_function_evaluations=NumEval

%cost function
function [z]=cost(x)
w_m1 = x(1);

w_m2 = x(1);

w_m3 = x(2);

w_m4 = x(2);

w_m5 = x(3)+120;

w_m6 = x(4);

w_m7 = x(5);

w_m8 = x(3);


f1 = fopen('width.txt','w');

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

% Execute ELDO to do AC Simulation 
eldo_run=1;
while eldo_run==1

    eldo_run= system('eldo ./temp.cir');
end

disp(' ELDO chalgaya !!!')

% Create outputs.txt file

disp(' vmax and f3db.txt :CREATION COMPLETE!!')
system('./extractdata');%SCRIPT

%extraction of required data
disp('outputs.txt nikal rahe hain!!!')
disp('Vmax')
width_outputfile = fopen('vmax.txt','r');
parameters1=fscanf(width_outputfile,'%f\n');
%parameters(2)=fscanf(width_outputfile,'%f\n');
fclose(width_outputfile);
clear width_outputfile; 

% disp('F3DB')
% width_outputfile = fopen('f3db.txt','r');
% parameters2=fscanf(width_outputfile,'%f\n');
% %parameters(2)=fscanf(width_outputfile,'%f\n');
% fclose(width_outputfile);
% clear width_outputfile;


disp('param1');
disp(parameters1);
%disp(parameters2);
temp1 = parameters1;     %*parameters2;


%extracting dimensions
%disp('extracting dimensions');
%outputfile = fopen('chk1.txt','r');
%parameter1=fscanf(outputfile,'%f\n');
%fclose(outputfile);
%clear outputfile; 

%outputfile = fopen('chk2.txt','r');
%parameter2=fscanf(outputfile,'%f\n');
%fclose(outputfile);
%clear outputfile; 

%disp('param2');
%disp(parameter1);
%disp(parameter2);
%temp2=parameter1*parameter2;

%z=temp1*temp2;
z=temp1;
disp(z);
disp('END!!!');



%disp('Extracting from vmax and f3db !!!')
%width_outputfile = fopen('./vmax.txt','r');
%parameters(1)=fscanf(width_outputfile,'%f\n');
%fclose(width_outputfile);
%clear width_outputfile; 

%width_outputfile = fopen('./f3db.txt','r');
%parameters(2)=fscanf(width_outputfile,'%f\n');
%parameters(3)=fscanf(width_outputfile,'%f\n');
%fclose(width_outputfile);
%clear width_outputfile; 
%nargin

%disp(length(parameters));
%disp('param length\n');
%disp(parameters);
%z = parameters(1)*parameters(2)*parameters(3);







%%% --------------------------------------------------%%%
%%% Do not modify the following codes unless you want %%%
%%% to improve its performance etc                    %%%
% -------------------------------------------------------
% ===Start of the Firefly Algorithm Implementation ======
%         Lb = lower bounds/limits
%         Ub = upper bounds/limits
%   para == optional (to control the Firefly algorithm)
% Outputs: nbest   = the best solution found so far
%          fbest   = the best objective value
%      NumEval = number of evaluations: n*MaxGeneration
% Optional:
% The alpha can be reduced (as to reduce the randomness)
% ---------------------------------------------------------

% Start FA
function [nbest,fbest,NumEval]...
           =ffa_mincon(u0, Lb, Ub, para)
% Check input parameters (otherwise set as default values)
if nargin<4, para=[1 1 0.25 0.20 1]; end
if nargin<3, Ub=[]; end
if nargin<2, Lb=[]; end
if nargin<1
disp('Usuage: FA_mincon(@cost,u0,Lb,Ub,para)');
end

% n=number of fireflies
% MaxGeneration=number of pseudo time steps
% ------------------------------------------------
% alpha=0.25;      % Randomness 0--1 (highly random)
% betamn=0.20;     % minimum value of beta
% gamma=1;         % Absorption coefficient
% ------------------------------------------------
n=para(1);  MaxGeneration=para(2);
alpha=para(3); betamin=para(4); gamma=para(5);

% Total number of function evaluations
NumEval=n*MaxGeneration;

% Check if the upper bound & lower bound are the same size
if length(Lb) ~=length(Ub)
    disp('Simple bounds/limits are improper!');
    return
end

% Calcualte dimension
d=length(u0);

% Initial values of an array
zn=ones(n,1)*10^100;
% ------------------------------------------------
% generating the initial locations of n fireflies
[ns,Lightn]=init_ffa(n,d,Lb,Ub,u0);

% Iterations or pseudo time marching
for k=1:MaxGeneration     %%%%% start iterations

% This line of reducing alpha is optional
 %alpha=alpha_new(alpha,MaxGeneration);

% Evaluate new solutions (for all n fireflies)
for i=1:n
    z=cost(ns(i,:));
   zn(i)=cost(ns(i,:)); %fhandle is pointing towards the cost function 
   Lightn(i)=zn(i);
end

% Ranking fireflies by their light intensity/objectives

if( z > 60)
[Lightn,Index]=sort(zn);
ns_tmp=ns;
for i=1:n
 ns(i,:)=ns_tmp(Index(i),:);
end
end

%% Find the current best
nso=ns; Lighto=Lightn;
nbest=ns(1,:); Lightbest=Lightn(1);

% For output only
fbest=Lightbest;

% Move all fireflies to the better locations
[ns]=ffa_move(n,d,ns,Lightn,nso,Lighto,nbest,...
      Lightbest,alpha,betamin,gamma,Lb,Ub);

end   %%%%% end of iterations

% -------------------------------------------------------
% ----- All the subfunctions are listed here ------------
% The initial locations of n fireflies
function [ns,Lightn]=init_ffa(n,d,Lb,Ub,u0)
  % if there are bounds/limits,
if length(Lb)>0
   for i=1:n
   ns(i,:)=Lb+(Ub-Lb).*rand(1,d);
   end
else
   % generate solutions around the random guess
   for i=1:n
   ns(i,:)=u0+randn(1,d);
   end
end

% initial value before function evaluations
Lightn=ones(n,1)*10^100;

% Move all fireflies toward brighter ones
function [ns]=ffa_move(n,d,ns,Lightn,nso,Lighto,...
             nbest,Lightbest,alpha,betamin,gamma,Lb,Ub)
% Scaling of the system
scale=abs(Ub-Lb);

% Updating fireflies
for i=1:n
% The attractiveness parameter beta=exp(-gamma*r)
   for j=1:n
      r=sqrt(sum((ns(i,:)-ns(j,:)).^2));
      % Update moves
      if Lightn(i)<Lighto(j) % Brighter and more attractive-inequality change as we want minima
        beta0=1; beta=(beta0-betamin)*exp(-gamma*r.^2)+betamin;
        tmpf=alpha.*(rand(1,d)-0.5).*scale;
         ns(i,:)=ns(i,:).*(1-beta)+nso(j,:).*beta+tmpf;
       end
   end % end for j

end % end for i

% Check if the updated solutions/locations are within limits
[ns]=findlimits(n,ns,Lb,Ub);

% This function is optional, as it is not in the original FA
% The idea to reduce randomness is to increase the convergence,
% however, if you reduce randomness too quickly, then premature
% convergence can occur. So use with care.
%function alpha=alpha_new(alpha,NGen)
% alpha_n=alpha_0(1-delta)^NGen=10^(-4);
% alpha_0=0.9
%delta=1-(10^(-4)/0.9)^(1/NGen);
%alpha=(1-delta)*alpha;

% Make sure the fireflies are within the bounds/limits
function [ns]=findlimits(n,ns,Lb,Ub)
for i=1:n
     % Apply the lower bound
  ns_tmp=ns(i,:);
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);

  % Apply the upper bounds
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move
  ns(i,:)=ns_tmp;
end
%paste coeffs to cir file later
