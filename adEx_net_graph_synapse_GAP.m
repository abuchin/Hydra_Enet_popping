
%% INTEGRATION PARAMETERS
clear;       % remove previous varibles
close all;

T=1001;        % total time, mS
dt=0.1;        % time step, ms

%time=(1:1:(T/dt)).*dt;      % TIME VECTOR
Df=1/dt;            % Delta function approximation

Tframe=0;    % movie
dTframe=5; % ms
frame=5;
%%

%% NETWORK PARAMETERS
N=100;        % Number of neurons
N_index=1:1:N;
%%

%% CONNECTIVITY
%
% connection strength
gEE_mean=400;  % mA/cm^2 (current synapses)  # 33.35
gEE_sigma=0;

% synaptic timescale
tau_EE=5.4;      % ms, AMPA 5.4

load('Adjacency_20X5.mat');

S_EE=random('Normal',gEE_mean,gEE_sigma,N,N);
S_EE=A_tube.*S_EE;    % connectivity matrix

% GAP junctions, adds them to all neurons
g_GAP=0;      % 0.0001

% trick to efficiently solve the border conditions
l=zeros(N,4);
for i=1:1:N 
   l(i,1:length(find(S_EE(i,1:end)>0)))=find(S_EE(i,1:end)>0);
end

% number of connections per element
m=l>0;
m=sum(m,2);

load('Sensors.mat');
sensory_cells=S; % sensory cells all over the body

load('Actuators_mouth_20.mat')
actuator_cells=R;

% excitable cells are all the rest
excitable_cells=setdiff(N_index,actuator_cells');
excitable_cells=setdiff(excitable_cells,sensory_cells');

% Representative cell
repr=round(N/2);             % number of representative neuron

%%


%%  Graph plot details

A_graph=graph(A_tube);
% tube parameters
m=5;    % wideness
n=20;   % length

[Y,X]=meshgrid(1:n,1:m);        
x=reshape(X,[1,m*n]);
y=reshape(Y,[1,m*n]);
%%


%% NEURON PARAMETERS
% cortical neurons, Brette-Gerstner 2005
c(1:N)=281; 
gl(1:N)=30;
%random('Normal',30,5,1,N); % heterogenous population, different leaks
el(1:N)=-70.6; 
vt(1:N)=-50.4;
delta(1:N)=2; 
vreset(1:N)=-70.6;
a(1:N)=4;
tauw(1:N)=144;
b(1:N)=80;
% spike mark and number of spikes
vspike=10;
%%

%% VOLUME MODEL
tau_volume=200;  % ms, volume relaxation
volume_0=10;            % 1, basic volume
thr_volume=2;    % threshould for the volume
a_volume=80;    % maximal value for threshold
%%

%% ACTUATORS MODEL
tau_A=100;                               % ms, actuators time delay ()
w_A=zeros(N,1);
w_A(actuator_cells)=1;                  % actuators weights
%%

%% SENSORS MODEL
w_S=zeros(N,1);
w_S(sensory_cells)=2000;              % weights of cells feeling the volume
%%

%% STIMULATION
IN(1:N)=600;                 % mA/cm^2 basic input to all cells
%%

%% ICs
% Neural net
V(1:N)=el;   % mV
W(1:N)=0;
V_repr=zeros(1); % mV
W_repr=zeros(1); % pA
% actuators
rate_A=zeros(1);
% volume
volume=zeros(1);

% SYNAPTIC VARIABLES
I_EE=zeros(N,1);       % reccurent syn input
I_GAP=zeros(N,1);      % reccurent GAP input
I_sensory=zeros(N,1);  % sensory input

% FIRINGS
firings=[];            % spike timings
fired=[];              % indexes of spikes fired at t
fired_ones=zeros(N,1); % vector of fired elements expressed in 1s and 0s
V_sp=50*ones(N,1);     % vector of times of elements that did not spike

fired_delta(1:N)=0;   % vector of fired delta-function
%%



%% TIME INTEGRATION LOOP
figure('units','normalized','outerposition',[0 0 1 1]); % show figure window

for t=1:1:round(T/dt)


% Synaptic variable
I_EE=(-I_EE/tau_EE + Df.*sum(S_EE(:,fired),2) )*dt + I_EE;
% GAP junciton currents
I_GAP=g_GAP*(m.*V'-sum(repmat(V',1,N).*S_EE,1)');
% Sensory input
I_sensory=volume(t)*w_S;
% Voltage
V=(dt./c).*(-gl.*(V-el)+gl.*delta.*exp((V-vt)./delta)-W +IN +I_EE'+I_GAP'+I_sensory') + V;
% A rate
rate_A(t+1)=(dt/tau_A)*(-rate_A(t) +Df*sum(fired_ones.*w_A)) +rate_A(t);
% volume
volume(t+1)=(dt/tau_volume)*(-volume(t) +volume_0 -volume_step(rate_A(t),a_volume,thr_volume) ) +volume(t);

if volume(t+1)<0    % make sure the volume is not negative
    volume(t+1)=0;
end
%



% FIRINGS Proccessing, avoid W jumps!
fired=[];
fired=find(V>=vspike);

% update fired ones vector
fired_ones=zeros(N,1);
fired_ones(fired)=1;

% reset condition       
if isempty(fired)==0
    V(fired)=vreset(fired);
    W(fired)=W(fired)+b(fired);
else
    W=(dt./tauw).*(a.*(V-el)-W) + W;
end

% Representative cell
V_repr(t)=V(repr);
W_repr(t)=W(repr);
I_repr_ext(t)=IN(repr);
I_repr_EE(t)=I_EE(repr);
% add stick to the repr neuron
if V_repr(t)>=-40
    V_repr(t-1)=10;
    V_repr(t)=vreset(repr);
end

% Record the spikes
t_fired(1:length(fired))=t;

if isrow(fired)==0
     fired=fired';
end

if isrow(t_fired)==0
     t_fired=t_fired';
end

spikes=horzcat(t_fired',fired');
firings=vertcat(firings,spikes);
t_fired=[];
spikes=[];
fired_delta(:)=0;           % vector of current delta functions of all neurons
fired_delta(fired)=Df;      

%%

%% FINAL PLOT / ONLINE PLOT

if (t-1)*dt==Tframe

% indexes of sensory cells - green
[N_index,sensory_temp]=ismember(firings(:,2),sensory_cells);
sensory_cells_index=find(sensory_temp);

% indexes of actuators - red
[N_index,actuators_temp]=ismember(firings(:,2),actuator_cells);
actuator_cells_index=find(actuators_temp);    
    
% indexes of excitable units - blue
[N_index,excitable_temp]=ismember(firings(:,2),excitable_cells);
excitable_cells_index=find(excitable_temp);    


subplot(2,2,1);
if isempty(firings)==0
%plot(firings(:,1)*dt,firings(:,2),'.','MarkerSize',5); % plot them all
plot(firings(actuator_cells_index,1)*dt,firings(actuator_cells_index,2),'.','MarkerSize',8,'Color','r'); hold on;
plot(firings(excitable_cells_index,1)*dt,firings(excitable_cells_index,2),'.','MarkerSize',8,'Color','b');
plot(firings(sensory_cells_index,1)*dt,firings(sensory_cells_index,2),'.','MarkerSize',8,'Color','g'); hold on;


end

ylabel('Cell index');
ylim([0 N]);
xlim([0 t*dt]);
set(gca,'FontSize',20);             % set the axis with big font
title(sprintf('adEx population, T=%d ms',round(Tframe)));
set(gca,'FontSize',20);             % set the axis with big font
box off;

subplot(2,2,2);
plot((1:t)*dt,V_repr);
ylabel('Voltage (mV)');
xlim([0 t*dt]);
set(gca,'FontSize',20);             % set the axis with big font
title('Representative neuron');
set(gca,'FontSize',20);             % set the axis with big font
box off;

subplot(2,2,3);
h=plot(A_graph,'Xdata',x,'Ydata',y);    %handle for the graph plot
h.MarkerSize=7;
% actuators
highlight(h,actuator_cells,'NodeColor','r')
% sensors
highlight(h,sensory_cells,'NodeColor','g')
title('Blue - excitable, Green - sensors, Red - actuators');
set(gca,'Fontsize',15);
box off;
%{
imagesc((reshape(V,sqrt(N),sqrt(N)))',[-75 0]); % Voltage distribution
set(gca,'Ydir','normal');
set(gca,'FontSize',20);             % set the axis with big font
title('Voltage distribution (mV)');
colorbar;
box off;
%}

subplot(2,2,4);
s=plot(A_graph,'Xdata',x,'Ydata',y);    %handle for the graph plot
set(gca, 'CLim', [-100, 10]);
s.NodeCData = V;
s.MarkerSize=7;
title('Voltage distribution');
set(gca,'Fontsize',15);
colormap('jet');
colorbar;
box off;

%{
imagesc((reshape((IN +I_EE'),sqrt(N),sqrt(N)))',[0 gEE_mean*4]); % Voltage distribution
set(gca,'Ydir','normal');
set(gca,'FontSize',20);             % set the axis with big font
colorbar;
title('Synaptic input (\muA/cm^2)');
box off;
%}

%%

MOV(frame)=getframe(gcf);

frame=frame+1;          % counter for the movie
Tframe=Tframe + dTframe;
%%

end

end

figure
plot((1:t+1)*dt,volume)
set(gca,'FontSize',20);             % set the axis with big font
xlabel('Time (ms)')
ylabel('Volume (units)')



% Save the movie
%movie2avi(MOV,'adEX_popping_net.avi','fps',10,'quality',1);

