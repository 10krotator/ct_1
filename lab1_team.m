clear all
close all
clc

s=tf('s');
G=minreal(3*(-s+1)/((5*s+1)*(10*s+1)));
% G=3*(1+s/(-1))/((s/(1/5)+1)*(s/(1/10)+1))
disp('Properties of G')
[Gm,Pm,wp,wc]=margin(G)

figure()
bode(G)

%% 4.1.1. Use the procedure introduced in the basic course to construct a lead-
% lag controller which eliminates the static control error for a step response 
% in the reference signal.

% Desired values
wc_des=0.4;
Pm_des=30;

Pmc=atand(-wc_des)-atand(wc_des/(1/5))-atand(wc_des/(1/10));
dPm_des=Pm_des+(-180-Pmc)+5.7

% Lead compensator
N=2;
beta=1/N;
b=wc_des/sqrt(N);
TD=1/b;
Flead=(TD*s+1)/(beta*TD*s+1);
% K=1/(abs(evalfr(Flead,1j*wc_des))*abs(evalfr(G,1j*wc_des)));

% Lag compensator, M high but not Inf
a=0.1*wc_des;
TI=1/a;
M=1e30;
gamma=1/M;
Flag=(TI*s+1)/(TI*s+gamma);
K=1/(abs(evalfr(Flead,1j*wc_des))*abs(evalfr(Flag,1j*wc_des))*...
    abs(evalfr(G,1j*wc_des)));

F=K*Flead*Flag;
G0=F*G;
Gc=G0/(1+G0);
figure()
bode(G0)

disp('Properties of G0 in 4.1.1.')
[Gm,Pm,wp,wc]=margin(G0)

%% 4.1.2. Determine the bandwidth of the closed-loop system and the resonance
% peak MT. Also, determine the rise time and the overshoot for step changes in the
% reference when the controller designed in 4.1.1. is used.
figure()
bode(Gc)
figure()
step(Gc)
disp('Properties of Gc')
wB=bandwidth(Gc)
wtest=linspace(0,1,1000);
[MT,idxT]=max(abs(freqresp(Gc,j*wtest)))
wT=wtest(idxT)

sinfoGc = stepinfo(Gc);
Tr=sinfoGc.RiseTime
M=sinfoGc.Overshoot     % [percent] 

%% 4.1.3. Modify the controller in 4.1.1. such that the phase margin increases
% to 50? while the cross-over frequency is unchanged. For the corresponding 
% closed-loop system, determine the bandwidth and resonance peak. Also, 
% determine the rise time and the overshoot of the step response.

% Desired values
wc_des=0.4;
Pm_des=50;

Pmc=atand(-wc_des)-atand(wc_des/(1/5))-atand(wc_des/(1/10));
dPm_des=Pm_des+(-180-Pmc)+5.7

% Lead compensator
N=4;
beta=1/N;
b=wc_des/sqrt(N);
TD=1/b;
Flead=(TD*s+1)/(beta*TD*s+1);
% K=1/(abs(evalfr(Flead,1j*wc_des))*abs(evalfr(G,1j*wc_des)));

% Lag compensator, M high but not Inf
a=0.1*wc_des;
TI=1/a;
M=1e30;
gamma=1/M;
Flag=(TI*s+1)/(TI*s+gamma);
K=1/(abs(evalfr(Flead,1j*wc_des))*abs(evalfr(Flag,1j*wc_des))*...
    abs(evalfr(G,1j*wc_des)));

F=K*Flead*Flag;
G0=F*G;
Gc=G0/(1+G0);
figure()
bode(G0)

disp('Properties of G0 in 4.1.3.')
[Gm,Pm,wp,wc]=margin(G0)

figure()
bode(Gc)
figure()
step(Gc)
disp('Properties of Gc in 4.1.3.')
wB=bandwidth(Gc)
wtest=linspace(0,1,1000);
[MT,idxT]=max(abs(freqresp(Gc,j*wtest)))
wT=wtest(idxT)

sinfoGc = stepinfo(Gc);
Tr=sinfoGc.RiseTime
M=sinfoGc.Overshoot     % [percent] 

%% 4.2.1. For which frequencies is control action needed? Control is needed
% at least at frequencies where |Gd(j?)| > 1 in order for disturbances to 
% be attenuated. Therefore the cross-over frequency must be large enough. 
% First, try to design Fy such that L(s) ? ?c/s and plot the closed-loop 
% transfer function from d to y and the corresponding step response. 
% (A simple way to find L = ?c/s is to let Fy = G?1?c/s. However, this 
% controller is not proper. A procedure to fix this is to ?add? a number of
% poles in the controller such that it becomes proper. How should these 
% poles be chosen?)

close all

G=minreal(20/((s+1)*((s/20)^2+s/20+1)));
Gd=minreal(10/(s+1));

% figure()
% bode(Gd)
% |Gd|=1 for w=sqrt(99)
% [Gm,Pm,wp,wc]=margin(Gd)

% Rule of thumb: Choose poles at (at least) 4 times the crossover frequency
% of the disturbance. Here: 8 times.
wc_des=2.5*sqrt(99);
Fy=minreal(1/(s/wc_des+1)^2*(G^(-1)*Gd));
% figure()
% bode(G*Fy)

% 4.2.2. Let us now reconstruct Fy according to the instructions above. We
% will start with the disturbance attenuation. In a second step, 
% adjustments can be made on Fr to obtain the desired reference tracking 
% properties. Start by choosing Fy according to (1). Try different 
% approximations of the product G?1Gd in the controller, and choose ?I 
% large enough so that step disturbances are attenuated according to the
% specifications.

wI=10;
Fy=minreal((s+wI)/s*Fy);
% figure()
% bode(G*Fy)

% 4.2.3. To fulfill the reference tracking specifications, we can combine 
% lead lag control and prefiltering of the reference signal. First, try to 
% add lead action to Fy to reduce the overshoot. Then it can be necessary 
% to add prefilter action to fulfill all specifications. Note that Fr 
% should be as simple as possible (why?). Also, remember to check the size 
% of the control signal (u = FyFrSr ? FyGdSd)! Typically a low pass filter 
% is chosen, for example
% Fr =
% 1
% 1 + ?s
% .
% Hint: Consider the signals r and d independently.

tau=0.1;
Fr=minreal(1/(1+tau*s));

% S=(1+Fy*G)^(-1);
% figure()
% bode(Fy*Fr*G*S)
% figure()
% bode(Fy*Gd*S)

% L=G*Fy;
% Gcr=Fr*L/(1+L);
% Gcd=Gd/(1+L);

% figure()
% bode(Fy*G)
% figure()
% bode(Gcr)
% figure()
% bode(Gcd)
% 
% figure()
% step(Gcr)
% figure()
% step(Gcd)
% stepinfo(Gcr)

% Lead compensator
N=4;
beta=1/N;
b=wc_des/sqrt(N);
TD=1/b;
Flead=minreal((TD*s+1)/(beta*TD*s+1));
K=1/(abs(evalfr(Flead,1j*wc_des))*abs(evalfr(Fy,1j*wc_des))*...
    abs(evalfr(G,1j*wc_des)));
% K=1;

Fy=minreal(K*Flead*Fy);
L=minreal(G*Fy);
Gcr=minreal(Fr*L/(1+L));
Gcd=minreal(Gd/(1+L));
figure()
bode(L)
figure()
bode(Gcr)
figure()
bode(Gcd)

figure()
step(Gcr)
figure()
step(Gcd)
stepinfo(Gcr)
% u = Fr*Fy*r - Fy*y
% y = L*Fr/(1+L)*r
% u = Fr*Fy/(1+L)*r
figure()
step(Fr*Fy/(1+L))

S=1/(1+L);
T=L/(1+L);

figure()
bodemag(S)
hold on; grid on;
bodemag(T)