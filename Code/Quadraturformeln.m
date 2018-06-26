%%This programm compares uniform, trapez and Simpson-Sceme for %%
%%approximation of solutions to the Fredholm integral equation.%%
%%written by Tim Jaschek as a part of his bachelor thesis%%

%%Used to generate data for Tabular 6.1 and 6.2 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Import the class Kernels which contains some Kernels and 
%integration scemes.
Kernels;

%Parameter for the Number of approximation steps
N = 45;

%Generation of different Kernels
BrownianMotion = Kernels.KMat(1,N);
BrownianBridge = Kernels.KMat(2,N);
ExponentialKer = Kernels.KMat(3,N);

%BROWNIAN MOTION
%Solve Fredhol integralequality with different Scemes
[lambda1,Phi1] = Kernels.uniform_Sceme(BrownianMotion);
[lambda2,Phi2] = Kernels.trapez_Sceme(BrownianMotion);
[lambda3,Phi3] = Kernels.simpson_Sceme(BrownianMotion);
%Compute analytic solutions for first Eigenvalues
lambda = [lambda1(1) lambda2(1) lambda3(1)];
Phi = [Phi1(:,1) Phi2(:,1) Phi3(:,1)];
la = (2/pi)^2;
ph = zeros(N+2,1);
for i=1:N+2
   ph(i) = sqrt(2)* sin(0.5*pi*((i-1)/(N+2))); 
end
%plot(linspace(0,1,N+2),ph,linspace(0,1,N+2),Phi(:,1))
%Compute the error terms
absolute_error_lambda = abs(la(1)-lambda)
relative_error_lambda = abs(la(1)-lambda)/la(1)*100
absolute_error_phi = zeros(1,3);
relative_error_phi = zeros(1,3);
for i=1:3
     absolute_error_phi(i) = max(abs(ph-Phi(:,i)));
     relative_error_phi(i) = absolute_error_phi(i)/max(abs(Phi(:,i)));
end
absolute_error_phi
relative_error_phi*100

%BROWNIAN BRIDGE
%Solve Fredhol integralequality with different Scemes
[lambda1,Phi1] = Kernels.uniform_Sceme(BrownianBridge);
[lambda2,Phi2] = Kernels.trapez_Sceme(BrownianBridge);
[lambda3,Phi3] = Kernels.simpson_Sceme(BrownianBridge);
%Compute analytic solutions for first Eigenvalues
lambda = [lambda1(1) lambda2(1) lambda3(1)];
Phi = [Phi1(:,1) Phi2(:,1) Phi3(:,1)];
la = (1/pi)^2;
ph = zeros(N+2,1);
for i=1:N+2
   ph(i) = sqrt(2)* sin(pi*((i-1)/(N+2))); 
end
plot(linspace(0,1,N+2),ph,linspace(0,1,N+2),Phi(:,1))
%Compute the error terms
absolute_error_lambda = abs(la(1)-lambda)
relative_error_lambda = abs(la(1)-lambda)/la(1)*100
absolute_error_phi = zeros(1,3);
relative_error_phi = zeros(1,3);
for i=1:3
     absolute_error_phi(i) = max(abs(ph-Phi(:,i)));
     relative_error_phi(i) = absolute_error_phi(i)/max(abs(Phi(:,i)));
end
absolute_error_phi
relative_error_phi*100



