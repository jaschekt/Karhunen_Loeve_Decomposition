%%This Programm generates a Plot for different Kernels and partial%%
%%sums of Mercer's series for the covariance function%%
%%written by Tim Jaschek as a part of his bachelor thesis%%

%%This Programm is used to generate FIGURE 1 in the thesis%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the class Kernels
Kernels;

%Parameter for accuracy
N=16;

figure
for i=1:3;
    Mat = Kernels.KMat(i,N);
    [lambda,Phi] = Kernels.trapez_Sceme(Mat);
    for j=1:3;
        K=Kernels.MercerApprox(lambda,Phi,j+(j-1)^2);
        subplot(4,3,i+3*(j-1));
        surfc(linspace(0,1,N+2),linspace(0,1,N+2),K,'edgealpha','1');
        if j == 1
            if i ==1
                title('K(s,t)=min(s,t)');
                zlabel('n=1');
            elseif i ==2
                title('K(s,t)=min(s,t) - st');
            else
                title('K(s,t)=exp(-|s-t|)');
            end
        elseif j== 2
            if i == 1
                zlabel('n=3');
            end
        else
            if i == 1
                zlabel('n=7');
            end
        end
    end
    subplot(4,3,i+9);
    surfc(linspace(0,1,N+2),linspace(0,1,N+2),Mat,'edgealpha','1');
    if i == 1
        zlabel('analytic');
    end
end



