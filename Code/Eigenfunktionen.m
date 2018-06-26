%%This Programm generates a Plot eigenvalues and eigenfunctions of 
%%the induced integral operators of different kernels
%%written by Tim Jaschek as a part of his bachelor thesis%%

%%This Programm is used to generate FIGURE 2 in the thesis%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the class Kernels
Kernels;

%Parameter for accuracy
N=50;

figure
for i=1:3;
    Mat = Kernels.KMat(i,N);
    [lambda,Phi] = Kernels.trapez_Sceme(Mat);
    subplot(3,3,[1+3*(i-1) 2+3*(i-1)]);
    for j=1:6
        hold on;
        plot(linspace(0,1,N+2),Phi(:,j));
    end
    if i ==1
        title('First 6 Eigenfunctions');
        ylabel('K(s,t)=min(s,t)');
    elseif i == 2
        ylabel('K(s,t)=min(s,t) - st');
    else
        ylabel('K(s,t)=exp(-|s-t|)');
    end
    hold off;
    subplot(3,3,3*i);
    plot(linspace(1,10,10),lambda(1:10),'o','color','red');
    if i ==1
        title('First 10 Eigenvalues');
    end
end



