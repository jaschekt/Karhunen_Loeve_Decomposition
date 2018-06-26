%%This class generates a plot of a Brownian sheet%%
%%written by Tim Jaschek as a part of his bachelor thesis%%

%%Used to generate FIGURE 7 %%
%%...to generate it, type the following in your MATLAB command:
%%Sheet;
%%Sheet.plotit();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Sheet
   properties (Constant)
       N = 100;
       n = 20;
   end
   methods (Static)
       function lam = sqlambda(i,j)
           lam = 4 / ((2*j-1)*(2*i-1) * pi^2);
       end
       function ph = phi(i,j)
           ph = zeros(Sheet.N,Sheet.N);
           for k=1:Sheet.N
               for l=1:Sheet.N
                    ph(k,l) = 2*sin((i-0.5)*pi*k/Sheet.N)*sin((j-0.5)*pi*l/Sheet.N);
               end
           end
       end
       function plotit()
           BS = zeros(Sheet.N,Sheet.N);
           BS2 = zeros(Sheet.N,Sheet.N);
           xi = randn(1,Sheet.n^2);
           figure
           for i=1:5
               for j=1:5
                   lam = Sheet.sqlambda(i,j);
                   phi = Sheet.phi(i,j);
                   BS2 = BS2 + lam*phi*xi(Sheet.n*(i-1)+j);
               end
               i
           end
           subplot(3,1,1);
           surf(linspace(1/Sheet.N,1,Sheet.N),linspace(1/Sheet.N,1,Sheet.N),BS2,'edgealpha','0');
           colormap jet
           tic;
           for i=1:Sheet.n
               for j=1:Sheet.n
                   lam = Sheet.sqlambda(i,j);
                   phi = Sheet.phi(i,j);
                   BS = BS + lam*phi*xi(Sheet.n*(i-1)+j);
               end
               i
           end
           toc
           subplot(3,1,[2,3])
           surf(linspace(1/Sheet.N,1,Sheet.N),linspace(1/Sheet.N,1,Sheet.N),BS,'edgealpha','0');
           colormap jet
       end
   end
end
