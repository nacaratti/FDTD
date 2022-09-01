%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       FDTD - 1D - TNM-5845                            %
% Exercício: Escrever um código para a propagação de um pulso Gaussiano %
% Davi Pontes Nacaratti                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

ep0 = 1;                              % Permissividade
mu0 = 1;                              % Permeabilidade magnetica
%c = 1/sqrt(ep0*mu0);                 % Velocidade da luz no vacuo
c=1;                                  
lambda0 = 1.064e-6;                   % Comprimento de onda laser Yb
n = 1.45;                             % Indice de refracao            
L = 30e-6;                            % Comprimento em z
N = 100;                              % Numero de series em z
ep = ones(N,1);                       % Vetor de Permissividade         
mu = ones(N,1);                       % Vetor de Permeabilidade magnetica      
dz = 0.1;                             % Incremento em z
z = (0:(N-1))*dz;                     % Vetor z
dt = 0.01;                            % Incremento t              
Z = sqrt(mu0/ep0);                    % Impedancia E'=E/Z
Hy = zeros(N,1);                      % Zerando o vetor Hy
Ex = zeros(N,1);                      % Zerando o vetor Ex

% Condição para o incremento deltaz
if (lambda0/dz) >= (10*n)
    fprintf('Condição N=(lambda/h)>=10n satisfeita!!');
else
    fprintf('Condição N=(lambda/h)>=10n NAO satisfeita!!');
end

% Editando o perfil de permissividade
% for i = 1:20
%     ep(i,1) = 1.0; 
% end  
%   for i = 20+1:50
%       ep(i,1) = 2.5; 
%   end  
% figure(1);
% plot(ep);
% title('Perfil permissividade');

% Criando as matrizes
Df = spdiags([-mu mu],0:1,N,N)/dz;
Df(N,1) = Df(1,2);
Db = spdiags([-ep ep],-1:0,N,N)/dz; 
Db(1,N) = Db(2,1);
figure(2);
subplot(1,2,1);
imagesc(Df);
subplot(1,2,2);
imagesc(Db);

% Pulso Gaussiano (hard source)
Ex = exp(-((z'-5)).^2);

Nt = 1000;

% Main Looping FDTD using diag matrix
  for i=1:Nt
      Hy = Hy - Df*dt*Ex;
      Ex = Ex - Db*dt*Hy;
      figure(3);
      subplot(2,1,1);
      plot(z,Ex,'b');           
      %axis([0 N -1 1]);
      grid on
      xlabel('z');
      ylabel('Ex');
      subplot(2,1,2);
      plot(z,Hy,'r');      
      grid on
      %axis([0 N -1 1]);
      xlabel('z');
      ylabel('Hy');
      drawnow;
  end
 
 % Fim do programa
