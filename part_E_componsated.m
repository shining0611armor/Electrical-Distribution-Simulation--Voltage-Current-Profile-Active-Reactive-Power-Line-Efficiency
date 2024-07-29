
%-------------------part E
%-------------------- ********compensation in series and parallel form

z=0.0165+0.3306*1i; %z in kilometer
y=(4.674*(10^-6))*1i; %y in kilometer
%----componsation procedure
y_prime=y*(1-0.8) ;%   80% compensated in parallel
Z_prime_iamg=imag(z)*(1-0.2) ;%   20% compensated in series
Z_prime_real=real(z) ;
Z_prime=Z_prime_real+Z_prime_iamg*1i;
Zc=sqrt(Z_prime/y_prime);%  in ohm
gama=sqrt(Z_prime*y_prime);% in  one over kilometer unit
%-------------------------
VLL_rated=400*10^3;
VLoad_ln=VLL_rated/sqrt(3); % load volatage in line to nutral
L=400;% inkilometer
PL=600*10^6;%active power
cosphi=0.99;
angle_of_cosphi=-acosd(cosphi); %power factor lag (determined in it's sign)
gamma_L=gama*L;% gamma & length production in pu
vr=VLoad_ln;% load voltage line to nutral
ir_mag=(PL/(abs(cosphi)*VLoad_ln*3));% load current
ir=(ir_mag * cos(angle_of_cosphi*(pi/180)) + ir_mag * sin(angle_of_cosphi*(pi/180))*1i); % load current in Cartesian



X=0:0.01:L;% variable lengh in the lenght is zero in vicinity of  the load
gama_X =gama*X;

%---------------------------- calculation of ABCD parameter

A=cosh(gama_X);
B=Zc*sinh(gama_X);
C=(1/Zc)*sinh(gama_X);
D=A;

%---------------------------- calculation of current and voltage over entire line
%Vs=voltage over entire line
%Is=current over entire line

vs=A*vr + B*ir;% sending line_to_nutral voltage at sending end
is=C*vr + D*ir;

vs_mag=abs(vs);
is_mag=abs(is);


is_angle=zeros(1,length(X));
vs_angle=zeros(1,length(X));
for i=1:length(X)          % assiging in the loop, for the sake of inner matrix dimensions must agree.
is_angle(i)=atand(imag(is(i))/real(is(i)));% construction of is angle matrix
vs_angle(i)=atand(imag(vs(i))/real(vs(i)));% construction of Vs angle matrix
end

%----------------------------- calculation of active and reactive power over the line

angle_VS_IS =vs_angle-is_angle;
cosine_Vs_Is_angle=cosd(angle_VS_IS);
sine_Vs_Is_angle=sind(angle_VS_IS);


P=zeros(1,length(X));% construction of P matrix
Q=zeros(1,length(X));% construction of Q matrix
for i=1:length(X)           % assiging in the loop, for the sake of inner matrix dimensions must agree.
P(i)=3*vs_mag(i)*is_mag(i)*cosine_Vs_Is_angle(i);
Q(i)=3*vs_mag(i)*is_mag(i)*sine_Vs_Is_angle(i);
end


%------------------ bazdeh khat
PS=P(length(X)); %sening end power
randeman=(P/PS)*100; % each power over power in teh seding end (it plased on lenghth(X) argument)

figure('Name','partE - computer assignment - componsated','NumberTitle','off')

subplot(3,2,1);plot(X,vs_mag);set(gca,'XDir','rev');% just reverse x axis in standard form    
title('voltage profile')

subplot(3,2,2);plot(X,is_mag);set(gca,'XDir','rev')% just reverse x axis in standard form 
title('current profile')

subplot(3,2,3);plot(X,P);set(gca,'XDir','rev');% just reverse x axis in standard form    
title('active power over line')

subplot(3,2,4);plot(X,Q);set(gca,'XDir','rev')% just reverse x axis in standard form 
title('reactive power over line')


subplot(3,2,5);plot(X,randeman);set(gca,'XDir','rev')% just reverse x axis in standard form 
title(' bazdeh khat dar noghate mokhtalef khat')

