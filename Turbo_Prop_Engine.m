clear all; clc;
h = 0;
T04 = 1389;
Q=43960000;
gm=1.4;
m=13.61;
M = 0.7;
a=6.5;
Tssl=288.15;
R=287;
Pa = 101325;
Cwp=1.0079;
Ta =Tssl - a*h;
SS = sqrt(gm*R*Ta);
Toub=zeros(1,24);
TSFC = zeros(1,24);
Tta = Ta*(1+(((gm-1)/2)*(M^2))); 
Pta = Pa*(1+(((gm-1)/2)*(M^2)))^(gm/(gm-1)); 
Cp=gm*R/(gm-1);
Ua = M*SS;
W_prop = Cwp*m*Cp*Ta;

Pi_c = 1:1:24;
T05 = zeros(1,24);
tou5 = zeros(1,24);
Cwe = zeros(1,24);
B = zeros(1,24);
SFC = zeros(1,24);
Th_sf = zeros(1,24);
Th = zeros(1,24);
Cwjet = zeros(1,24);
T03 = zeros(1,24);
mf = zeros(1,24);
P03 = zeros(1,24);
P04 = zeros(1,24);
P05 = zeros(1,24);
Me = zeros(1,24);
Te = zeros(1,24);
Ue = zeros(1,24); 
P_nz = zeros(1,24);
n_prop = zeros(1,24);
n_th = zeros(1,24);
n_tol = zeros(1,24);
P = zeros(1,24);

for i= 1:1:24
 T05(i) = T04 - ((W_prop + (m*Cp*Tta*(Pi_c(i)^(2/7)-1)))/(m*Cp)); 
 tou5(i) = (Tta/Ta)*(Pi_c(i)^(2/7))*(T05(i)/T04);
 B(i) = ((T04/Ta) * (Ta/Tta) * ((tou5(i)-1) / (Pi_c(i)^(2/7)))) / ((Tta/Ta)-1);
 Cwjet(i) = (gm-1)*(M^2)*(sqrt(B(i))-1);
 Cwe(i) = Cwp + real(Cwjet(i));
 Th(i) = Cwe(i)*m*Cp*Ta/Ua;
 P(i) = Cwe(i)*m*Cp*Ta;
 Th_sf(i) = Th(i)/m;
 T03(i) = Tta*Pi_c(i)^(2/7);
 mf(i) = m*Cp*(T04-T03(i))/Q;
 TSFC(i) = mf(i)/Th(i);
 
 %propulsive effeciency
 %Compressor
 P03(i) = Pta*Pi_c(i);
 T03(i) = Tta*(Pi_c(i))^((gm-1)/gm);
 
 %Burner
 P04(i) = P03(i);
 
 %Turbine
 P05(i) = P04(i) * (T05(i)/(T04))^(gm/(gm-1));
 
 %nozzel
 P_nz(i) = ((P05(i)/Pa)^(0.2857) - 1);
 Me(i) = sqrt( (2/(gm-1)) * P_nz(i) );
 Te(i) = T05(i) / (1+ ((gm-1)/2)*Me(i)^2);
 Ue(i) = Me(i)*sqrt(gm*R*Te(i));
 
 %propulsive effeciency
 %n_prop(i) = (P(i))/((m+mf(i))*Ue(i)^2 - (m)*Ua^2); 
 n_prop(i) = 2*Ua/(Ua + Ue(i));
 
 %Thermal effeciency
 n_th(i)= P(i)/(Q*mf(i));
 
 %Thermal effeciency
 n_tol(i) = n_prop(i)*n_th(i);
 
end

plot(Pi_c(1,2:end),TSFC(1,2:end),'-*')
xlabel('Compression Pressure Ratio') 
ylabel('TSFC') 
figure;
plot(Pi_c(1,2:end),Th_sf(1,2:end),'-o')
xlabel('Compression Pressure Ratio') 
ylabel('Specific Thrust(N/kg)') 

figure;
plot(Pi_c(1,3:end),n_th(1,3:end),'-s')
xlabel('Compression Pressure Ratio') 
hold on
plot(Pi_c(1,3:end),n_prop(1,3:end),'+-')
plot(Pi_c(1,3:end),n_tol(1,3:end),'diamond-')
hold off
legend({'Thermal effeciency','propulsive effeciency','Total effeciency'},'Location','best')

