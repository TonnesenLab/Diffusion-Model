clc
clear all
load('Tensormap_TM05.mat')

n=4;

theta=THETA{n};
angles=theta.*180/pi;
angles(angles>180)=angles(angles>180)-180;

puntas=PUNTAS{n};
mag=sqrt(puntas(:,1).^2 + puntas(:,2).^2);
mag=mag.*0.06;

angulo_medio=angulo_medio_axial(angles);