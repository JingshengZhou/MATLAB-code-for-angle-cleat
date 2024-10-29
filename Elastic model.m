clc
clear all
%geometry
gb=70; %beam bolt gauge
gc=50; %column bolt gauge
t=9.;  %leg thickness
r=9.84; %root radius
b=160; %width of the angle
%material
E=205000; %elastic stiffness
fy=338;   %yield stress
fff=100;  %friction force
my=b*fy*t^2/6; %cross-sectional yield moment

I=b*t^3/12; %moment of inertia of a rectangle

g1=gb-t-r*0.8-24/4; %gb,m
g2=gc-t-r*0.8-24/4; %gc,m
z=t/2+r*0.8;        %z
stif=zeros(1,10000);
dert=zeros(1,10000);
L=zeros(1,10000);

a=1;

for thetar=0.001:0.001:0.006 %give a rotation of root 
    for dertapcatry=0:0.00001:2 %select a slipage of bolt at point A
        dertac=z*(1-cos(thetar)+sin(thetar))-dertapcatry; %horizontal displacement of point C
        mc=2*E*I/g1*(2*thetar+3*dertac/g1); %moment at point C
        md=2*E*I/g1*(thetar+3*dertac/g1);   %momnet at point D
        R=(mc+md)/g1;                       %force R
        if dertapcatry<0.000000000001       %distinguish the correctness of the slipage of bolt at point A
            if R<fff
                dertapca=0;
            end
        else
            if dertapcatry<2
                abs(R-fff)/R;
                if abs(R-fff)/fff<0.01
                    dertapca=dertapcatry;
                end
            else
                if R>fff
                    dertapca=2;
                end
            end
        end
    end
    dertac=z*(1-cos(thetar)+sin(thetar))-dertapca;
    mc=2*E*I/g1*(2*thetar+3*dertac/g1);
    md=2*E*I/g1*(thetar+3*dertac/g1);
    R=(mc+md)/g1;
    for dertabtry=0:0.00000001:2                %select a vertical displacement of point B
        ma=2*E*I/g2*(-thetar+3*dertabtry/g2);   %moment at point A
        mb=2*E*I/g2*(-2*thetar+3*dertabtry/g2); %moment at point V
        f1=(ma+mb)/g2;                          %total force
        f2=(mc-mb+R*z*(cos(thetar)+sin(thetar)))/z/(cos(thetar)-sin(thetar)); %total force
        if (f1-f2)/f1<0.01  %distinguish the correctness of the selected vertical displacement of point B
            f=f1;
            dertab=dertabtry;
        end
    end
    ma=2*E*I/g2*(-thetar+3*dertab/g2);
    mb=2*E*I/g2*(-2*thetar+3*dertab/g2);
   
    derta=dertab+z*(cos(thetar)+sin(thetar)-1);%total deformation of the angle
    kangle=2*f/derta/1000; %stiffness of the angle
    f/1000;
    %distinguish the formation of the plastic hinge:
    if ma>my  
        f*2/1000
        f=0;
    end    
    if mb>my
        f*2/1000
        f=0;
    end
    if mc>my
        f*2/1000
        f=0;
    end
    if md>my
        f*2/1000
        f=0;
    end
    stif(a)=k;
    L(a)=f/1000;

    dert(a)=derta;
    a=a+1;
end                    
   
plot(dert,L);

