clc
clear all
%geometry
gb=50;
gc=50;
t=9.;
r=9.84;
b=160;
%material
E=205000;
Eh=1703;
Eu=537;
fy=337.89;
fu=600.11;
ff=799.8;
xy=fy/E;
xh=0.015;
xm=0.1689;
xu=0.54;
I=b*t^3/12;
fff=1;
elastic=0.2;
%calculation

cy=2*xy/t;
ch=2*xh/t;
cm=2*xm/t;
cu=2*xu/t;

g1=gb-t-r*0.8-6;
g2=gc-t-r*0.8-6;
z=t/2+r*0.8;

% moment
my=b*fy*t^2/6;

mh=0.5*my*(3-(xy/xh)^2);
mm=0.5*my*(3-(xy/xm)^2)+my/2*Eh/E*(xm-xh)/cy*(1-xh/xm)*(2+xh/xm);
mu=0.5*my*(3-(xy/xu)^2)+my/2*Eh/E*(xu-xh)/cy*(1-xh/xu)*(2+xh/xu)-my/2*(Eh-Eu)/E*(xu-xm)/xy*(1-xm/xu)*(2+xm/xu);

%stiffness of spring
kh=(mh-my)/2/(xh-xy);
km=(mm-mh)/2/(xm-xh);
ku=(mu-mm)/2/(xu-xm);

%settings
a=1;
L=zeros(1,10000);
dert=zeros(1,10000);
cd1=1;
cd2=1;
cd3=1;

for thetar=0.2:0.1:0.2
    for rfatry=1:0.01:1
        for rfbtry=1:0.01:1
            for rfctry=1:0.01:1
                for rfdtry=0.01:0.01:1
                    for thetactry=thetar:0.001:(xu-xh)*2
                        thetad=thetactry-thetar;
                        mc1=calmoment(thetactry,rfctry,my,xm,xh,Eh,E,xy,Eu);
                        md=calmoment(thetad,rfdtry,my,xm,xh,Eh,E,xy,Eu);
                        dertac=g1*sin(thetad);
                        for dertapcatry=elastic:0.0001:2-elastic
                            thetaa=acos(1-(dertac+dertapcatry-z*(1-cos(thetar)+sin(thetar)))/g2);
                            dertab=g2*sin(thetaa);
                            thetab=thetaa-thetar;
                            ma=calmoment(thetaa,rfatry,my,xm,xh,Eh,E,xy,Eu);
                            mb=calmoment(thetab,rfbtry,my,xm,xh,Eh,E,xy,Eu);
                            f=((mc1+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                            n1=((mc1+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);

                            if dertapcatry<0.000000001
                                if n1<fff
                                    dertapca=0;
                                end
                            else
                                if dertapcatry<2-elastic
                                    abs(n1-fff)/fff;
                                    if abs(n1-fff)/fff<0.01
                                        dertapca=dertapcatry;
                                    end
                                else
                                    if n1>fff
                                        dertapca=2-elastic;
                                    end
                                end
                            end
                        end
                        thetaa=acos(1-(dertac+dertapca-z*(1-cos(thetar)+sin(thetar)))/g2);
                        thetab=thetaa-thetar;
                        ma=calmoment(thetaa,rfatry,my,xm,xh,Eh,E,xy,Eu);
                        mb=calmoment(thetab,rfbtry,my,xm,xh,Eh,E,xy,Eu);
                        dertab=g2*sin(thetaa);
                        f=((mc1+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                        n=((mc1+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                        mc2=mb+f*z*(cos(thetar)-sin(thetar))-n*z*(cos(thetar)+sin(thetar));
                        abs(mc2-mc1)/mc1;
                        if abs(mc2-mc1)/mc2<0.01
                            if abs(mc2-mc1)/mc1<0.01
                                thetactry;
                                if cd1<1.1
                                    thetactry
                                    thetac=thetactry;
                                    cd1=2;
                                end
                            end
                        end
                    end       
                    thetad=thetac-thetar;
                    mc=calmoment(thetac,rfctry,my,xm,xh,Eh,E,xy,Eu);
                    md=calmoment(thetad,rfdtry,my,xm,xh,Eh,E,xy,Eu);
                    dertac=g1*sin(thetad);
                    for dertapcatry=2-elastic:0.0001:2-elastic
                        thetaa=acos(1-(dertac+dertapcatry-z*(1-cos(thetar)+sin(thetar)))/g2);
                        dertab=g2*sin(thetaa);
                        thetab=thetaa-thetar;
                        ma=calmoment(thetaa,rfatry,my,xm,xh,Eh,E,xy,Eu);
                        mb=calmoment(thetab,rfbtry,my,xm,xh,Eh,E,xy,Eu);
                        f=((mc+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                        n1=((mc+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                        if dertapcatry<0.000000001
                            if n1<fff
                                dertapca=0;
                            end
                        else
                            if dertapcatry<2-elastic
                                abs(n1-fff)/fff;
                                if abs(n1-fff)/fff<0.01
                                    dertapca=dertapcatry;
                                end
                            else
                                if n1>fff
                                    dertapca=2-elastic;
                                end
                            end
                        end
                    end
                    thetaa=acos(1-(dertac+dertapca-z*(1-cos(thetar)+sin(thetar)))/g2);
                    thetab=thetaa-thetar;
                    ma=calmoment(thetaa,rfatry,my,xm,xh,Eh,E,xy,Eu);
                    mb=calmoment(thetab,rfbtry,my,xm,xh,Eh,E,xy,Eu);
                    dertab=g2*sin(thetaa);
                    f=((mc+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                    n=((mc+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                    stressdb=calstress(thetad,xm,xh,Eh,Eu,fy,fu);                    
                    stressdt=f/t/b/rfdtry;
                    stressdp=stressdb+stressdt;
                    stressdn=stressdb-stressdt;
                    straindp=calstrain(stressdp,xm,xh,Eh,Eu,fy,fu);
                    straindn=calstrain(stressdn,xm,xh,Eh,Eu,fy,fu);
                    %rfd1=(1-(straindp-straindn)/2);
                    rfd1=1/(straindp+straindn)*(1/exp(straindp)-1/exp(straindp*2+straindn));
                    if abs(rfd1-rfdtry)/rfdtry<0.01
                        rfd=rfdtry;
                    end
                end
                for thetactry=thetar:0.001:(xu-xh)*2
                    thetad=thetactry-thetar;
                    mc1=calmoment(thetactry,rfctry,my,xm,xh,Eh,E,xy,Eu);
                    md=calmoment(thetad,rfd,my,xm,xh,Eh,E,xy,Eu);
                    dertac=g1*sin(thetad);
                    for dertapcatry=2-elastic:0.0001:2-elastic
                        thetaa=acos(1-(dertac+dertapcatry-z*(1-cos(thetar)+sin(thetar)))/g2);
                        dertab=g2*sin(thetaa);
                        thetab=thetaa-thetar;
                        ma=calmoment(thetaa,rfatry,my,xm,xh,Eh,E,xy,Eu);
                        mb=calmoment(thetab,rfbtry,my,xm,xh,Eh,E,xy,Eu);
                        f=((mc1+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                        n1=((mc1+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);

                        if dertapcatry<0.000000001
                            if n1<fff
                                dertapca=0;
                            end
                        else
                            if dertapcatry<2-elastic
                                abs(n1-fff)/fff;
                                if abs(n1-fff)/fff<0.01
                                    dertapca=dertapcatry;
                                end
                            else
                                if n1>fff
                                    dertapca=2-elastic;
                                end
                            end
                        end
                    end
                    thetaa=acos(1-(dertac+dertapca-z*(1-cos(thetar)+sin(thetar)))/g2);
                    thetab=thetaa-thetar;
                    ma=calmoment(thetaa,rfatry,my,xm,xh,Eh,E,xy,Eu);
                    mb=calmoment(thetab,rfbtry,my,xm,xh,Eh,E,xy,Eu);
                    dertab=g2*sin(thetaa);
                    f=((mc1+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                    n=((mc1+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                    mc2=mb+f*z*(cos(thetar)-sin(thetar))-n*z*(cos(thetar)+sin(thetar));
                    abs(mc2-mc1)/mc1;
                    if abs(mc2-mc1)/mc2<0.01
                        if abs(mc2-mc1)/mc1<0.01
                            thetactry;
                            if cd2<1.1
                                thetactry
                                thetac=thetactry;
                                cd2=2;
                            end
                        end
                    end
                end       
                thetad=thetac-thetar;
                mc=calmoment(thetac,rfctry,my,xm,xh,Eh,E,xy,Eu);
                md=calmoment(thetad,rfd,my,xm,xh,Eh,E,xy,Eu);
                dertac=g1*sin(thetad);
                for dertapcatry=2-elastic:0.0001:2-elastic
                    thetaa=acos(1-(dertac+dertapcatry-z*(1-cos(thetar)+sin(thetar)))/g2);
                    dertab=g2*sin(thetaa);
                    thetab=thetaa-thetar;
                    ma=calmoment(thetaa,rfatry,my,xm,xh,Eh,E,xy,Eu);
                    mb=calmoment(thetab,rfbtry,my,xm,xh,Eh,E,xy,Eu);
                    f=((mc+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                    n1=((mc+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                    if dertapcatry<0.000000001
                        if n1<fff
                            dertapca=0;
                        end
                    else
                        if dertapcatry<2-elastic
                            abs(n1-fff)/fff;
                            if abs(n1-fff)/fff<0.01
                                dertapca=dertapcatry;
                            end
                        else
                            if n1>fff
                                dertapca=2-elastic;
                            end
                        end
                    end
                end
                thetaa=acos(1-(dertac+dertapca-z*(1-cos(thetar)+sin(thetar)))/g2);
                thetab=thetaa-thetar;
                ma=calmoment(thetaa,rfatry,my,xm,xh,Eh,E,xy,Eu);
                mb=calmoment(thetab,rfbtry,my,xm,xh,Eh,E,xy,Eu);
                dertab=g2*sin(thetaa);
                f=((mc+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                n=((mc+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);

                stresscb=calstress(thetac,xm,xh,Eh,Eu,fy,fu);
                stressct=f/t/b/rfctry;
                stresscp=stresscb+stressct;
                stresscn=stresscb-stressct;
                straincp=calstrain(stresscp,xm,xh,Eh,Eu,fy,fu);
                straincn=calstrain(stresscn,xm,xh,Eh,Eu,fy,fu);
                rfc1=(1-(straincp-straincn)/2);
                if abs(rfc1-rfctry)/rfctry<0.01
                    rfc=rfctry;
                end
            end
            for thetactry=thetar:0.001:(xu-xh)*2
                thetad=thetactry-thetar;
                mc1=calmoment(thetactry,rfc,my,xm,xh,Eh,E,xy,Eu);
                md=calmoment(thetad,rfd,my,xm,xh,Eh,E,xy,Eu);
                dertac=g1*sin(thetad);
                for dertapcatry=2-elastic:0.0001:2-elastic
                    thetaa=acos(1-(dertac+dertapcatry-z*(1-cos(thetar)+sin(thetar)))/g2);
                    dertab=g2*sin(thetaa);
                    thetab=thetaa-thetar;
                    ma=calmoment(thetaa,rfatry,my,xm,xh,Eh,E,xy,Eu);
                    mb=calmoment(thetab,rfbtry,my,xm,xh,Eh,E,xy,Eu);
                    f=((mc1+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                    n1=((mc1+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);

                    if dertapcatry<0.000000001
                        if n1<fff
                            dertapca=0;
                        end
                    else
                        if dertapcatry<2-elastic
                            abs(n1-fff)/fff;
                            if abs(n1-fff)/fff<0.01
                                dertapca=dertapcatry;
                            end
                        else
                            if n1>fff
                                dertapca=2-elastic;
                            end
                        end
                    end
                end
                thetaa=acos(1-(dertac+dertapca-z*(1-cos(thetar)+sin(thetar)))/g2);
                thetab=thetaa-thetar;
                ma=calmoment(thetaa,rfatry,my,xm,xh,Eh,E,xy,Eu);
                mb=calmoment(thetab,rfbtry,my,xm,xh,Eh,E,xy,Eu);
                dertab=g2*sin(thetaa);
                f=((mc1+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                n=((mc1+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                mc2=mb+f*z*(cos(thetar)-sin(thetar))-n*z*(cos(thetar)+sin(thetar));
                abs(mc2-mc1)/mc1;
                if abs(mc2-mc1)/mc2<0.01
                    if abs(mc2-mc1)/mc1<0.01
                        thetactry;
                        if cd3<1.1
                            thetactry
                            thetac=thetactry;
                            cd3=2;
                        end
                    end
                end
            end       
            thetad=thetac-thetar;
            mc=calmoment(thetac,rfc,my,xm,xh,Eh,E,xy,Eu);
            md=calmoment(thetad,rfd,my,xm,xh,Eh,E,xy,Eu);
            dertac=g1*sin(thetad);
            for dertapcatry=2-elastic:0.0001:2-elastic
                thetaa=acos(1-(dertac+dertapcatry-z*(1-cos(thetar)+sin(thetar)))/g2);
                dertab=g2*sin(thetaa);
                thetab=thetaa-thetar;
                ma=calmoment(thetaa,rfatry,my,xm,xh,Eh,E,xy,Eu);
                mb=calmoment(thetab,rfbtry,my,xm,xh,Eh,E,xy,Eu);
                f=((mc+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                n1=((mc+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                if dertapcatry<0.000000001
                    if n1<fff
                        dertapca=0;
                    end
                else
                    if dertapcatry<2-elastic
                        abs(n1-fff)/fff;
                        if abs(n1-fff)/fff<0.01
                            dertapca=dertapcatry;
                        end
                    else
                        if n1>fff
                            dertapca=2-elastic;
                        end
                    end
                end
            end
            thetaa=acos(1-(dertac+dertapca-z*(1-cos(thetar)+sin(thetar)))/g2);
            thetab=thetaa-thetar;
            ma=calmoment(thetaa,rfatry,my,xm,xh,Eh,E,xy,Eu);
            mb=calmoment(thetab,rfbtry,my,xm,xh,Eh,E,xy,Eu);
            dertab=g2*sin(thetaa);
            f=((mc+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
            n=((mc+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);

            stressbb=calstress(thetab,xm,xh,Eh,Eu,fy,fu);                
            stressbt=n/t/b/rfbtry;
            stressbp=stressbb+stressbt;
            stressbn=stressbb-stressbt;
            strainbp=calstrain(stressbp,xm,xh,Eh,Eu,fy,fu);
            strainbn=calstrain(stressbn,xm,xh,Eh,Eu,fy,fu);
            rfb1=(1-(strainbp-strainbn)/2);
            if abs(rfb1-rfbtry)/rfbtry<0.01
                rfb=rfbtry;
            end
        end
        for thetactry=thetar:0.001:(xu-xh)*2
            thetad=thetactry-thetar;
            mc1=calmoment(thetactry,rfc,my,xm,xh,Eh,E,xy,Eu);
            md=calmoment(thetad,rfd,my,xm,xh,Eh,E,xy,Eu);
            dertac=g1*sin(thetad);
            for dertapcatry=2-elastic:0.0001:2-elastic
                thetaa=acos(1-(dertac+dertapcatry-z*(1-cos(thetar)+sin(thetar)))/g2);
                dertab=g2*sin(thetaa);
                thetab=thetaa-thetar;
                ma=calmoment(thetaa,rfatry,my,xm,xh,Eh,E,xy,Eu);
                mb=calmoment(thetab,rfb,my,xm,xh,Eh,E,xy,Eu);
                f=((mc1+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
                n1=((mc1+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);

                if dertapcatry<0.000000001
                    if n1<fff
                        dertapca=0;
                    end
                else
                    if dertapcatry<2-elastic
                        abs(n1-fff)/fff;
                        if abs(n1-fff)/fff<0.01
                            dertapca=dertapcatry;
                        end
                    else
                        if n1>fff
                            dertapca=2-elastic;
                        end
                    end
                end
            end
            thetaa=acos(1-(dertac+dertapca-z*(1-cos(thetar)+sin(thetar)))/g2);
            thetab=thetaa-thetar;
            ma=calmoment(thetaa,rfatry,my,xm,xh,Eh,E,xy,Eu);
            mb=calmoment(thetab,rfb,my,xm,xh,Eh,E,xy,Eu);
            dertab=g2*sin(thetaa);
            f=((mc1+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
            n=((mc1+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
            mc2=mb+f*z*(cos(thetar)-sin(thetar))-n*z*(cos(thetar)+sin(thetar));
            abs(mc2-mc1)/mc1;
            if abs(mc2-mc1)/mc2<0.01
                if abs(mc2-mc1)/mc1<0.01
                    thetactry;
                    if cd4<1.1
                        thetactry
                        thetac=thetactry;
                        cd4=2;
                    end
                end
            end
        end       
        thetad=thetac-thetar;
        mc=calmoment(thetac,rfc,my,xm,xh,Eh,E,xy,Eu);
        md=calmoment(thetad,rfd,my,xm,xh,Eh,E,xy,Eu);
        dertac=g1*sin(thetad);
        for dertapcatry=2-elastic:0.0001:2-elastic
            thetaa=acos(1-(dertac+dertapcatry-z*(1-cos(thetar)+sin(thetar)))/g2);
            dertab=g2*sin(thetaa);
            thetab=thetaa-thetar;
            ma=calmoment(thetaa,rfatry,my,xm,xh,Eh,E,xy,Eu);
            mb=calmoment(thetab,rfb,my,xm,xh,Eh,E,xy,Eu);
            f=((mc+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
            n1=((mc+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
            if dertapcatry<0.000000001
                if n1<fff
                    dertapca=0;
                end
            else
                if dertapcatry<2-elastic
                    abs(n1-fff)/fff;
                    if abs(n1-fff)/fff<0.01
                        dertapca=dertapcatry;
                    end
                else
                    if n1>fff
                        dertapca=2-elastic;
                    end
                end
            end
        end
        thetaa=acos(1-(dertac+dertapca-z*(1-cos(thetar)+sin(thetar)))/g2);
        thetab=thetaa-thetar;
        ma=calmoment(thetaa,rfatry,my,xm,xh,Eh,E,xy,Eu);
        mb=calmoment(thetab,rfb,my,xm,xh,Eh,E,xy,Eu);
        dertab=g2*sin(thetaa);
        f=((mc+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
        n=((mc+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);

        stressab=calstress(thetaa,xm,xh,Eh,Eu,fy,fu);
        stressat=n/t/b/rfatry;
        stressap=stressab+stressat;
        stressan=stressab-stressat;
        strainap=calstrain(stressap,xm,xh,Eh,Eu,fy,fu);
        strainan=calstrain(stressan,xm,xh,Eh,Eu,fy,fu);
        rfa1=(1-(strainap-strainan)/2);
        if abs(rfa1-rfatry)/rfatry<0.01
            rfa=rfatry;
  
        end
    end
    for thetactry=thetar:0.001:(xu-xh)*2
        thetad=thetactry-thetar;
        mc1=calmoment(thetactry,rfc,my,xm,xh,Eh,E,xy,Eu);
        md=calmoment(thetad,rfd,my,xm,xh,Eh,E,xy,Eu);
        dertac=g1*sin(thetad);
        for dertapcatry=2-elastic:0.0001:2-elastic
            thetaa=acos(1-(dertac+dertapcatry-z*(1-cos(thetar)+sin(thetar)))/g2);
            dertab=g2*sin(thetaa);
            thetab=thetaa-thetar;
            ma=calmoment(thetaa,rfa,my,xm,xh,Eh,E,xy,Eu);
            mb=calmoment(thetab,rfb,my,xm,xh,Eh,E,xy,Eu);
            f=((mc1+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
            n1=((mc1+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);

            if dertapcatry<0.000000001
                if n1<fff
                    dertapca=0;
                end
            else
                if dertapcatry<2-elastic
                    abs(n1-fff)/fff;
                    if abs(n1-fff)/fff<0.01
                        dertapca=dertapcatry;
                    end
                else
                    if n1>fff
                        dertapca=2-elastic;
                    end
                end
            end
        end
        thetaa=acos(1-(dertac+dertapca-z*(1-cos(thetar)+sin(thetar)))/g2);
        thetab=thetaa-thetar;
        ma=calmoment(thetaa,rfa,my,xm,xh,Eh,E,xy,Eu);
        mb=calmoment(thetab,rfb,my,xm,xh,Eh,E,xy,Eu);
        dertab=g2*sin(thetaa);
        f=((mc1+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
        n=((mc1+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
        mc2=mb+f*z*(cos(thetar)-sin(thetar))-n*z*(cos(thetar)+sin(thetar));
        abs(mc2-mc1)/mc1;
        if abs(mc2-mc1)/mc2<0.01
            if abs(mc2-mc1)/mc1<0.01
                thetactry;
                if cd5<1.1
                    thetactry
                    thetac=thetactry;
                    cd5=2;
                end
            end
        end
    end       
    thetad=thetac-thetar;
    mc=calmoment(thetac,rfc,my,xm,xh,Eh,E,xy,Eu);
    md=calmoment(thetad,rfd,my,xm,xh,Eh,E,xy,Eu);
    dertac=g1*sin(thetad);
    for dertapcatry=2-elastic:0.0001:2-elastic
        thetaa=acos(1-(dertac+dertapcatry-z*(1-cos(thetar)+sin(thetar)))/g2);
        dertab=g2*sin(thetaa);
        thetab=thetaa-thetar;
        ma=calmoment(thetaa,rfa,my,xm,xh,Eh,E,xy,Eu);
        mb=calmoment(thetab,rfb,my,xm,xh,Eh,E,xy,Eu);
        f=((mc+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
        n1=((mc+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
        if dertapcatry<0.000000001
            if n1<fff
                dertapca=0;
            end
        else
            if dertapcatry<2-elastic
                abs(n1-fff)/fff;
                if abs(n1-fff)/fff<0.01
                    dertapca=dertapcatry;
                end
            else
                if n1>fff
                    dertapca=2-elastic;
                end
            end
        end
    end
    thetaa=acos(1-(dertac+dertapca-z*(1-cos(thetar)+sin(thetar)))/g2);
    thetab=thetaa-thetar;
    ma=calmoment(thetaa,rfa,my,xm,xh,Eh,E,xy,Eu);
    mb=calmoment(thetab,rfb,my,xm,xh,Eh,E,xy,Eu);
    dertab=g2*sin(thetaa);
    f=((mc+md)*dertab+(ma+mb)*g2*cos(thetad))/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);
    n=((mc+md)*g2*cos(thetaa)+(ma+mb)*dertac)/(g1*g2*cos(thetad)*cos(thetaa)-dertab*dertac);

    derta=dertab+z*(sin(thetar)+cos(thetar)-1)+g2*(cos(thetad)-1);
    L(a)=f/1000;
    f/1000
    dert(a)=derta;
    derta

    stressba=calstress(thetaa,xm,xh,Eh,Eu,fy,fu);
    stressbb=calstress(thetab,xm,xh,Eh,Eu,fy,fu);
    stressbc=calstress(thetac,xm,xh,Eh,Eu,fy,fu);
    stressbd=calstress(thetad,xm,xh,Eh,Eu,fy,fu);
    stressta=n/b/t/rfa;
    stressal=stressta+stressba;
    if stressal>ff
        L(a)=0;
    end
    stresstb=n/b/t/rfb;
    stressbl=stresstb+stressbb;
    if stressbl>ff
        L(a)=0;
    end  
    stresstc=f/b/t/rfc;
    stresscl=stresstc+stressbc;
    if stresscl>ff
        L(a)=0;
    end               
                  
    stresstd=f/b/t/rfd;
    stressdl=stresstd+stressbd;
    if stressdl>ff
        L(a)=0;
    end             
    a=a+1;
end
plot(dert,L);

function moment =calmoment(theta,thickness,my,xm,xh,Eh,E,xy,Eu)
    strain=theta/2;
    if theta<xm*2
        moment=0.5*my*thickness^2*(3-(xy/(strain+xh))^2)+my*thickness^2/2*Eh/E*(strain)/xy*(1-xh/(strain+xh))*(2+xh/(strain+xh));
    else
        moment=0.5*my*thickness^2*(3-(xy/(strain+xh))^2)+my*thickness^2/2*Eh/E*(strain)/xy*(1-xh/(strain+xh))*(2+xh/(strain+xh))-my*thickness^2/2*(Eh-Eu)/E*((strain+xh)-xm)/xy*(1-xm/(strain+xh))*(2+xm/(strain+xh));
    end
end

function stress =calstress(theta,xm,xh,Eh,Eu,fy,fu)
    strain=theta/2;
    if strain<xm
        stress=fy+(strain)*Eh;
    else
        stress=fu+(strain+xh-xm)*Eu;
    end
end

function strain =calstrain(stress,xm,xh,Eh,Eu,fy,fu)
    if stress<fu
        strain=xh+(stress-fy)/Eh;
    else
        strain=xm+(stress-fu)/Eu;
    end
end

