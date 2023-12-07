
%{
CLASS: HUG
    Hugoniot with linear Us-Up.

    INPUTS:
        r0: double
            Initial density

        C:  double
            Hugoniot Us-up intercept

        s: double
            Hugoniot Us-up slope

    PROPERTIES:
        P: Lambda(up)
            Pressure as a function of particle velocity (up)

        Us: lambda(up)
            Shock Velocity as a function of up

        v: lambda(up)
            spcific volume as a function of up            

%}
classdef Hug2
    properties
        P;Us;v;
        C1;s1;C2;s2;r0;
        name;color;
    end
    methods
        function self=Hug2(r0,C1,s1,C2,s2,name,color)
            self.C1=C1;self.C2=C2;
            self.s1=s1;self.s2=s2;
            self.r0=r0;
            self.name=name;
            self.color=color;
            self.Us=@(u) min([C1+s1*u,C2+s2*u]);
            self.P=@(u) r0*self.Us(u).*u.*(u>0)+r0*C1*u.*(u<=0);
            self.v=@(u) (1-u/self.Us(u))/r0;
        end
    end
end