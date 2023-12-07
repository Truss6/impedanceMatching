
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
classdef Hugf
    properties
        P;Us;v;Ur;e;
        r0;
        name;color;
    end
    methods
        function self=Hugf(r0,Us,name,color)
            self.r0=r0;
            self.name=name;
            self.color=color;
            if color==false
                self.color=[rand rand rand];
            end
            self.Us=Us;
            self.P=@(u) r0*Us(u).*u*(u>0)+r0*Us(0)*u*(u<=0);
            self.v=@(u) (1-u/self.Us(u))/r0;
            v0=1/r0;
            self.e=@(u) (1-u/self.Us(u));
%             self.Ur=@(u) sqrt(((self.v(u).^2)*(C^2).*(v0+s*(v0-self.v(u))))./((v0-s*(v0-self.v(u))).^3));
        end
    end
end