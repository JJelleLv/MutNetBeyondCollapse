/* all Jacobians and higher derivatives for nomal models
   runwhich makes it possible to run part of the Jacobians (no check for dependency)
   vector with 7 elements 1=yes 0=no  order: Jacobian, Jacobianp Hessian Hessianp der3 der4 der5*/
getJacobi# := proc(maxtime#,runwhich#)
    local starttime#, neqs#, nvars#, npars#, i1#, i2#, i3#, i4#, i5#, i6#;
begin
    starttime#:=rtime();
    neqs#:=nops(eqns#);
    nvars#:=nops(vars#);
    npars#:=nops(pars#);
    Jacobian#:=NIL;
    Jacobianp#:=NIL;
    Hessian#:=NIL;
    Hessianp#:=NIL;
    der3#:=NIL;
    der4#:=NIL;
    der5#:=NIL;

    //Jacobian
    if runwhich#[1]=1 then
        Jacobian#:=array(1..neqs#,1..nvars#);
        for i1# from 1 to neqs# do
            for i2# from 1 to nvars# do
                Jacobian#[i1#,i2#]:=diff(eqns#[i1#],vars#[i2#]);
            end_for;
        end_for;
        if rtime()-starttime#>maxtime# then
            return(0);
        end_if;
    end_if;

    //Jacobianp
    if runwhich#[2]=1 then
        Jacobianp#:=array(1..neqs#,1..npars#);
        for i1# from 1 to neqs# do
            for i2# from 1 to npars# do
                Jacobianp#[i1#,i2#]:=diff(eqns#[i1#],pars#[i2#]);
            end_for;
        end_for;
        if rtime()-starttime#>maxtime# then
            return(0);
        end_if;
    end_if;

    //Hessian - skips symmetrical part
    if runwhich#[3]=1 then
        Hessian#:=array(1..neqs#,1..nvars#,1..nvars#);
        for i1# from 1 to neqs# do
            for i2# from 1 to nvars# do
                for i3# from i2# to nvars# do
                    Hessian#[i1#,i2#,i3#]:=diff(Jacobian#[i1#,i2#],vars#[i3#]);
                end_for;
            end_for;
        end_for;
        if rtime()-starttime#>maxtime# then
            return(0);
        end_if;
    end_if;

    //Hessianp
    if runwhich#[4]=1 then
        Hessianp#:=array(1..neqs#,1..nvars#,1..npars#);
        for i1# from 1 to neqs# do
            for i2# from 1 to nvars# do
                for i3# from 1 to npars# do
                    Hessianp#[i1#,i2#,i3#]:=diff(Jacobian#[i1#,i2#],pars#[i3#]);
                end_for;
            end_for;
        end_for;
        if rtime()-starttime#>maxtime# then
            return(0);
        end_if;
    end_if;

    //der3 - skips symmetrical part
    if runwhich#[5]=1 then
        der3#:=array(1..neqs#,1..nvars#,1..nvars#,1..nvars#);
        for i1# from 1 to neqs# do
            for i2# from 1 to nvars# do
                for i3# from i2# to nvars# do
                    for i4# from i3# to nvars# do
                        der3#[i1#,i2#,i3#,i4#]:=diff(Hessian#[i1#,i2#,i3#],vars#[i4#]);
                    end_for;
                end_for;
            end_for;
        end_for;
        if rtime()-starttime#>maxtime# then
            return(0);
        end_if;
    end_if;

    //der4 - skips symmetrical part
    if runwhich#[6]=1 then
        der4#:=array(1..neqs#,1..nvars#,1..nvars#,1..nvars#,1..nvars#);
        for i1# from 1 to neqs# do
            for i2# from 1 to nvars# do
                for i3# from i2# to nvars# do
                    for i4# from i3# to nvars# do
                        for i5# from i4# to nvars# do
                            der4#[i1#,i2#,i3#,i4#,i5#]:=diff(der3#[i1#,i2#,i3#,i4#],vars#[i5#]);
                        end_for;
                    end_for;
                end_for;
            end_for;
        end_for;
        if rtime()-starttime#>maxtime# then
            return(0);
        end_if;
    end_if;

    //der5 - skips symmetrical part
    if runwhich#[7]=1 then
        der5#:=array(1..neqs#,1..nvars#,1..nvars#,1..nvars#,1..nvars#,1..nvars#);
        for i1# from 1 to neqs# do
            for i2# from 1 to nvars# do
                for i3# from i2# to nvars# do
                    for i4# from i3# to nvars# do
                        for i5# from i4# to nvars# do
                            for i6# from i5# to nvars# do
                                der5#[i1#,i2#,i3#,i4#,i5#,i6#]:=diff(der4#[i1#,i2#,i3#,i4#,i5#],vars#[i6#]);
                            end_for;
                        end_for;
                    end_for;
                end_for;
            end_for;
        end_for;
    end_if;
    return(1)
end_proc;

/*implicitJacs(dvars#) makes 3 Jacobians for implicit dae models
 we can use the general formula for implicit differentiation of a
 function R(x,y)=0

 dy/dx= -(dR/dx)/(dR/dy)
 comes from the generalized chain rule
 0 = R(x,y)
 d 0/dx = d(R(x,y)/dx
 0= d(R(x,y)/dx

 Another method is to calculate the implicit differntial and solve unknowns:
 R:=subs(R,y=y(x));
 dRdx :=diff(R,x);
 dydt:=solve(dRdx=0,diff(y(x),y))
 same result as above*/

implicitJacs# := proc(dvars#)
    local neqs#, i1#, i2#, i3#;
begin
    neqs#:=nops(eqns#);
    Jacobian#:=array(1..neqs#,1..neqs#);
    Jacobian_y#:=array(1..neqs#,1..neqs#);
    Jacobian_yp#:=array(1..neqs#,1..neqs#);
    for i1# from 1 to neqs# do
        for i2# from 1 to neqs# do
            Jacobian_y#[i1#,i2#]:=diff(eqns#[i1#],vars#[i2#]);
            Jacobian_yp#[i1#,i2#]:=diff(eqns#[i1#],dvars#[i2#]);
        end_for;
    end_for;
    for i1# from 1 to neqs# do
        for i2# from 1 to neqs# do
            // %note that we need to take the diagonal of Jac_yp [dx/dx' dy/dx';dxd/dy' dy/dy']
            Jacobian#[i1#,i2#]:=-Jacobian_y#[i1#,i2#]/Jacobian_yp#[i1#,i1#]
        end_for;
    end_for;
    return(1);
end_proc:

Sensitivities# := proc()
begin
    neqs#:=nops(eqns#);
    neqs#:=nops(pars#);
    Sensitivp#:=array(1..neqs#,1..neqs#);
    for i2# from 1 to npars# do
        eqns1#=eqns#;
        for i1# from 1 to neqs# do
            subs(eqns1#,vars#[i1]=vars#[i1](pars#[i2]));
        end_for;
        for i1# from 1 to neqs# do
            Sensitivp#[i1#,i2#]:=diff(eqns1#[i1#],pars#[i2#]);
        end_for;
    end_for;
    return(1);
end_proc


