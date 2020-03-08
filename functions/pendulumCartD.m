function Yk1 = pendulumCartD(Yk0,uk,dk,Ts,p)
    M = 2;
    delta = Ts/M;
    Yk1 = Yk0;
    for ct=1:M
        Yk1 = Yk1 + delta*pendulumCart(Yk1,uk,dk,p);
    end