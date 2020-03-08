function dy = pendulumCart(Y,u,d,A, B)
    %slo by to i jako dx = ax +Bu a pak dx(3)+=d(1) atd..  
    dy = [        
        Y(3); %Toto už je derivace x
        Y(4); %Toto už je derivace alfa
        d(1)+ A(3,1)*Y(1)+ A(3,2)*Y(2)+ A(3,3)*Y(3)+ A(3,4)*Y(4)+ B(3,1)*u;
        d(2) + A(4,1)*Y(1)+ A(4,2)*Y(2)+ A(4,3)*Y(3)+ A(4,4)*Y(4)+ B(3,1)*u;
        ];
end