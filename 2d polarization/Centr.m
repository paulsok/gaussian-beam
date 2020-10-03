function [centr]= Centr(E0) %incident beam center of gravity calculation
        [N M]=size(E0);
        full=0;
        centrmass(1:2)=0;
        for i=1:N
            for j=1:N
                centrmass(1)=centrmass(1)+abs(E0(i,j))^2*i;
                centrmass(2)=centrmass(2)+abs(E0(i,j))^2*j;
                full=full+abs(E0(i,j))^2;
            end;
        end;
        centr(1)=centrmass(1)/full;
        centr(2)=centrmass(2)/full;
