function k = wrapper_spatialStroop_CE_SOA(n, w)
% w(1)-w(4)-A_0
% w(5)-w(8)-kc
% w(9)-w(12)-ki
% w(13)-w(16)-da
% w(17)-w(20)-dc
% w(21)-w(24)-a
% w(25)-w(28)-Ter
% w(29)-std
% w(30)-w(33)-c
    if n==1 % null model
       k = spatialStroop_CE_SOA(w(1),w(1),w(1),w(1),w(2),w(2),w(2),w(2),w(3),w(3),w(3),w(3),w(4),w(4),w(4),w(4),w(5),w(5),w(5),w(5),...
            w(6),w(6),w(6),w(6),w(7),w(7),w(7),w(7),w(8),0,0,0,0);  
    elseif n==2 % Ter
        k = spatialStroop_CE_SOA(w(1),w(1),w(1),w(1),w(2),w(2),w(2),w(2),w(3),w(3),w(3),w(3),w(4),w(4),w(4),w(4),w(5),w(5),w(5),w(5),...
            w(6),w(6),w(6),w(6),w(7),w(8),w(9),w(10),w(11),0,0,0,0);
    elseif n==3 % c, Ter
        k = spatialStroop_CE_SOA(w(1),w(1),w(1),w(1),w(2),w(2),w(2),w(2),w(3),w(3),w(3),w(3),w(4),w(4),w(4),w(4),w(5),w(5),w(5),w(5),...
            w(6),w(6),w(6),w(6),w(7),w(8),w(9),w(10),w(11),0,0,w(12),w(12));
    elseif n==4 % da, Ter
        k = spatialStroop_CE_SOA(w(1),w(1),w(1),w(1),w(2),w(2),w(2),w(2),w(3),w(3),w(3),w(3),w(4),w(5),w(6),w(7),w(8),w(8),w(8),w(8),...
            w(9),w(9),w(9),w(9),w(10),w(11),w(12),w(13),w(14),0,0,0,0);
    elseif n==5 % A0, Ter
        k = spatialStroop_CE_SOA(w(1),w(2),w(3),w(4),w(5),w(5),w(5),w(5),w(6),w(6),w(6),w(6),w(7),w(7),w(7),w(7),w(8),w(8),w(8),w(8),...
            w(9),w(9),w(9),w(9),w(10),w(11),w(12),w(13),w(14),0,0,0,0);
    elseif n==6 % da (-100), Ter
        k = spatialStroop_CE_SOA(w(1),w(1),w(1),w(1),w(2),w(2),w(2),w(2),w(3),w(3),w(3),w(3),w(4),w(4),w(4),w(5),w(6),w(6),w(6),w(6),...
            w(7),w(7),w(7),w(7),w(8),w(9),w(10),w(11),w(12),0,0,0,0);
    end
end