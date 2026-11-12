function k = wrapper_flanker_CE_SOA(n, w)
% w(1)-w(4)-A_0
% w(5)-w(8)-kc
% w(9)-w(12)-ki
% w(13)-w(16)-da
% w(17)-w(20)-dc
% w(21)-w(24)-a
% w(25)-w(28)-Ter
% w(29)-std
% w(30)-w(33)-c
    if  n==1 % A0, kc, ki, da, dc, Ter
        k = flanker_CE_SOA(w(1),w(2),w(3),w(4),w(5),w(6),w(7),w(8),w(9),w(10),w(11),w(12),w(13),w(14),w(15),w(16),w(17),w(18),w(19),...
            w(20),w(21),w(21),w(21),w(21),w(22),w(23),w(24),w(25),w(26),0,0,0,0);
    elseif n==2 % A0, kc, ki, da, dc, a, Ter
        k = flanker_CE_SOA(w(1),w(2),w(3),w(4),w(5),w(6),w(7),w(8),w(9),w(10),w(11),w(12),w(13),w(14),w(15),w(16),w(17),w(18),w(19),...
            w(20),w(21),w(22),w(23),w(24),w(25),w(26),w(27),w(28),w(29),0,0,0,0);
    elseif n==3   % A0, kc, ki, da, dc, c, Ter
        k = flanker_CE_SOA(w(1),w(2),w(3),w(4),w(5),w(6),w(7),w(8),w(9),w(10),w(11),w(12),w(13),w(14),w(15),w(16),w(17),w(18),w(19),...
            w(20),w(21),w(21),w(21),w(21),w(22),w(23),w(24),w(25),w(26),w(27),0,w(27),w(27));
    elseif n==4 % A0, kc, ki, da, dc, c (negative SOAs), Ter
        k = flanker_CE_SOA(w(1),w(2),w(3),w(4),w(5),w(6),w(7),w(8),w(9),w(10),w(11),w(12),w(13),w(14),w(15),w(16),w(17),w(18),w(19),...
            w(20),w(21),w(21),w(21),w(21),w(22),w(23),w(24),w(25),w(26),0,0,w(27),w(27));
    end
end