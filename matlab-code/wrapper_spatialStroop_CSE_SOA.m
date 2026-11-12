function k = wrapper_spatialStroop_CSE_SOA(n, w)
% w(1)-w(8)-A_0
% w(9)-w(10)-kc
% w(11)-w(12)-ki
% w(13)-w(20)-da
% w(21)-dc
% w(22)-w(23)-a
% w(24)-w(27)-Ter
% w(28)-std
% w(29)-c
    if n==1 % null model
        k = spatialStroop_CSE_SOA(w(1),w(1),w(1),w(1),w(1),w(1),w(1),w(1),w(2),w(2),w(3),w(3),w(4),w(4),...
            w(4),w(4),w(4),w(4),w(4),w(4),w(5),w(6),w(6),w(7),w(8),w(9),w(10),w(11),0);
    elseif n==2 % restricted model -- A0
        k = spatialStroop_CSE_SOA(w(1),w(2),w(2),w(2),w(1),w(2),w(2),w(2),w(3),w(3),w(4),w(4),w(5),w(5),...
            w(5),w(5),w(5),w(5),w(5),w(5),w(6),w(7),w(7),w(8),w(9),w(10),w(11),w(12),0);
    elseif n==3 % restricted model -- A0, da (-100 ms)
        k = spatialStroop_CSE_SOA(w(1),w(2),w(2),w(2),w(1),w(2),w(2),w(2),w(3),w(3),w(4),w(4),w(5),w(5),...
            w(5),w(6),w(5),w(5),w(5),w(6),w(7),w(8),w(8),w(9),w(10),w(11),w(12),w(13),0);
    elseif n==4 % restricted model -- A0, da (-50-ms, -100-ms previ)
        k = spatialStroop_CSE_SOA(w(1),w(2),w(2),w(2),w(1),w(2),w(2),w(2),w(3),w(3),w(4),w(4),w(5),w(5),...
            w(5),w(5),w(6),w(6),w(7),w(7),w(8),w(9),w(9),w(10),w(11),w(12),w(13),w(14),0);
    elseif n==5 % restricted model -- A0, da (-100-ms & -50-ms, -100-ms previ)
        k = spatialStroop_CSE_SOA(w(1),w(2),w(2),w(2),w(1),w(2),w(2),w(2),w(3),w(3),w(4),w(4),w(5),w(5),...
            w(5),w(6),w(6),w(6),w(7),w(8),w(9),w(10),w(10),w(11),w(12),w(13),w(14),w(15),0);
    elseif n==6 % restricted model -- A0, da (previ)
        k = spatialStroop_CSE_SOA(w(1),w(2),w(2),w(2),w(1),w(2),w(2),w(2),w(3),w(3),w(4),w(4),w(5),w(5),...
            w(5),w(6),w(7),w(7),w(8),w(8),w(9),w(10),w(10),w(11),w(12),w(13),w(14),w(15),0);
    elseif n==7 % restricted model -- A0, da (-100 ms & previ)
        k = spatialStroop_CSE_SOA(w(1),w(2),w(2),w(2),w(1),w(2),w(2),w(2),w(3),w(3),w(4),w(4),w(5),w(5),...
            w(5),w(6),w(7),w(7),w(8),w(9),w(10),w(11),w(11),w(12),w(13),w(14),w(15),w(16),0);
    elseif n==8 % da 
        k = spatialStroop_CSE_SOA(w(1),w(2),w(2),w(2),w(1),w(2),w(2),w(2),w(3),w(3),w(4),w(4),w(5),w(6),...
            w(7),w(8),w(9),w(10),w(11),w(12),w(13),w(14),w(14),w(15),w(16),w(17),w(18),w(19),0);
    end
end