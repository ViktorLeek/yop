for k=1:10
figure(1)
subplot(411); hold on
sol.plot(t, w_ice);
subplot(412); hold on
sol.plot(t, p_im);
subplot(413); hold on
sol.plot(t, p_em);
subplot(414); hold on
sol.plot(t, w_tc);

figure(2)
subplot(311); hold on
sol.stairs(t, u_f)
sol.plot(t, y.u_f_max);
subplot(312); hold on
sol.stairs(t, u_wg)
subplot(313); hold on
sol.stairs(t, P_gen)

sol.value(int(y.cylinder.fuel_massflow));
end