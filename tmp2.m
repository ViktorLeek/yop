figure(1);
subplot(311); hold on
sol.plot(t, x, 'mag', 10);

subplot(312); hold on
sol.plot(t, v);
td = sol.value(int(abs(v)));
text(0.3, 0.5, ['Traveled distance is ', num2str(td)], 'FontSize', 14)

subplot(313); hold on
sol.stairs(t, a);
J_min = sol.value(0.5*int(a^2));
text(0.35, -2, ['Minimum cost ', num2str(J_min)], 'FontSize', 14)