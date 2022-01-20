figure(2); 
subplot(511); hold on
sol.plot(t/3600, x(1))
xlabel('Time [h]'); ylabel('Altitude [m]')
subplot(512); hold on
sol.plot(t/3600, x(4)*3.6)
xlabel('Time [h]'); ylabel('True airspeed [km/h]')
subplot(513); hold on
sol.plot(t/3600, r2d(x(5)))
xlabel('Time [h]'); ylabel('Flight path angle [deg]')
subplot(514); hold on
sol.plot(t/3600, r2d(x(6)))
xlabel('Time [h]'); ylabel('Tracking angle [deg]')
subplot(515); hold on
sol.plot(t/3600, x(7))
xlabel('Time [h]'); ylabel('mass [kg]')

figure(3); 
subplot(311); hold on
sol.plot(t/3600, r2d(u(1)))
xlabel('Time [h]'); ylabel('Angle of attack [deg]')
subplot(312); hold on
sol.plot(t/3600, r2d(u(2)))
xlabel('Time [h]'); ylabel('Roll angle [deg]')
subplot(313); hold on
sol.stairs(t/3600, u(3)*100)
xlabel('Time [h]'); ylabel('Throttle setting [%]')