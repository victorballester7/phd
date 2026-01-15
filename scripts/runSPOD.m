u_mean = mean(u, 1);
v_mean = mean(v, 1);

u_fluct = u - u_mean;
v_fluct = v - v_mean;

Z = [u_fluct, v_fluct];

point_example = 480;
disp("Plotting v of point at x = " + num2str(x(point_example)) + ", y = " + num2str(y(point_example)));
figure;
plot(time, v_fluct(:, point_example));
xlabel('Time');
ylabel('v fluctuation');
title('v fluctuation at example point');


%%%%%% doing spod
[L, P, F] = spod(Z, 2500, w);

[value, indx] = max(L(:,1));

[sortedLVals, sortedLIdx] = sort(L(:,1), 'descend');

% plot dominant frequencies with their energy
disp('Dominant frequencies and their energies:');
for i = 1:7
    disp(['Frequency (f) = ', num2str(F(sortedLIdx(i))), ' with energy = ', num2str(L(sortedLIdx(i),1))]);
end
figure;
plot(F, L, 'o-');

% plotting main mode
tri = delaunay(x,y);    
mainModeRe = real(P(sortedLIdx(1), :, 1));
mainModeRe_u = mainMode_re(1:length(mainMode_re)/2);
mainModeRe_v = mainMode_re(length(mainMode_re)/2+1:end);
figure;
trisurf(tri, x, y, mainModeRe_u); shading interp;
xlabel('X');
ylabel('Y');
title('SPOD Main Mode - u component');
%change default view
view(2);
figure;
trisurf(tri, x, y, mainModeRe_v); shading interp;
xlabel('X');
ylabel('Y');
title('SPOD Main Mode - v component');
view(2);


