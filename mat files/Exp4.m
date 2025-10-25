% Gauss_Seidel : LOAD FLOW PROGRAM

% Line data: lineno, from bus, to bus, resistance, reactance, line charging, tap ratio
data = [
    1 1 2 0.02 0.06 0.030 1.0
    2 1 3 0.08 0.24 0.025 1.0
    3 2 3 0.06 0.18 0.020 1.0
    4 2 4 0.06 0.18 0.020 1.0
    5 2 5 0.04 0.12 0.015 1.0
    6 3 4 0.01 0.03 0.010 1.0
    7 4 5 0.08 0.24 0.025 1.0
];

% Bus data: PGEN, QGEN, PLOAD, QLOAD, Vsp, Angle, BusCode
% BusCode: 3=Slack, 2=PVBUS, 1=PQBUS
bus = [
    0.00 0.00 0.00 0.00 1.06 0.00 3
    0.40 0.30 0.20 0.10 1.00 0.00 1
    0.00 0.00 0.45 0.15 1.00 0.00 1
    0.00 0.00 0.40 0.05 1.00 0.00 1
    0.00 0.00 0.60 0.10 1.00 0.00 1
];

% Extract line data
ln = data(:,1);
ns = data(:,2);
nf = data(:,3);
r = data(:,4);
x = data(:,5);
ych = data(:,6);
tap = data(:,7);

nline = length(data(:,2));
n = max(max(ns), max(nf));

% Extract bus data
pg = bus(:,1);
qg = bus(:,2);
pl = bus(:,3);   % Active power load (column vector)
ql = bus(:,4);   % Reactive power load (column vector)
es = bus(:,5);
an = bus(:,6);
code = bus(:,7);

eps = 1e-4;    % Tolerance for convergence
maxit = 30;      % Maximum number of iterations

% Formation of bus admittance matrix Ybus
yb = zeros(n, n); % Initialize Ybus matrix

sz = zeros(nline,1); wz = zeros(nline,1); % initialize optional arrays
for p = 1:nline
    z(p) = r(p) + 1j*x(p);
    ybr(p) = 1 / z(p);
    n1 = ns(p);
    n2 = nf(p);
    ys(p) = 1j * ych(p);

    if tap(p) ~= 1
        t = 1 / tap(p);
        sz(p) = sz(p) + ybr(p) * (t - 1) * t;
        wz(p) = wz(p) + ybr(p) * (1 - t);
        ybr(p) = ybr(p) * t;
    end

    yb(n1, n1) = yb(n1, n1) + ybr(p) + sz(p) + ys(p);
    yb(n2, n2) = yb(n2, n2) + ybr(p) + wz(p) + ys(p);
    yb(n1, n2) = yb(n1, n2) - ybr(p);
    yb(n2, n1) = yb(n2, n1) - ybr(p);
end

% Initialize voltages and net injections as column vectors
ps = zeros(n,1); qs = zeros(n,1);
vo = zeros(n,1); vn = zeros(n,1);

for p = 1:n
    ps(p) = pg(p) - pl(p);
    qs(p) = qg(p) - ql(p);
    vo(p) = es(p)*cos(an(p)) + 1j*es(p)*sin(an(p)); % column
    vn(p) = vo(p);
end

% Iterative Gauss-Seidel load flow
for icount = 1:maxit
    for p = 2:n
        a = (ps(p) - 1j*qs(p)) / conj(vo(p));
        sum1 = 0 + 1j*0;
        sum2 = 0 + 1j*0;

        for q = 1:p-1
            sum1 = sum1 + yb(p,q)*vn(q);
        end
        for q = p+1:n
            sum2 = sum2 + yb(p,q)*vo(q);
        end

        vn(p) = (a - sum1 - sum2) / yb(p,p);

        % Acceleration factor
        vn(p) = vo(p) + 1.4 * (vn(p) - vo(p));

        if code(p) == 2 % PVBUS: voltage magnitude fixed
            delta = angle(vn(p));
            vn(p) = es(p)*cos(delta) + 1j*es(p)*sin(delta);
        end
    end

    dv = vn - vo;
    dvm = max(abs(dv));
    if dvm <= eps
        break;
    end

    vo = vn;
end

% Check convergence
if icount < maxit
    fprintf('\nLoad flow converged in %d iterations using GAUSS-SEIDEL method.\n', icount);
else
    fprintf('\nLoad flow did NOT converge in %d iterations using GAUSS-SEIDEL method.\n', maxit);
end

% Print bus voltages and angles
disp('BUSNO  BUSMAG    BUSANGLE (degrees)');
con = 180/pi;
for p = 1:n
    fprintf('%4d %8.4f %12.4f\n', p, abs(vo(p)), angle(vo(p))*con);
end

% Calculate bus power injections after load flow
pcal = zeros(n,1);
qcal = zeros(n,1);
for p = 1:n
    sum1 = 0 + 1j*0;
    for q = 1:n
        sum1 = sum1 + yb(p,q)*vo(q);
    end
    s = conj(vo(p)) * sum1;
    pcal(p) = real(s);
    qcal(p) = -imag(s);
end

% Print power generation, load, voltages and angles
fprintf('\nLOAD_FLOW SOLUTION (All quantities except angles in P.U.)\n');
fprintf(' BUS  PG       QG       PL       QL       VOLT     ANG(deg)\n');
for p = 1:n
    fprintf('%4d %8.3f %8.3f %8.3f %8.3f %8.4f %10.4f\n', ...
        p, pcal(p)+pl(p), qcal(p)+ql(p), pl(p), ql(p), abs(vo(p)), angle(vo(p))*con);
end

fprintf('\n');
fprintf('Total Real Power Generation = %8.4f P.U.\n', sum(pcal + pl));
fprintf('Total Reactive Power Generation = %8.4f P.U.\n', sum(qcal + ql));
fprintf('Total Real Power Load = %8.4f P.U.\n', sum(pl));
fprintf('Total Reactive Power Load = %8.4f P.U.\n', sum(ql));

% Calculate and print line flows
fprintf('\nLINE FLOWS IN THE SYSTEM (P.U.)\n');
fprintf(' From To   Real Power   Reactive Power\n');

total_loss = 0 + 1j*0;

for p = 1:nline
    n1 = ns(p);
    n2 = nf(p);
    zline = r(p) + 1j*x(p);
    sh1 = 1j * ych(p);

    temp1 = (vo(n1) - vo(n2))/zline + vo(n1)*sh1;
    temp2 = (vo(n2) - vo(n1))/zline + vo(n2)*sh1;

    sp1 = vo(n1) * conj(temp1);
    sp2 = vo(n2) * conj(temp2);

    total_loss = total_loss + sp1 + sp2;

    fprintf('%4d %4d %12.3f %12.3f\n', n1, n2, real(sp1), imag(sp1));
    fprintf('%4d %4d %12.3f %12.3f\n\n', n2, n1, real(sp2), imag(sp2));
end

fprintf('Total system loss:\n');
fprintf('Real Power Loss = %8.4f P.U.\n', real(total_loss));
fprintf('Reactive Power Loss = %8.4f P.U.\n', imag(total_loss));
