m1=2.3;
m2=2.1;
d1=3.2;
d2=8.6;
k=400;
km=2.95;
ky1=280;
ky2=280;

if exist('m1_m', 'var') == 0
	m1_m = m1;
	m2_m = m2;
	d1_m = d1;
	d2_m = d2;
end

res=inputdlg({
	['m1 (default ' num2str(m1) ')'],
	['m2 (default ' num2str(m2) ')'],
	['d1 (default ' num2str(d1) ')'],
	['d2 (default ' num2str(d2) ')']}, 'Process parameters:', 1, {
	num2str(m1_m),
	num2str(m2_m),
	num2str(d1_m),
	num2str(d2_m)});

if length(res) ~= 0
	m1_m=eval(res{1});
	m2_m=eval(res{2});
	d1_m=eval(res{3});
	d2_m=eval(res{4});

	A=[0 1 0 0;-k/m1_m -d1_m/m1_m k/m1_m 0; 0 0 0 1; k/m2_m 0 -k/m2_m -d2_m/m2_m];
	B=[0 km/m1_m 0 0]';
	C_proc=[ky1 0 0 0; 0 0 ky2 0];
	C1=C_proc(1,:);
	C2=C_proc(2,:);

	% process model on state-space form
	P_sim=ss(A,B,C2,0);
end
