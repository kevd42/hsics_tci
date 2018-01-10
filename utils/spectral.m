function s = spectral(m)

if nargin < 1, m = size(get(gcf,'colormap'),1); end

base = [
0.4667 0.0000 0.5333;
0.5333 0.0000 0.6000;
0.0000 0.0000 0.6667;
0.0000 0.0000 0.8667;
0.0000 0.4667 0.8667;
0.0000 0.6000 0.8667;
0.0000 0.6667 0.6667;
0.0000 0.6667 0.5333;
0.0000 0.6000 0.0000;
0.0000 0.7333 0.0000;
0.0000 0.8667 0.0000;
0.0000 1.0000 0.0000;
0.7333 1.0000 0.0000;
0.9333 0.9333 0.0000;
1.0000 0.8000 0.0000;
1.0000 0.6000 0.0000;
1.0000 0.0000 0.0000;
0.8667 0.0000 0.0000;
0.8000 0.0000 0.0000;
];


s = zeros(m,3);

for i = 1:3
    s(:,i) = linterp(linspace(1,m,19),base(:,i),1:m);
end
