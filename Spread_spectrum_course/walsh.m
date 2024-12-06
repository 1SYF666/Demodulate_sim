%编写Walsh码函数
function walsh=walsh(N)
M=ceil(log2(N));   %ceil(a):“天花板”，取比a大的最小整数，即朝正无穷方向取整。
wn=0;
for i=1:M
    w2n=[wn,wn;wn,~wn];
    wn=w2n;
end
walsh=wn;