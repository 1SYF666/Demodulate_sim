%��дWalsh�뺯��
function walsh=walsh(N)
M=ceil(log2(N));   %ceil(a):���컨�塱��ȡ��a�����С�����������������ȡ����
wn=0;
for i=1:M
    w2n=[wn,wn;wn,~wn];
    wn=w2n;
end
walsh=wn;