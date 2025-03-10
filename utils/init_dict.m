function D = init_dict(para)
M = para.m; K = para.K;
s = max(M);
D =[];

for i = 1:length(M)
    d = randn(M(i),M(i),K(i));
    d  = (d)./sqrt(sum(d.^2,1:2));
    d = padarray(d,[s-M(i)  s-M(i)],'post');
    D = cat(4,D,d);
end

end