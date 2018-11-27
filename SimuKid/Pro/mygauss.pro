function  mygauss, x, p

return, p[0]+p[1]*exp(-(x-p[3])^2/(2*p[2]^2)) + x*p[4]

end
