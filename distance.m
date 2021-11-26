function d = distance(F, f1,f2)


f2t_F_f1= f2'*F*f1;

d = abs(f2t_F_f1);


end