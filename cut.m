## Created: 2022-07-01
% cut function
function retval = cut (x, y)
 % x=1.5, y = infsup(-1,1)
if inIR(infsup(x,x), y) > 0
  retval =x;
end  
if x > sup(y) 
  retval = sup(y);
end  
if x < inf(y) 
  retval = inf(y);
end 
end
