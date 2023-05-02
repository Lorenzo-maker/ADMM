function updatePlot_localPOE(Glocal, q, transforms)
%
%
%

Glocalnum = full(Glocal(q));
Glocalnum = reshape(Glocalnum, 4, 4, length(q));

for i = 1:length(q)
   
    transforms{i}.Matrix = Glocalnum(:,:,i);
    
end



end