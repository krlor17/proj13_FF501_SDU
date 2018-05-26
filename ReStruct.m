function matrix = ReStruct(vec,tpoints,tf)
%this function takes a column vector $vec and returns a matrix $matrix
%the matrix are of the dimensions tpoints x #runs
%if $tf=0 (true) then the values of $matrix are the log of the $vec values
%if $tf=1 (false) then the values of $matrix are the same as the $vec values
%if $tf is neither 1 or 0, $matrix contains only zeroes (error) 
	%% we consider only half the points, as we are not
	%% interested in the parity partner
	t=round(tpoints/2);
	len=round(length(vec)/tpoints);
	matrix=zeros(t, len);
	if tf==0
		for i = 1:t
    			for k= 1: (tpoints*len)
        			if mod(k,tpoints) == mod(i,tpoints)         % If t_k = t_i,
	 				% then put respective data in matrix entry [i, k mod 128 +1]
					matrix(i,ceil(k/tpoints))=log(vec(k));
				end
    			end
		end
	elseif tf==1
	for i = 1:t
    			for k= 1: (tpoints*len)
        			if mod(k,tpoints) == mod(i,tpoints)         % If t_k = t_i,
	 				% then put respective data in matrix entry [i, k mod 128 +1]
					matrix(i,ceil(k/tpoints))=vec(k);
				end
    			end
		end
	end
end
