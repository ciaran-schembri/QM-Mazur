for D in [6,15] do
  for N in [4,6] do 
    for del in [1,2] do 
	  rhocirc:=AllEnhancedSubgroups(D,del,N : PQMtorsion:=true, verbose:=true, minimal:=false, lowgenus:=true);
	  print "";
	end for;
  end for;
end for;
