




MakeGeneraTable:= function(D,del,N)
  filename:=Sprintf("QM-Mazur/genera-tables/genera-D%o-deg%o-N%o.m",D,del,N);
  r:=Open(filename,"r");
  possible_tors:=[   [3], [2,3], [3,3], [4], [2,4], [2,2,2], [2,2,3],[3,4],[4,4], [2,2,4],[2,3,3] ];
  possible_endodeg:=  [" C2 ", " C2^2 ", " D2 ", " D3 ", " S3 ", " D4 ", " D6 "];

  lines:=[];
  while true do
    line :=Gets(r);
    if IsEof(line) then
      break;
    end if;
    
    if "<" in line then 
      split:=Split(line,"<");
      if split[1] notin lines then 
        //split[1];
        Append(~lines,split[1]);
      end if;
    end if;
  end while;

  for line in lines do
    split:=Split(line,"||");
    tors:=eval(split[4]);
    endodeg:=split[5];
    if tors in possible_tors and endodeg in possible_endodeg then
      //exists(e){ x : x in possible_endodeg | IsIsomorphic(endodeg,x) };  
      line;
    end if;
  end for;
  return "";
end function;


D=6;
del:=1;
N:=4;
MakeGeneraTable(D,del,N);
