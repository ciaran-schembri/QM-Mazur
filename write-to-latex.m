




MakeGeneraTable:= function(D,del,N : genus_range:=[0,1], IsSplit:=true)
  filename:=Sprintf("QM-Mazur/genera-tables/genera-D%o-deg%o-N%o.m",D,del,N);
  r:=Open(filename,"r");
  possible_tors:=[   [3], [2,3], [3,3], [4], [2,4], [2,2,2], [2,2,3],[3,4],[4,4], [2,2,4],[2,3,3] ];
  possible_endogroup:=  [ " C2 ", " C2^2 ", " D2 ", " D3 ", " S3 ", " D4 ", " D6 "];
  endogroup_latex:=     [  "C_2",  "D_2",   "D_2", "D_3",  "D_3", "D_4", "D_6" ];

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


  tors:= [ eval(Split(line,"||")[4]) : line in lines ];
  ParallelSort(~tors,~lines);

  for line in lines do
    split:=Split(line,"||");
    tors:=eval(split[4]);
    if tors ne [] then 
      torslatex1:= Sprintf("(\\Z/%o\\Z)", tors[1]);
      if #tors ge 2 then 
        torslatex2:=(&cat[ Sprintf(" \\times (\\Z/%o\\Z) ", tors[i]) : i in [2..#tors] ]);
      else 
        torslatex2:="";
      end if;
      torslatex:=torslatex1 cat torslatex2;
    else 
      torslatex:="1";
    end if;
    endogroup:=split[5];
  
    genus:=Integers()!eval(split[1]);
    Hsplit:=eval(split[7]);
    if tors in possible_tors and endogroup in possible_endogroup then
      //exists(e){ x : x in possible_endodeg | IsIsomorphic(endodeg,x) };  
      ind:=Index(possible_endogroup,endogroup);
      latex_group:=endogroup_latex[ind];
      if genus in genus_range and Hsplit eq IsSplit then 
        printf " $ %o $ & $%o$ & $%o$ & $ %o $ & $%o$ & %o \\\\\n", genus, torslatex, latex_group, split[6], split[2], split[7];  
        //line;
      end if;
    end if;
  end for;
  return "";
end function;


D:=6;
del:=1;
N:=4;
MakeGeneraTable(D,del,N);



for D in [6] do 
  for N in [4,6] do 
    for del in [1,2] do 
      "\\\\";
      printf "\\multicolumn{6}{c}{Discriminant of $\\calO$ is $%o$, $N=%o$ and degree of polarization is $%o$} \\\\ \n", D, N, del;
      "\\\\";
      MakeGeneraTable(D,del,N : genus_range:=[0], IsSplit:=true);
      MakeGeneraTable(D,del,N : genus_range:=[0], IsSplit:=false);
      MakeGeneraTable(D,del,N : genus_range:=[1], IsSplit:=true);
      MakeGeneraTable(D,del,N : genus_range:=[1], IsSplit:=false);
      MakeGeneraTable(D,del,N : genus_range:=[2], IsSplit:=true);
      MakeGeneraTable(D,del,N : genus_range:=[2], IsSplit:=false);
    end for;
  end for;
end for;
