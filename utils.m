


intrinsic getLines(file::MonStgElt) -> Any
  {turn text file into lines}
  F := Open(file, "r");
  Ls := [];
  while true do
    s := Gets(F);
    if IsEof(s) then break; end if;
    Append(~Ls,s);
  end while;
  return Ls;
end intrinsic;

intrinsic read_data(file::MonStgElt) -> Any
  {read the data file and make into a list with each line a record}

  RF := recformat< n : Integers(),
  IgusaClebsch,
  TorsionHeursiticTwist
   >;
  s := rec< RF | >;
  list:=[];
  lines:=getLines(file);
  for line in lines do
    items:=Split(line,"|");
    s`IgusaClebsch:=eval(items[1]);
    s`TorsionHeursiticTwist:=Integers()!eval(items[2]);
    Append(~list,s);
  end for;

  return list;
end intrinsic;

intrinsic HeuristicTorsionOrder4(file::MonStgElt) -> Any
  {return the curves who have an order 4 torsion subgroups from the heuristics.}
  lines:=read_data(file);
  data:=[];
  for s in lines do
    if IsDivisibleBy(s`TorsionHeursiticTwist,4) then
      Append(~data,s);
    end if;
  end for;
  return data;
end intrinsic;
