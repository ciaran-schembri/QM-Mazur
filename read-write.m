
intrinsic Read(line::MonStgElt) -> Rec
  {turn the line of data into a record}

  split:=Split(line,"||");

  genus:=Integers()!eval(split[1]);
  fuchsindex:=eval(split[2]);
  torsioninvariants:=eval(split[4]);
  endogroup:=split[5];
  AutmuOnorms:=eval(split[6]);
  Hsplit:=eval(split[7]);

  RF := recformat< n : Integers(),
  genus,
  fuchsindex,
  torsioninvariants,
  endogroup,
  AutmuOnorms,
  Hsplit
  >;

  s := rec< RF | >;
  s`genus:=genus;
  s`fuchsindex:=fuchsindex;
  s`torsioninvariants:=torsioninvariants;
  s`endogroup:=endogroup;
  s`AutmuOnorms:=AutmuOnorms;
  s`Hsplit:=Hsplit;
  
  return s;
end intrinsic;


intrinsic GeneraTableToRecords(D::RngIntElt,del::RngIntElt,N::RngIntElt) -> Any 
  {}
  filename:=Sprintf("QM-Mazur/genera-tables/genera-D%o-deg%o-N%o.m",D,del,N);
  r:=Open(filename,"r");

  records:=[];
  while true do
    line :=Gets(r);
    if IsEof(line) then
      break;
    end if;

    if "<" in line then 
      s:=Read(line);
      Append(~records,s);
    end if;
  end while;

  return records;
end intrinsic;


intrinsic AbelianInvariantsToLatex(T::SeqEnum) -> MonStgElt 
 {}
  list:=[  [2],[2,2],[3],[2,3],[4], [2,4], [2,2,2], [3,3], [2,2,3],[3,4],[4,4], [2,2,4], [2,3,3] ];
  list_latex:=[  "(\\Z/2\\Z)", "(\\Z/2\\Z)^2", "(\\Z/3\\Z)", "(\\Z/2\\Z) \\times (\\Z/3\\Z)" ,  "(\\Z/4\\Z)", 
  	"(\\Z/2\\Z) \\times (\\Z/4\\Z)", "(\\Z/2\\Z)^3", "(\\Z/3\\Z)^2", "(\\Z/2\\Z)^2 \\times (\\Z/3\\Z)", 
  	"(\\Z/3\\Z) \\times (\\Z/4\\Z)", "(\\Z/4\\Z)^2", "(\\Z/2\\Z)^2 \\times (\\Z/4\\Z)", "(\\Z/2\\Z) \\times (\\Z/3\\Z)^2"];
  
  index:=Index(list,T);
  return list_latex[index];
 end intrinsic;


intrinsic GroupToLatex(G::MonStgElt) -> MonStgElt
  {}

  possible_endogroup:=  [ " C2 ", " C2^2 ", " D2 ", " D3 ", " S3 ", " D4 ", " D6 "];
  endogroup_latex:=     [  "C_2",  "D_2",   "D_2", "D_3",  "D_3", "D_4", "D_6" ];

  index:=Index(possible_endogroup,G);
  return endogroup_latex[index];
end intrinsic;



