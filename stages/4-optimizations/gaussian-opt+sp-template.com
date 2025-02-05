%chk=<name>
%mem=5GB

#p pbepbe/def2svp opt freq empiricaldispersion=gd3bj

OPT <name>

 <charge> 1
<xyz>

--Link1--
%Chk=<name>

#p pbe1pbe/def2tzvp geom=allcheckpoint guess=read empiricaldispersion=gd3bj pop=(allorbitals,threshorbitals=1,npa,full) volume

SP <name>
