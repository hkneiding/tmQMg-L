%chk=<name>
%mem=5GB

#p pbe1pbe/def2tzvp empiricaldispersion=gd3bj pop=(allorbitals,threshorbitals=1,npa,full) integral=noxctest

<name>

 <charge> 1
<xyz>

