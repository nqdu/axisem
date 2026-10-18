#!/bin/csh -f

if ( ${#argv} < 1 || "$1" == "-h" ) then
    echo "USAGE:"
    echo "  Argument 1:  name for the mesh"
    echo ""
    exit
endif

set homepath = `echo $PWD`
set meshname = "$1"

if ( ${#argv} == 2 ) then
  set xmeshrundir = $2
else
  set xmeshrundir = "."
endif

set bgmodel = `grep "^BACKGROUND_MODEL" "$xmeshrundir/inparam_mesh" | awk '{print $2}'`
if ( $bgmodel == 'external') then
  set extmodel = `grep "^EXT_MODEL" "$xmeshrundir/inparam_mesh" | awk '{print $2}'`
  # Use the external file's stem while retaining a requested period suffix.
  set extmodel = `echo "$extmodel" | sed 's/\"//g'`
  if ( "$meshname" == 'external' || "$meshname" =~ external.* ) then
    set model_name = "$extmodel:t:r"
    set name_suffix = "$meshname:s/external//"
    set meshname = "${model_name}${name_suffix}"
  endif
endif
set meshpath = "../SOLVER/MESHES/$meshname"
echo "Moving mesh from $xmeshrundir to $meshpath"

# test if mesher finished without problems:
if (`grep 'DONE WITH MESHER' $xmeshrundir/OUTPUT | wc -l` != '1') then
  echo "ERROR: MESHER did not finish yet or encountered a problem."
  echo "       Check 'OUTPUT' file for more details."
  exit
else
  echo 'MESHER finished smoothly...'
endif

if ( ! -d ../SOLVER/MESHES ) then
  mkdir ../SOLVER/MESHES
endif

if ( ! -d $meshpath ) then
  mkdir $meshpath
else
  echo "ERROR: the mesh folder " $meshpath " already exists!"
  exit
endif

# convert relative to absolute path:
cd $meshpath
set meshpath = `echo $PWD`
cd $homepath

echo "Moving mesh to" $meshpath

mv $xmeshrundir/meshdb.dat* $meshpath
mv $xmeshrundir/mesh_params.h $meshpath
mv $xmeshrundir/OUTPUT $meshpath
#mv Diags $meshpath
cp -p $xmeshrundir/inparam_mesh $meshpath

echo $bgmodel
if ( $bgmodel == 'external') then
  set model_path = "$extmodel"
  if ( "$extmodel" !~ /* ) set model_path = "$xmeshrundir/$extmodel"
  echo "Copying external model file $model_path"
  if ( -f "$model_path") then
    cp -p "$model_path" $meshpath/external_model.bm
  else
    echo "external model file $model_path does not exist. Did you delete it? Why?"
    exit
  endif
endif

mv $xmeshrundir/Diags/* $meshpath

echo "Contents in" $meshpath ":"
ls $meshpath
cd $meshpath
echo "Total size: `du -sh` "
echo "DONE."
