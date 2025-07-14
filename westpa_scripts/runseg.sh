#!/bin/bash


if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF

ln -sv $WEST_SIM_ROOT/common_files/abl1_i1_eq_drude.pdb .

if [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_CONTINUES" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/Production_drude.py > Production_drude.py
  cp $WEST_PARENT_DATA_REF/seg.rst ./seg.rst
elif [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_NEWTRAJ" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/Production_drude.py > Production_drude.py 
  cp $WEST_SIM_ROOT/common_files/basis.rst ./seg.rst
fi


# Run the dynamics with OpenMM
python Production_drude.py -psf $WEST_SIM_ROOT/common_files/abl1_i1_eq_drude.psf \
                           -crd $WEST_SIM_ROOT/common_files/abl1_i1_eq_drude.pdb \
                           -toppar $WEST_SIM_ROOT/common_files/toppar_drude.str \
                           -state seg.rst \
                           -temp 310 \
                           -runtime  100 \
                           -nstride  1   \
                           -savefreq 100 \
                           -dt 1 \
                           -outname seg


#Calculate pcoord with MDAnalysis
python3 $WEST_SIM_ROOT/common_files/get_dist.py

cat dih1.dat > $WEST_DIH1_RETURN
cat dih2.dat > $WEST_DIH2_RETURN

# paste angle1.dat dist.dat dist2.dat > $WEST_PCOORD_RETURN
paste dist.dat dist2.dat angle1.dat > $WEST_PCOORD_RETURN
# cp abl1_i1_eq_drude.pdb $WEST_TRAJECTORY_RETURN
# cp seg.dcd $WEST_TRAJECTORY_RETURN
	
# cp abl1_i1_eq_drude.pdb $WEST_RESTART_RETURN
# cp seg.pdb $WEST_RESTART_RETURN/parent.pdb

# cp seg.log $WEST_LOG_RETURN

# Clean up
rm -f Production_drude.py seg.dcd dih1.dat dih2.dat dist.dat dist2.dat angle1.dat
