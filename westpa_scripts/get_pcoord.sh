#!/bin/bash
set -x

cd $WEST_STRUCT_DATA_REF
python $WEST_SIM_ROOT/common_files/init_pcoord.py
# paste $WEST_STRUCT_DATA_REF/pcoord1.dat $WEST_STRUCT_DATA_REF/pcoord2.dat $WEST_STRUCT_DATA_REF/pcoord3.dat > $WEST_PCOORD_RETURN 
paste $WEST_STRUCT_DATA_REF/pcoord2.dat $WEST_STRUCT_DATA_REF/pcoord1.dat $WEST_STRUCT_DATA_REF/pcoord3.dat > $WEST_PCOORD_RETURN
