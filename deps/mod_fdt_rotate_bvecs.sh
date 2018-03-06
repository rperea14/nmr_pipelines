#!/bin/bash
#Created by Rodrigo Perea
#Goal:  Rotate the *.bvecs into fslstd space based on the informatin in the
#       *.matrix (or 1st argument passed). 
#       Caution: This should only be access through the
#       matlab diffusion pipeline or run it at your own risk.
#       Dependencies: DWI_HAB pipelines and FSL tools (for ivscale)
#       ! This script only checks sign changes and doesn't take into account shearing, scaling, rotating orientation (as is not needed for HAB)


#$1 --> should be the input *.bvec
bvecs_IN=$1
#$2 --> should be the *.matrix
matrix_IN=$2
#$3 --> should be the *.(output) bvec
bvecs_OUT=$3

xrot=$( avscale ${matrix_IN} | grep Rotation -A 1 | tail -n -1 | awk '{print $1}' )
yrot=$( avscale ${matrix_IN} | grep Rotation -A 2 | tail -n -1 | awk '{print $2}' )
zrot=$( avscale ${matrix_IN} | grep Rotation -A 3 | tail -n -1 | awk '{print $3}' )

#echo $xrot 
#echo $yrot 
#echo $zrot

if [ -f $bvecs_OUT ] ; then
        rm $bvecs_OUT 
fi
echo "In $bvecs_IN ..."
while read -r cur_XYZ ; do
        #echo "Values: $cur_XYZ" 
        #For every line in bvecs:
        cur_X=`echo ${cur_XYZ} | awk '{ print $1 }'  `
        cur_Y=`echo ${cur_XYZ} | awk '{ print $2 } '` 
        cur_Z=`echo ${cur_XYZ} | awk '{ print $3 } '` 
        #Multiply it with its xyz_rot values:
        cur_Xmod=`echo "${cur_X} * $xrot" | bc -l | awk   '{printf "%.5f\n", $1}'  `
        cur_Ymod=`echo "${cur_Y} * $yrot" | bc -l | awk '{printf "%.5f\n", $1}' `
        cur_Zmod=`echo "${cur_Z} * $zrot" | bc -l | awk '{printf "%.5f\n", $1}'`
        #echo "New Values: $cur_Xmod $cur_Ymod $cur_Zmod"
        echo "$cur_Xmod $cur_Ymod $cur_Zmod" >> $bvecs_OUT
done < $bvecs_IN



