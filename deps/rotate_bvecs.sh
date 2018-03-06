#!/bin/bash
#
#
# Usage:
#  ./rotate_bvecs <x> <y> <z> <affXFM file>
#

x=$1;
y=$2;
z=$3;
affXFMfile=$4;



m11=`sed -n "1p" ${affXFMfile} | awk '{print $1}'`;
m11=`echo "${m11}" | sed 's/e+/*10^/g;s/ /*/' | sed 's/e-/*10^-/g;s/ /*/' | bc -l`;
m12=`sed -n "1p" ${affXFMfile} | awk '{print $2}'`; 
m12=`echo "${m12}" | sed 's/e+/*10^/g;s/ /*/' | sed 's/e-/*10^-/g;s/ /*/' | bc -l`;
m13=`sed -n "1p" ${affXFMfile} | awk '{print $3}'`; 
m13=`echo "${m13}" | sed 's/e+/*10^/g;s/ /*/' | sed 's/e-/*10^-/g;s/ /*/' | bc -l`;
m21=`sed -n "2p" ${affXFMfile} | awk '{print $1}'`; 
m21=`echo "${m21}" | sed 's/e+/*10^/g;s/ /*/' | sed 's/e-/*10^-/g;s/ /*/' | bc -l`;
m22=`sed -n "2p" ${affXFMfile} | awk '{print $2}'`; 
m22=`echo "${m22}" | sed 's/e+/*10^/g;s/ /*/' | sed 's/e-/*10^-/g;s/ /*/' | bc -l`;
m23=`sed -n "2p" ${affXFMfile} | awk '{print $3}'`; 
m23=`echo "${m23}" | sed 's/e+/*10^/g;s/ /*/' | sed 's/e-/*10^-/g;s/ /*/' | bc -l`;
m31=`sed -n "3p" ${affXFMfile} | awk '{print $1}'`; 
m31=`echo "${m31}" | sed 's/e+/*10^/g;s/ /*/' | sed 's/e-/*10^-/g;s/ /*/' | bc -l`;
m32=`sed -n "3p" ${affXFMfile} | awk '{print $2}'`;
m32=`echo "${m32}" | sed 's/e+/*10^/g;s/ /*/' | sed 's/e-/*10^-/g;s/ /*/' | bc -l`;
m33=`sed -n "3p" ${affXFMfile} | awk '{print $3}'`; 
m33=`echo "${m33}" | sed 's/e+/*10^/g;s/ /*/' | sed 's/e-/*10^-/g;s/ /*/' | bc -l`;

rX=`echo "scale=6; ($m11 * $x) + ($m12 * $y) + ($m13 * $z);" | bc -l`;
rY=`echo "scale=6; ($m21 * $x) + ($m22 * $y) + ($m23 * $z);" | bc -l`;
rZ=`echo "scale=6; ($m31 * $x) + ($m32 * $y) + ($m33 * $z);" | bc -l`;
#echo $(echo "scale=6; ${rX};" | bc -l)  $(echo "scale=6; ${rY};" | bc -l) $(echo "scale=6; ${rZ};" | bc -l); 

round()
{
echo $(printf %.$2f $(echo "scale=$2;(((10^$2)*$1)+0.5)/(10^$2)" | bc))
};

echo $(round $rX 6)  $(round $rY 6) $(round $rZ 6)