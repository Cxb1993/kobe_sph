#!/bin/bash - 
#===============================================================================
#
#          FILE: go.sh
# 
#         USAGE: ./go.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: YOUR NAME (), 
#  ORGANIZATION: 
#       CREATED: 10/30/2013 07:22:41 PM JST
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

rm data?????
rm sph?????.png

./sph

#for num in $(seq 0 300)
for num in $(seq 0 600)
do
  fnum=$(printf %05d $num)
  file=data$fnum
  echo $file
  #povray +W800 +H600 -K$fnum +Osph$fnum.png -Dfalse sph.pov
  povray -K$fnum +Osph$fnum.png -Dfalse sph.pov
done
