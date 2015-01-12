#!/bin/bash
comp=$1
type=$2

echo "#Number Ele Time" > timing${comp}${type}
grep "total computation " ${1}_* >> timing${comp}${type}
sed -i "s/${comp}_${type}_//g" timing${comp}${type} 
sed -i 's/total computation/ /g' timing${comp}${type}
sed -i 's/_/ /g' timing${comp}${type}
sed -i 's/://g' timing${comp}${type}
sed -i 's/sec//g' timing${comp}${type}
sed -i 's/  */ /g' timing${comp}${type}
