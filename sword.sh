#!/bin/sh
#

word=$1

for a in src/*;do echo $a;cat $a|grep -i $1; done | less
for a in include/*;do echo $a;cat $a|grep -i $1; done | less
for a in *;do echo $a;cat $a|grep -i $1; done | less


