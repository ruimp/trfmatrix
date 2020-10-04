#!/bin/bash
for i in `seq 20 5 40`
do
	python interf.py $i &
done
