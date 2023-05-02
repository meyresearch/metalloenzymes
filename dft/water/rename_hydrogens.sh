#!/bin/bash

for i in $(seq 1 5)
do
	sed -i 's/HO/HW/g' WT$i.mol2
done
