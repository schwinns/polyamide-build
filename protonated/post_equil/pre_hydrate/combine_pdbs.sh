#!/bin/bash

grep -v CONECT hydrated_close.pdb > tmp.pdb
grep -v END tmp.pdb > hydrated_initial.pdb
grep ATOM reservoirs.pdb >> hydrated_initial.pdb
echo END >> hydrated_initial.pdb

rm tmp.pdb