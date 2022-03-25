#!/bin/bash

echo '$ head -2 PA_converted.gro > PA_hydrated.gro'
head -2 PA_converted.gro > PA_hydrated.gro
echo '$ grep PA PA_converted.gro >> PA_hydrated.gro'
grep PA PA_converted.gro >> PA_hydrated.gro
echo '$ grep SOL reservoirs.gro >> PA_hydrated.gro'
grep SOL reservoirs.gro >> PA_hydrated.gro
echo '$ tail -1 PA_converted.gro >> PA_hydrated.gro'
tail -1 PA_converted.gro >> PA_hydrated.gro
