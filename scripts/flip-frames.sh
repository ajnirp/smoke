#!/bin/bash

for f in frames/*.ppm
do
  convert -flip $f $f
done