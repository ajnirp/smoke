#!/bin/bash

avconv -i frames/frame_%04d.ppm -r 30 -b 65536k smoke.mp4