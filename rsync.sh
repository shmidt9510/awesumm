#!/bin/bash

dir=awesumm/WorkCases/olegrog

rsync -auv --delete -e ssh hilbert:$dir/* $HOME/$dir/
