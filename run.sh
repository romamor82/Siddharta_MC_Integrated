#!/bin/bash

# 1 -> number of events
# 2 -> config file (CARD.dat)
# 3 -> out file

nohup ./sidd_mc $1 $2 > $3 &

