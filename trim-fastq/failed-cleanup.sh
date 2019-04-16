#!/bin/bash

grep FAILED trace.* | awk '{print $2}' | xargs -I % rm -rf work/%*