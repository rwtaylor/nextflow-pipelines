#!/bin/bash
rsync -avh --progress --exclude stages outputs/* ryan@bio-dap15.stanford.edu:/zstor/2016-awd/alignments/align-2

