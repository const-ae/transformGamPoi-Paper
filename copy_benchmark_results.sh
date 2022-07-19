#!/bin/bash

# Remove old folder content
rm -r benchmark
# Make new folder
mkdir benchmark
mkdir benchmark/renv

# Copy relevant files (i.e., exclude renv library because it is too big)
scp -r ahlmanne@seneca:~/projects/transformGamPoi-Paper2/{README.md,job_overview.yaml,renv.lock,.Rprofile} benchmark/.
scp -r ahlmanne@seneca:~/projects/transformGamPoi-Paper2/src benchmark/src
scp -r ahlmanne@seneca:~/projects/transformGamPoi-Paper2/renv/activate.R benchmark/renv/activate.R
scp -r ahlmanne@seneca:~/projects/transformGamPoi-Paper2/output benchmark/output
