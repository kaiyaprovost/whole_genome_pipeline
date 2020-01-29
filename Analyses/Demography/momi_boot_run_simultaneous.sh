#!/bin/bash

## runs a job file $1 100 times, each python file does 10 boots so init 10 jobs
for i in {1..10}; do
	qsub $1
done
