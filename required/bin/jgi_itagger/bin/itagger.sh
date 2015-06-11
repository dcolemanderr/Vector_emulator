#!/bin/bash -l

module load itagger
exec itagger.pl $*
