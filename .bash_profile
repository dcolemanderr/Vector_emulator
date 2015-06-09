#begin .profile

#Default umask: user rw, group r, other none
umask 022

alias lsl='ls -l -a -G'
alias genepool='ssh colemand@genepool.nersc.gov'
alias cgrl='ssh colemand@hpc.brc.berkeley.edu'

export CLICOLOR=1
export LSCOLORS=gxBxhxDxfxhxhxhxhxcxcx

#export EDITOR=/usr/bin/nano

PATH=$PATH:$HOME/bin
export PATH

