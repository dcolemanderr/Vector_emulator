#begin .profile

#Default umask: user rw, group r, other none
umask 022

alias lsl='ls -l -a -G'

export CLICOLOR=1
export LSCOLORS=gxBxhxDxfxhxhxhxhxcxcx

export PS1="\u@\w > "

PATH=$PATH:$HOME/bin
PATH=$PATH:$HOME/lib
PATH=$PATH:$HOME/scripts
#PATH=$PATH:$HOME/ActivePerl-5.18/lib/
#PATH=$HOME/jgi_itagger/bin/:$PATH
#PATH=$HOME/ActivePerl-5.18/bin/:$PATH

#export PERL5LIB=/global/home/users/colemand/ActivePerl-5.18/lib/

export PATH
export TMPDIR=/global/scratch/
export RDP_JAR_PATH=/required/rdp_classifier_2.5/rdp_classifier-2.5.jar
#export EDITOR=/usr/bin/nano

PATH=$PATH:$HOME/bin
export PATH

