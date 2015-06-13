#begin .profile

#Default umask: user rw, group r, other none
umask 022

alias lsl='ls -l -a -G'

export CLICOLOR=1
export LSCOLORS=gxBxhxDxfxhxhxhxhxcxcx

export PS1="\u@\w > "

PATH=$PATH:$HOME/bin
PATH=$PATH:$HOME/lib

#export PERL5LIB=/home/vagrant/required/bin/perlbrew/perls/perl-5.18.4/lib
sudo su root 
source /root/perl5/perlbrew/etc/bashrc
exit

export TMPDIR=/global/scratch/
export RDP_JAR_PATH=/required/rdp_classifier_2.5/rdp_classifier-2.5.jar
export PERLBREW_ROOT=/home/vagrant/required/perl5/perlbrew
#export EDITOR=/usr/bin/nano

export PATH

