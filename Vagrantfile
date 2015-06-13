# -*- mode: ruby -*-
# vi: set ft=ruby :

# All Vagrant configuration is done below. The "2" in Vagrant.configure
# configures the configuration version (we support older styles for
# backwards compatibility). Please don't change it unless you know what
# you're doing.

Vagrant.configure(2) do |config|

# Every Vagrant development environment requires a box. You can search for
# boxes at https://atlas.hashicorp.com/search.

config.vm.box = "scientific-linux-6"

# The url from where the 'config.vm.box' box will be fetched if it
# doesn't already exist on the user's system.

config.vm.box_url = "http://lyte.id.au/vagrant/sl6-64-lyte.box"

# Create a forwarded port mapping which allows access to a specific port
# within the machine from a port on the host machine. In the example below,
# accessing "localhost:8080" will access port 80 on the guest machine.
# config.vm.network "forwarded_port", guest: 80, host: 8080

# Create a private network, which allows host-only access to the machine
# using a specific IP.
# config.vm.network "private_network", ip: "192.168.33.10"

# Create a public network, which generally matched to bridged network.
# Bridged networks make the machine appear as another physical device on
# your network.
# config.vm.network "public_network"

# Share an additional folder to the guest VM. The first argument is
# the path on the host to the actual folder. The second argument is
# the path on the guest to mount the folder. And the optional third
# argument is a set of non-required options.

config.vm.synced_folder "/Users/colemanderr/repo/test-sci-linux-vm/required", "/home/vagrant/required"

# Provider-specific configuration so you can fine-tune various
# backing providers for Vagrant. These expose provider-specific options.
# Example for VirtualBox:
#
# config.vm.provider "virtualbox" do |vb|
#   # Display the VirtualBox GUI when booting the machine
#   vb.gui = true
#
#   # Customize the amount of memory on the VM:
#   vb.memory = "1024"
# end
#
# View the documentation for the provider you are using for more
# information on available options.

# Define a Vagrant Push strategy for pushing to Atlas. Other push strategies
# such as FTP and Heroku are also available. See the documentation at
# https://docs.vagrantup.com/v2/push/atlas.html for more information.
# config.push.define "atlas" do |push|
#   push.app = "YOUR_ATLAS_USERNAME/YOUR_APPLICATION_NAME"
# end

# Enable provisioning with a shell script. Additional provisioners such as
# Puppet, Chef, Ansible, Salt, and Docker are also available. Please see the
# documentation for more information about their specific syntax and use.

## Provisions ## 

#adding modified .bash profile to synced directory; loads upon vagrant ssh
#config.vm.provision "file", source: "./required/.bash_profile", destination: "/home/vagrant/.bash_profile"

#getting required files for itagger that cannot be downloaded with wget or installed with yum 
#config.vm.provision "file", source: "./required/get-pip.py", destination: "get-pip.py"

#Usearch is a C binary file, no make oir install required. Moved executable to /bin below in shell provision
#config.vm.provision "file", source: "./required/usearch7.0.1090_i86linux32", destination: "usearch"

#tarred files that are not always downloadable: connection refused.
#config.vm.provision "file", source: "./required/libgtextutils-0.6.tar.bz2", destination: "libgtextutils-0.6.tar.bz2"
#config.vm.provision "file", source: "./required/fastx_toolkit-0.0.12.tar.bz2", destination: "fastx_toolkit-0.0.12.tar.bz2"
#config.vm.provision "file", source: "./required/FLASH-1.2.11.tar.gz", destination: "FLASH-1.2.11.tar.gz"

#Shell provision for majority of itagger dependencies
config.vm.provision "shell", inline: <<-SHELL

#replaces default .bash_profile loaded with vagrant ssh
sudo cp /home/vagrant/required/.bash_profile /home/vagrant/

#makes all dependencies executable
chmod -R 755 /home/vagrant/required/*

## Editor ## 
# adds nano as editor
sudo yum install -y nano.x86_64

## PERL ##
# Adds latest version of perlbrew and perl v5.18.4, takes a while.
#export PERLBREW_ROOT=/home/vagrant/required/perl5/perlbrew
curl -L http://install.perlbrew.pl | bash
#sudo yum install -y perl-cpan
#sudo cpan App::perlbrew
#perlbrew init
#sudo chown -R vagrant: /home/vagrant/required/perl5/perlbrew
#source /home/vagrant/required/perl5/perlbrew/etc/bashrc
source ~/perl5/perlbrew/etc/bashrc
#sudo chown -R root: /home/vagrant/required/perl5/perlbrew
#chmod 710 -R /home/vagrant/required/perl5/perlbrew
#source /home/vagrant/required/bin/perlbrew/etc/bashrc
perlbrew install perl-5.18.4
perlbrew switch perl-5.18.4

## Cutadapt ##
# adds devel version of python necessary for gcc compile of cutadapt
sudo yum install -y python-devel.x86_64
#installing pip, and then cutadapt
python /home/vagrant/required/bin/get-pip.py    
pip install cutadapt

## Fastx_toolkit ##
# installing fastx_toolkit, and (first) it's dependencies : instructions from http://hannonlab.cshl.edu/fastx_toolkit/install_centos.txt
# wget http://cancan.cshl.edu/labmembers/gordon/files/libgtextutils-0.6.tar.bz2
# wget http://cancan.cshl.edu/labmembers/gordon/files/fastx_toolkit-0.0.12.tar.bz2 
sudo yum install pkgconfig.x86_64 gcc.x86_64 gcc-c++.x86_64 wget.x86_64
#tar -xjf libgtextutils-0.6.tar.bz2
cd /home/vagrant/required/bin/libgtextutils-0.6
./configure
make
sudo make install
cd ~
#rm libgtextutils-0.6.tar.bz2
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH
#tar -xjf fastx_toolkit-0.0.12.tar.bz2 ßß
cd /home/vagrant/required/bin/fastx_toolkit-0.0.12
./configure
make
sudo make install
cd ~
#rm fastx_toolkit-0.0.12.tar.bz2

## FLASH ##
# getting devel version of zlib for installing FLASH
sudo yum install -y zlib-devel.x86_64
# #getting FLASh 1.2.11
# downloaded from http://sourceforge.net/projects/flashpage/files/FLASH-1.2.11.tar.gz
#tar xfz FLASH-1.2.11.tar.gz
cd /home/vagrant/required/bin/FLASH-1.2.11
make
sudo cp flash /usr/local/bin
cd ~
#rm FLASH-1.2.11.tar.gz

## Java ##
sudo yum install -y java-1.7.0-openjdk.x86_64



SHELL

end
