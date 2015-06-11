#!/bin/sh

DATE=`date +%Y-%m-%d-%H:%M:%S`
COMMANDS_HSI=./hpss_backup_commands_hpss_hsi_itags_$DATE.sh
COMMANDS_HTAR=./hpss_backup_commands_hpss_htar_itags_$DATE.sh

jgi_itaggerRunHPSSBackup.pl $COMMANDS_HSI $COMMANDS_HTAR $@

## Then source the file containing hsi commands.
. $COMMANDS_HSI
. $COMMANDS_HTAR

