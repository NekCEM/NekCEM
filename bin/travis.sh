#!/usr/env bash

WHITELIST="cem_drive.F cem_maxwell.F"


if [[ $TESTS == 1 ]]; then
    bin/runtests --np 2
elif [[ $STYLE == 1 ]]; then
    rm -f style.log
    for f in $WHITELIST; do
	bin/stylecheck src/$f >> style.log
    done
    find tests -name "*.usr" -exec bin/stylecheck {} \; >> style.log
    if [ -s style.log ]; then
	cat style.log
	rm style.log
	exit 1
    else
	rm style.log
	exit 0
    fi
fi
