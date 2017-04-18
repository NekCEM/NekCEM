#!/usr/env bash

WHITELIST="drive.F maxwell.F utilities.F"


if [[ $TESTS == 1 ]]; then
    pytest --np 2
    if [[ $? != 0 ]]; then
	cat tests/build.log
	exit 1
    fi
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
