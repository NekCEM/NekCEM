#!/usr/env bash


if [[ $TESTS == 1 ]]; then
    pytest --build_command makenekmpi --arch linux-gnu-mpi --np 2
    if [[ $? != 0 ]]; then
	cat tests/build.log
	exit 1
    fi
elif [[ $STYLE == 1 ]]; then
    find tests -name "*.usr" -exec bin/stylecheck {} \; > style.log
    if [ -s style.log ]; then
	cat style.log
	exit 1
    else
	exit 0
    fi
fi
