#!/usr/env bash

SRC_WHITELIST="cem_drive.F cem_maxwell.F"
TST_WHITELIST="2dboxpec/2dboxpec.usr 2dboxper/2dboxper.usr
2dboxpml/2dboxpml.usr 2ddielectric/2ddielectric.usr
2dgraphene/2dgraphene.usr 3dboxpec/3dboxpec.usr 3dboxper/3dboxper.usr
3dboxpml/3dboxpml.usr 3ddielectric/3ddielectric.usr
3dgraphene/3dgraphene.usr acoustic2d/pec.usr cylwave/cylwave.usr
drude/drude.usr lorentz/lorentz.usr poisson2d/2dboxpec.usr
poisson2d/2dboxper.usr poisson2d/circpec.usr poisson3d/3dboxpec.usr
poisson3d/3dboxper.usr poisson3d/cylinder.usr"


if [[ $TESTS == 1 ]]; then
    bin/runtests --np 2
elif [[ $STYLE == 1 ]]; then
    rm -f style.log
    for f in $SRC_WHITELIST; do
	bin/stylecheck src/$f >> style.log
    done
    for f in $TST_WHITELIST; do
	bin/stylecheck tests/$f >> style.log
    done
    if [ -s style.log ]; then
	cat style.log
	rm style.log
	exit 1
    else
	rm style.log
	exit 0
    fi
fi
