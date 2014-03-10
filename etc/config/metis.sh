#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM.
#
#     OpenFOAM is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
#
# File
#     config/metis.sh
#
# Description
#     Setup file for metis include/libraries.
#     Sourced during wmake process only.
#
# Note
#     A csh version is not needed, since the values here are only sourced
#     during the wmake process
#------------------------------------------------------------------------------

export METIS_VERSION=metis-5.1.0
export METIS_ARCH_PATH=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$METIS_VERSION

if [ -n "$WM_USE_MACPORT" ]
then
    if [ -e "/opt/local/include/metis.h" ]
    then
	export METIS_ARCH_PATH=/opt/local
	unset METIS_VERSION
    else
	echo "No metis in MacPorts. Install metis with 'port install metis'"
    fi
fi

# -----------------------------------------------------------------------------
