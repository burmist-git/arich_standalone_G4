#!/bin/bash

########################################################################
#                                                                      #
# Copyright(C) 2017 - LBS - (Single person developer.)                 #
# Wed Mar 28 14:23:47 CEST 2018                                        #
# Autor: Leonid Burmistrov                                             #
#                                                                      #
# File description:                                                    #
#                  This script setup standalone G4 ARICH simulation.   #
#                                                                      #
# Input paramete:                                                      #
#                   NON.                                               #
#                                                                      #
# This software is provided "as is" without any warranty.              #
#                                                                      #
########################################################################

setup_standalone_ARICHG4_bash () {

    useroot60806;
    #usegeant4_10_01_p02;
    usegeant41004

    printenv | grep G4
    printenv | grep ROOT

    #xmlarichdataHome="$PWD/xmlarichdata/"
    #xmlarichdataHome="/home/burmist/CADMesh_v1.1_standAlone/xmlarichdata/"
    xmlarichdataHome="/home/burmist/xmlarichdata/"
    export XMLARICHDATAHOME=$xmlarichdataHome
    #export XMLARICHDATA_BASF2_INCLUDE="/home/b-lab050/KEK/releases/head_16.05.2018/include/"
    export XMLARICHDATA_BASF2_INCLUDE="/home/burmist/home2/KEK/releases/head_01.10.2018/include/"
    export XMLARICHDATA_Fully_Compiled="yes"
    arich_buildHome="$PWD/arich-build/"
    aerogel_transmittance_length_rayleigh_scatter="$xmlarichdataHome/rootdata/aerogel_xml_ver3/"
    rayleigh_scatter_l="aerogel_transmittance_vs_photon_wavelength_aerogel_xml_ver3.root"
    ln -s -f $xmlarichdataHome/indata/qe_vs_lam.dat $arich_buildHome/qe_vs_lam.dat
    ln -s -f $xmlarichdataHome/indata/n_A_Lam.dat $arich_buildHome/n_A_Lam.dat
    ln -s -f $xmlarichdataHome/indata/n_B_Lam.dat $arich_buildHome/n_B_Lam.dat
    ln -s -f $aerogel_transmittance_length_rayleigh_scatter/$rayleigh_scatter_l $arich_buildHome/$rayleigh_scatter_l
    #ln -s -f ../xmlarichdata/qe_vs_lam.dat ./anaarich/qe_vs_lam.dat

    export CADMESH_INSTALL="/home/burmist/CADMesh_v1.1/CADMesh-install/"
    export CADMESH_INCLUDE_DIRS="$CADMESH_INSTALL/include/"
    export CADMESH_LIBRARIES="$CADMESH_INSTALL/lib/"
    
    echo " "
    echo "XMLARICHDATAHOME           = $XMLARICHDATAHOME"
    echo "XMLARICHDATA_BASF2_INCLUDE = $XMLARICHDATA_BASF2_INCLUDE"    
    echo " "
    echo " "
    echo "CADMESH_INSTALL            = $CADMESH_INSTALL"
    echo "CADMESH_INCLUDE_DIRS       = $CADMESH_INCLUDE_DIRS"    
    echo "CADMESH_LIBRARIES          = $CADMESH_LIBRARIES"    
    echo " "
    

}

#diffL=`printenv | grep BELLE2 | wc -l`
#echo "diffL = $diffL"
#if [ $diffL -eq 0 ]
#then
#    setup_standalone_ARICHG4_bash
#else
#    echo "ERROR ---> please open new terminal without BELLE2 framework"
#fi

setup_standalone_ARICHG4_bash
