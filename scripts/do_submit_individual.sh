#!/bin/bash

Nx=$1
nx=$2
va=$3
tau=$4
lambda=$5
i=$6

input="seeds/seeds_Nx=${Nx}_nx=${nx}_va=${va}_tau=${tau}_lambda=${lambda}.txt"


sbatch -J "active_net_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_nx=${nx}_seed=${i}" -o "log/active_net_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_nx=${nx}_seed=${i}.o%j" -e "log/active_net_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_nx=${nx}_seed=${i}.e%j" scripts/submit_sims.sh "/home/laynefrechette/active-network/conf/active_net_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_nx=${nx}.conf" $i $input
