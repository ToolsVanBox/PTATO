#!/usr/bin/env bash

levenshtein ()
{
    local -r -- target=$1
    local -r -- given=$2
    local -r -- targetLength=${#target}
    local -r -- givenLength=${#given}
    local -- alt
    local -- cost
    local -- ins
    local -- gIndex=0
    local -- lowest
    local -- nextGIndex
    local -- nextTIndex
    local -- tIndex
    local -A -- leven

    while (( $gIndex <= $givenLength )); do

        tIndex=0
        while (( $tIndex <= $targetLength )); do
            (( $gIndex == 0 )) && leven[0,$tIndex]=$tIndex
            (( $tIndex == 0 )) && leven[$gIndex,0]=$gIndex

            (( tIndex++ ))
        done

        (( gIndex++ ))
    done

    gIndex=0
    while (( $gIndex < $givenLength )); do

        tIndex=0
        while (( $tIndex < $targetLength )); do
            [[ "${target:tIndex:1}" == "${given:gIndex:1}" ]] && cost=0 || cost=1

            (( nextTIndex = $tIndex + 1 ))
            (( nextGIndex = $gIndex + 1 ))

            (( del = leven[$gIndex,$nextTIndex] + 1 ))
            (( ins = leven[$nextGIndex,$tIndex] + 1 ))
            (( alt = leven[$gIndex,$tIndex] + $cost ))

            (( lowest = $ins <= $del ? $ins : $del ))
            (( lowest = $alt <= $lowest ? $alt : $lowest ))

            leven[$nextGIndex,$nextTIndex]=$lowest

            (( tIndex++ ))
        done

        (( gIndex++ ))
    done

    echo $lowest
}
