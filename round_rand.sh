#!/bin/bash

echo "#################### Begin round_rand.sh"
function array_shift {
    ARR=$1
    l=${#ARR[@]}
    for (( i=1; i<$l; i++ ))
    do
        j=$(($i - 1))
        ARR[$j]=${ARR[$i]}
    done
    unset ARR[$(($l - 1))]
}

function normalize_name {
    $(echo $1 | sed -r "s/(GCA_[0-9]*.[0-9]_)(.*)(protein.faa)/\1\3/g")
}

BASE=$1
PATHFASTA=$2
shift 2
echo "Base:"
echo ${BASE}
echo "Fasta:"
echo ${PATHFASTA}
LENGTH=$#
if [ ! -d "${BASE}/${LENGTH}" ]; then
    mkdir ${BASE}/${LENGTH}
fi
            
echo LENGTH = ${LENGTH}
echo "$@"

PATHHOMOLOGOUS="${BASE}/${LENGTH}/output_find_homologous"
PATHBLAST="${BASE}/${LENGTH}/output_blastp"
if [ ! -d "${PATHHOMOLOGOUS}" ]; then
    mkdir ${PATHHOMOLOGOUS}
    mkdir ${PATHHOMOLOGOUS}/round_1
    mkdir ${PATHHOMOLOGOUS}/round_1/output
fi
if [ ! -d "${PATHBLAST}" ]; then
    mkdir ${PATHBLAST}
    mkdir ${PATHBLAST}/REF_1
    mkdir ${PATHBLAST}/PARSER_1
fi
CURRENT_REF="$1"
./blastp.pl ${PATHBLAST}/REF_1 ${PATHBLAST}/PARSER_1 "$@"
./find_homologous.pl 1 ${PATHFASTA} ${PATHBLAST} ${PATHHOMOLOGOUS} ${PATHBLAST}/REF_1/*.txt

NORMALIZED=$(echo ${CURRENT_REF} | sed -r "s/(GCA_[0-9]*.[0-9]_)(.*)(protein.faa)/\1\3/g")
TOMOVE="${PATHHOMOLOGOUS}/round_1/output/$(basename ${NORMALIZED})"
echo "Current ref: ${CURRENT_REF}"
echo "To move ${TOMOVE}"
if [ -e "${TOMOVE}" ]; then
   echo "Moving ${TOMOVE} to .."
   mv ${TOMOVE} ${PATHHOMOLOGOUS}/round_1
fi

for ((i=2; i <=${LENGTH}; i++))
do  
  PREVROUND=$(($i - 1))
  PATHROUND="${PATHHOMOLOGOUS}/round_${PREVROUND}"
  NEWPATHROUND="${PATHHOMOLOGOUS}/round_$i"
  PATHREF="${PATHBLAST}/REF_$i"
  PATHPARSER="${PATHBLAST}/PARSER_$i"
  if [ ! -d "${NEWPATHROUND}" ]; then
      mkdir ${NEWPATHROUND}/
      mkdir ${NEWPATHROUND}/output
  fi
  if [ ! -d "${PATHREF}" ]; then
    mkdir ${PATHREF}
  fi
  if [ ! -d "${PATHPARSER}" ]; then
    mkdir ${PATHPARSER}
  fi
  CURRENT_REF_L=($(ls ${PATHROUND}/output/*.faa))
  CURRENT_REF=$(basename ${CURRENT_REF_L[0]})
  ./blastp.pl ${PATHREF} ${PATHPARSER} ${PATHROUND}/output/*.faa
  ./find_homologous.pl $i ${PATHFASTA} ${PATHBLAST} ${PATHHOMOLOGOUS} ${PATHREF}/*.txt
  TOMOVE="${NEWPATHROUND}/output/${CURRENT_REF}"
  echo "Move ${TOMOVE} to .."
  if [ -e "${TOMOVE}" ]; then
     echo "Moving ${TOMOVE} to .."
     mv ${TOMOVE} ${NEWPATHROUND}
  fi
done
