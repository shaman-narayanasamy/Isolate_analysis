#!/bin/bash -l

ROOTOUT=$GENOMES

ls $GENOMES | grep Isolate_ | grep -v "Predicted_Isolate_protein_DB" | while read ISOLATE
do
  SAMPLE=`echo $ISOLATE | cut -d$' ' -f1`
  INDIR="${ROOTOUT}/${SAMPLE}"
  OUTDIR="${ROOTOUT}/${SAMPLE}/NewAssembly"
  
  CMD="launchers/execution_rast.sh $SAMPLE $INDIR"

  #echo "$CMD_LAUNCH -n \"${SAMPLE}_newAssembly\" \"$CMD\""
  echo $CMD
  $CMD

done
