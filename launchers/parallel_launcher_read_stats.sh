#!/bin/bash -l

CMD_LAUNCH="oarsub --notify "mail:shaman.narayanasamy@uni.lu" -l core=1/nodes=1,walltime=1 -t idempotent -t besteffort"

#echo $CMD_LAUNCH

ROOTOUT=$GENOMES

ls $GENOMES | \grep Isolate_ | \grep -v "Predicted_Isolate_protein_DB" | while read ISOLATE
do
  SAMPLE=`echo $ISOLATE | cut -d$' ' -f1`
  INDIR="${ROOTOUT}/${SAMPLE}"
  OUTDIR="${ROOTOUT}/${SAMPLE}/NewAssembly"
  
  CMD="launchers/execution_read_stats.sh $SAMPLE $INDIR"

  # This is for testing purposes
  $CMD_LAUNCH -n "${SAMPLE}_readstats" "$CMD"

done
