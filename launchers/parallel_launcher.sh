#!/bin/bash -l

#CMD_LAUNCH="oarsub --notify "mail:shaman.narayanasamy@uni.lu" -l core=12/nodes=1,walltime=2 -p "memnode='48'" -t idempotent -t besteffort"
#CMD_LAUNCH="oarsub --notify "mail:shaman.narayanasamy@uni.lu" -l core=2/nodes=1,walltime=2 -p "memnode='48'" -t idempotent -t besteffort"
#CMD_LAUNCH="oarsub --notify "mail:shaman.narayanasamy@uni.lu" -l core=12/nodes=1,walltime=24 -t idempotent -t besteffort"

CMD_LAUNCH="oarsub --notify "mail:shaman.narayanasamy@uni.lu" -l core=1/nodes=1,walltime=24 -t idempotent -t besteffort"

echo $CMD_LAUNCH

ROOTOUT=$GENOMES

ls $GENOMES | \grep Isolate_ | \grep -v "Predicted_Isolate_protein_DB" | while read ISOLATE
do
  SAMPLE=`echo $ISOLATE | cut -d$' ' -f1`
  INDIR="${ROOTOUT}/${SAMPLE}"
  OUTDIR="${ROOTOUT}/${SAMPLE}/NewAssembly"
  
  CMD="launchers/execution.sh $SAMPLE $INDIR"

  #echo "oarsub --notify \"mail:shaman.narayanasamy@uni.lu\" -n \"${SAMPLE}_newAssembly\" \"$CMD\""
  #oarsub --notify \"mail:shaman.narayanasamy@uni.lu\" -n \"${SAMPLE}_newAssembly\" -S \"$CMD\"

  $CMD_LAUNCH -n "${SAMPLE}_newAssembly" "$CMD"

done
