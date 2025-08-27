picard -Xmx32g RevertSam \
  I=bams/E_coli_H10407_2.bam O=E_coli_H10407_2.unaligned.bam \
  REMOVE_ALIGNMENT_INFORMATION=true \
  RESTORE_ORIGINAL_QUALITIES=true \
  REMOVE_DUPLICATE_INFORMATION=true \
  SANITIZE=true \
  SORT_ORDER=queryname
picard -Xmx32g RevertSam \
  I=bams/E_coli_H10407_1.bam O=E_coli_H10407_1.unaligned.bam \
  REMOVE_ALIGNMENT_INFORMATION=true \
  RESTORE_ORIGINAL_QUALITIES=true \
  REMOVE_DUPLICATE_INFORMATION=true \
  SANITIZE=true \
  SORT_ORDER=queryname
