gatk=$1
${gatk} \
       LeftAlignAndTrimVariants \
        -R $2 \
        --variant $3 \
        -O $(basename $(basename $3 .gz) .vcf).leftAlignedAndTrimmed.vcf
