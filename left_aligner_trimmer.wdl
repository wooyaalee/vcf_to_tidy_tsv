task LeftAlignAndTrimTask{
    File inputVCF
    File? inputTBI
    File inputFA
    File inputFAI
    File inputFADICT

    Int diskGB

    String outbase = basename(basename(inputVCF, ".gz"), ".vcf")

    command{
        gatk-lite -T LeftAlignAndTrimVariants \
        -V ${inputVCF} \
        -R ${inputFA} \
        -o ${outbase}.leftAlignedAndTrimmed.vcf
    }

    runtime{
        docker : "erictdawson/gatk:3.8-1-0"
        cpu : 1
        memory : "3 GB"
        disks : "local-disk " + diskGB + " HDD"
        preemptible : 3
    }

    output{
        File alignedAndTrimmedVCF = "${outbase}.leftAlignedAndTrimmed.vcf"
    }


}

workflow GATKLeftAlignAndTrim{
    File inputVCF
    File? inputTBI
    File inputFA
    File inputFAI
    File inputFADICT

    Int diskGB = ceil(size(inputVCF, "GB") + size(inputTBI, "GB")) + 50

    call LeftAlignAndTrimTask{
        input:
            inputVCF=inputVCF,
            inputTBI=inputTBI,
            inputFA=inputFA,
            inputFAI=inputFAI,
            inputFADICT=inputFADICT,
            diskGB=diskGB
    }
    
}
