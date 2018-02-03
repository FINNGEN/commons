
task split_chunk {
	File bgen
	String chr
    Int start
    Int stop
    
	command {
	    bgenix -g ${bgen} -incl-range {$chr}:${start}-${stop} $ > ${basename(bgen)}_${chr}_${start}_${stop}.bgen
	}

    runtime {
        memory: "1G"
        sge_queue: "broad"
    }
	
	output {
		File outchunk="${basename(bgen)}_${chr}_${start}_${stop}.bgen"
	}

}

workflow split_to_chunks {

	File input_blocks
	
	Array[Array[String]] splits= read_tsv(input_blocks)

	scatter( split in splits) {
		call split_chunk {
			input: bgen=split[0],
			chr=split[1],
			start=split[2],
            stop=split[3]
		}
	} 
}
