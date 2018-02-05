
task split_chunk {
	File bgen
	## index file is not used but is specified here so that it gets localized.
	File indexfile
	String chr
  Int start
  Int stop

	String? sge_queue
	String memory

	command {
	    bgenix -g ${bgen} -incl-range ${chr}:${start}-${stop} > ${basename(bgen)}_${chr}_${start}_${stop}.bgen;
			bgenix -g ${basename(bgen)}_${chr}_${start}_${stop}.bgen -index
	}

  runtime {
    memory: "${memory}"
		sge_queue: "${sge_queue}"
  }

	output {
		File outchunk="${basename(bgen)}_${chr}_${start}_${stop}.bgen"
		File outchunkindex="${basename(bgen)}_${chr}_${start}_${stop}.bgen.bgi"
	}

}

workflow split_to_chunks {

	File input_blocks
	String? sge_queue
	String memory

	Array[Array[String]] splits= read_tsv(input_blocks)

	scatter( split in splits) {
		call split_chunk {
			input: bgen=split[0],
			indexfile = split[0]+".bgi",
			chr=split[1],
			start=split[2],
      stop=split[3],
			sge_queue=sge_queue,
			memory=memory
		}
	}
}
