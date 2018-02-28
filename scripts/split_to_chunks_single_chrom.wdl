
task split_chunk {
	File infile
	File indexfile
	File splitfile

	String memory
	Int ncpu
	Int disksize

	command {
		split_single_chrom.py ${infile} ${splitfile}
	}

  runtime {
    memory: "${memory}"
	cpu: "${ncpu}"
	disks: "local-disk ${disksize} SSD"
  }

	output {
		Array[File] outchunks=glob("*.bgen")
	}

}


task split_ranges {
	String splits

	command {
			echo '${splits}' | tr ',' '\n' > splitfile
	}

	output {
		File out = "splitfile"
	}

}

workflow split_to_chunks {

	File input_blocks
	String memory
	### vcf or bgen
	String filetype
	Int ncpu=8
	Int disksize=50

	Array[Array[String]] splits= read_tsv(input_blocks)

	String indexsuffix = if filetype=="vcf" then ".tbi" else ".bgi"


	scatter( split in splits) {

		call split_ranges {
			input: splits=split[1]
		}

		call split_chunk {
			input: infile=split[0],
			indexfile = split[0]+indexsuffix,
			splitfile=split_ranges.out,
			ncpu=ncpu,
			disksize=disksize,
			memory=memory
		}

	}
}
