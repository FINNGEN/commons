
task index {
	File inputfile

	command {
		tabix -o vcf ${inputfile}
	}

	output {
		File out = "${inputfile}.tbi"
	}

}

workflow index_vcf {
	File filelist

	Array[String] files = read_tsv(filelist)

	scatter(for file in files) {
		 call index_cmd  {
		 	input: inputfile=file
		 }
	}

}
