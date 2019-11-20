
task index {
	File inputfile
	Int disksize

	runtime {
		cpu:"2"
		disks: "local-disk 100 SSD"
	}

	command {
		tabix -p vcf ${inputfile}
	}

	output {
		File out = "${inputfile}.tbi"
	}

}

workflow index_vcf {
	File filelist
	Int disksize=100

	Array[String] files = read_lines(filelist)

	scatter(file in files) {

		 call index {
		 	input: inputfile=file,
		 	disksize = disksize
		 }
	}

}
