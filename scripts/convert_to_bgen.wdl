
task convert_file {
	File vcf
	String genofield
	String ofiletype
	Int precision
	Float input_rounding_error
	Int local_disk

	command {
		qctool -g ${vcf} -vcf-genotype-field ${genofield}  -filetype vcf -og ${basename(vcf)}.bgen -bgen-compression zlib -ofiletype ${ofiletype} -precision ${precision} -bgen-permitted-input-rounding-error ${input_rounding_error}
    bgenix -g ${basename(vcf)}.bgen -index
	}

	runtime {
		memory: "1G"
		sge_queue: "broad"
		disks: "local-disk ${local_disk} SSD"
	}

	output {
		File outbgen = "${basename(vcf)}.bgen"
		File index = "${basename(vcf)}.bgen.bgi"
	}

}

workflow convert_to_bgen {
	File files_to_conv
	String genofield ="GP"
	String ofiletype="bgen_v1.2"
	Int precision=8
  Float input_rounding_error=0.005
  ## disk size in GB.
  Int local_disk = 50

	Array[String] files = read_lines(files_to_conv)

	scatter( file in files) {

		call convert_file {
			input: vcf=file,
			genofield=genofield,
      ofiletype=ofiletype,
      precision=precision,
      input_rounding_error=input_rounding_error,
      local_disk=local_disk
		}
	}
}
