module Sample

  input :sample_type, :select, "Sample type (Tumor or Normal) to align", "Tumor", :select_options => %w(Tumor Normal)
  input :reference, :select, "Reference code", nil, :select_options => %w(b37 hg38 mm10), :nofile => true
  dep :ARGO_metadata
  dep_task :ARGO_BAM, ARGO, "dna-seq-alignment", :analysis_metadata => :placeholder do |sample,options,dependencies|
    options = Sample.add_sample_options sample, options
    sample_type = options[:sample_type]

    if sample_type == "Normal"
      study = Sample.sample_study sample
      normal_sample = Sample.matched_normal sample, study
      sample_files = Sample.sample_files normal_sample
    else
      sample_files = Sample.sample_files sample
    end

    fastq_files = sample_files["FASTQ"] || sample_files["uBAM"]

    metadata_job = dependencies.flatten.select{|dep| dep.task_name.to_s == "ARGO_metadata" }.first

    options["publish_dir"] = "JOB_FILE:publish_dir"

    options["analysis_metadata"] = metadata_job.file(sample_type.downcase + '.json')
    options["sequencing_files"] = fastq_files.first.first.reverse.sub(/1/, '*').reverse

    reference = options["reference"] || 'hg38'
    options["ref_genome_fa"] = ARGO.link_FASTA(BWA.prepare_FASTA(HTS.helpers[:reference_file].call(reference)))
    options["ref_genome_fa"] = options["ref_genome_fa"].sub(/\.gz$/,'')

    {:inputs => options}
  end

  dep_task :ARGO_BAM_normal, Sample, :ARGO_BAM, :sample_type => "Normal"

  input :use_rbbt_aligner, :boolean, "Use rbbt aligner instead of ARGO", false
  dep :ARGO_BAM do |sample,options,dependencies|
    if options[:use_rbbt_aligner]
      {:task => :BAM, :jobname => sample}
    else
      {:inputs => options, :jobname => sample}
    end
  end
  task :indexed_BAM => :string do
    bam = dependencies.first
    bam = bam.path if Step === bam
    raise "File does not exist" unless File.exists?(bam)
    output = file('index')
    Samtools.prepare_BAM bam, output
    output.glob("*").each do |file|
      basename, _sep, extension = File.basename(file).partition(".")
      fixed_name = output[[self.name, extension] * "."]
      Misc.in_dir output do
        `ln -s #{File.basename(file)} #{fixed_name}`
      end
    end
    File.basename(self.path)
  end

  input :use_rbbt_aligner, :boolean, "Use rbbt aligner instead of ARGO", false
  dep :ARGO_BAM_normal do |sample,options,dependencies|
    if options[:use_rbbt_aligner]
      {:task => :BAM_normal, :jobname => sample}
    else
      {:inputs => options, :jobname => sample}
    end
  end
  task :indexed_BAM_normal => :string do
    bam = dependencies.first
    bam = bam.path if Step === bam
    raise "File does not exist" unless File.exists?(bam)
    output = file('index')
    Samtools.prepare_BAM bam, output
    output.glob("*").each do |file|
      basename, _sep, extension = File.basename(file).partition(".")
      fixed_name = output[[self.name, extension] * "."]
      Misc.in_dir output do
        `ln -s #{File.basename(file)} #{fixed_name}`
      end
    end
    File.basename(self.path)
  end

end
