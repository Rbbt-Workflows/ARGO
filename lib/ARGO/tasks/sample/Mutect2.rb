module Sample

  dep :ARGO_metadata
  dep :indexed_BAM do |sample,options,dependencies|
    options = Sample.add_sample_options sample, options
    sample_files = Sample.sample_files sample
    if sample_files && sample_files.include?("CRAM")
      options["Sample#ARGO_BAM"] = sample_files["CRAM"].first
      options["not_overriden"] = true
    end
    {:inputs => options, :jobname => sample}
  end
  dep :indexed_BAM_normal do |sample,options,dependencies|
    options = Sample.add_sample_options sample, options
    normal_sample = Sample.matched_normal sample, Sample.sample_study(sample)
    sample_files = Sample.sample_files normal_sample
    if sample_files && sample_files.include?("CRAM")
      options["Sample#ARGO_BAM"] = sample_files["CRAM"].first
      options["not_overriden"] = true
    end
    {:inputs => options, :jobname => sample}
  end
  dep_task :ARGO_mutect2_pre, ARGO, "gatk-mutect2-variant-calling", :tumour_aln_cram => :indexed_BAM, :normal_aln_cram => :indexed_BAM_normal do |sample,options,dependencies|
    options = Sample.add_sample_options sample, options

    metadata = dependencies.flatten.select{|d| d.task_name.to_s == "ARGO_metadata" }.first
    options[:tumour_aln_metadata] = metadata.file('tumor.json')
    options[:normal_aln_metadata] = metadata.file('normal.json')


    bam = dependencies.flatten.select{|d| d.task_name.to_s == "indexed_BAM" }.first
    bam_normal = dependencies.flatten.select{|d| d.task_name.to_s == "indexed_BAM_normal" }.first

    options[:tumour_aln_cram] = bam.file(File.join('index', File.basename(bam.dependencies.first.path)))
    options[:normal_aln_cram] = bam_normal.file(File.join('index', File.basename(bam_normal.dependencies.first.path)))

    options[:tumour_extra_info] = metadata.file('tumor_extra.tsv')
    options[:normal_extra_info] = metadata.file('normal_extra.tsv')

    options["publish_dir"] = "JOB_FILE:publish_dir"

    reference = options["reference"] || 'hg38'
    options["ref_fa"] = ARGO.link_FASTA(BWA.prepare_FASTA(HTS.helpers[:reference_file].call(reference)))
    options["ref_fa"] = options["ref_fa"].sub(/\.gz$/,'')

    panel_of_normals = options[:panel_of_normals] ||= case reference.to_s
          when 'b37'
            Organism["Hsa"].b37.known_sites["panel_of_normals.vcf"].produce.find
          when 'hg38'
            Organism["Hsa"].hg38.known_sites["panel_of_normals.vcf"].produce.find
          else
            raise "No default or Broad panel of normals for reference '#{reference}'"
          end

    panel_of_normals = nil if panel_of_normals == 'none'
    options[:panel_of_normals] = GATK.prepare_VCF panel_of_normals if panel_of_normals

    germline_resource ||= options[:germline_resource_vcfs] || :gnomad
    germline_resource = HTS.helpers[:vcf_file].call( reference, germline_resource) if germline_resource and ! File.exists?(germline_resource.to_s)
    options[:germline_resource_vcfs] = GATK.prepare_VCF germline_resource if germline_resource

    pileup_germline_resource = options[:contamination_variants] || :small_exac
    pileup_germline_resource = HTS.helpers[:vcf_file].call(reference, pileup_germline_resource) if pileup_germline_resource and ! File.exists?(pileup_germline_resource.to_s)
    options[:contamination_variants] = GATK.prepare_VCF pileup_germline_resource

    basedir = Rbbt.modules["gatk-mutect2-variant-calling"].find(:lib)
    options["mutect2_scatter_interval_files"] ||= File.join(basedir, "/assets/mutect2.scatter_by_chr/chr*.interval_list")

    options[:perform_bqsr] = "false" if options[:perform_bqsr].nil?

    {:inputs => options}
  end

  dep :ARGO_mutect2_pre
  extension "vcf"
  task :ARGO_mutect2_filters => :text do
    job = step(:ARGO_mutect2_pre).join
    vcf_files = job.file("publish_dir").glob("*/*/*.vcf.gz")
    first = vcf_files.shift
    Open.write(self.tmp_path) do |f|
      first.open do |sin|
        Misc.consume_stream(sin, false, f, false)
      end
      vcf_files.each do |v|
        io = CMD.cmd("zcat #{v} | grep -v ^# ", :pipe => true)
        Misc.consume_stream(io, false, f)
      end
    end
    nil
  end

  dep :ARGO_mutect2_filters
  extension "vcf"
  task :ARGO_mutect2 => :text do
    TSV.traverse step(:ARGO_mutect2_filters), :into => :stream, :type => :array do |line|
      next line if line[0] =~ /^#/
      
      chr = line.split("\t").first
      next unless chr =~ /^(chr)?[0-9MTXY]+$/
      next unless line =~ /PASS/

      line
    end
  end

end
