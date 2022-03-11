module ARGO
  module Sanger

    DATA_URLS=<<-EOF.split("\n")
core_ref_GRCh38_hla_decoy_ebv.tar.gz
VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz
SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz
CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz
qcGenotype_GRCh38_hla_decoy_ebv.tar.gz
    EOF

    DATA_URLS.each do |file|
      url = File.join("http://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/", file)
      Rbbt.claim Rbbt.share.databases.Sanger[file], :url, url
    end
    
  end
end

module Sample

  dep :ARGO_metadata
  dep :indexed_BAM do |sample,options,dependencies|
    options = Sample.add_sample_options sample, options
    begin
      sample_files = Sample.sample_files sample
      if sample_files && sample_files.include?("CRAM")
        options["Sample#ARGO_BAM"] = sample_files["CRAM"].first
        options["not_overriden"] = true
      end
    rescue
    end
    {:inputs => options, :jobname => sample}
  end
  dep :indexed_BAM_normal do |sample,options,dependencies|
    options = Sample.add_sample_options sample, options
    begin
      normal_sample = Sample.matched_normal sample, Sample.sample_study(sample)
      sample_files = Sample.sample_files normal_sample
      if sample_files && sample_files.include?("CRAM")
        options["Sample#ARGO_BAM"] = sample_files["CRAM"].first
        options["not_overriden"] = true
      end
    rescue
    end
    {:inputs => options, :jobname => sample}
  end
  input :type_of_sequencing, :select, "Whole genome or whole exome", nil, :select_options => %w(WGS WES panel)
  dep_task :ARGO_sanger_pre, ARGO, "sanger-wgs-variant-calling", :tumour_aln_cram => :indexed_BAM, :normal_aln_cram => :indexed_BAM_normal do |sample,options,dependencies|
    options = Sample.add_sample_options sample, options

    metadata = dependencies.flatten.select{|d| d.task_name.to_s == "ARGO_metadata" }.first
    options[:tumour_aln_metadata] = metadata.file('tumor.json')
    options[:normal_aln_metadata] = metadata.file('normal.json')

    bam_normal = dependencies.flatten.select{|d| d.task_name.to_s == "indexed_BAM_normal" }.first
    bam_normal_name = bam_normal.step(:indexed_BAM).name
    bam = [dependencies.flatten.select{|d| d.task_name.to_s == "indexed_BAM" } - [bam_normal]].flatten.first

    extension = '.cram'
    options[:tumour_aln_cram] = bam.file('index')[bam.name + extension]
    options[:normal_aln_cram] = bam_normal.file('index')[bam_normal_name + extension]

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

    options["ref_genome_fa"] = ARGO.link_FASTA(BWA.prepare_FASTA(HTS.helpers[:reference_file].call(reference)))
    options["ref_genome_fa"] = options["ref_genome_fa"].sub(/\.gz$/,'')
    options["ref_genome_fa"] = options["ref_genome_fa"].sub(/\.gz$/,'')


    options["ref_genome_tar-sangerWgsVariantCaller"] = Rbbt.share.databases.Sanger["core_ref_GRCh38_hla_decoy_ebv.tar.gz"].produce.find
    options["ref_snv_indel_tar-sangerWgsVariantCaller"] = Rbbt.share.databases.Sanger["SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz"].produce.find
    options["vagrent_annot-sangerWgsVariantCaller"] = Rbbt.share.databases.Sanger["VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz"].produce.find
    options["ref_cnv_sv_tar-sangerWgsVariantCaller"] = Rbbt.share.databases.Sanger["CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz"].produce.find
    options["qcset_tar-sangerWgsVariantCaller"] = Rbbt.share.databases.Sanger["qcGenotype_GRCh38_hla_decoy_ebv.tar.gz"].produce.find

    options["ref_genome_tar-sangerWxsVariantCaller"] = Rbbt.share.databases.Sanger["core_ref_GRCh38_hla_decoy_ebv.tar.gz"].produce.find
    options["ref_snv_indel_tar-sangerWxsVariantCaller"] = Rbbt.share.databases.Sanger["SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz"].produce.find
    options["vagrent_annot-sangerWxsVariantCaller"] = Rbbt.share.databases.Sanger["VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz"].produce.find
    options["ref_cnv_sv_tar-sangerWxsVariantCaller"] = Rbbt.share.databases.Sanger["CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz"].produce.find
    options["qcset_tar-sangerWxsVariantCaller"] = Rbbt.share.databases.Sanger["qcGenotype_GRCh38_hla_decoy_ebv.tar.gz"].produce.find
    options["ref_genome_fa-generateBas"] = options["ref_genome_fa"]

    type_of_sequencing = options[:type_of_sequencing] || "WGS"

    case type_of_sequencing.to_s.upcase
    when "WGS"
      {:inputs => options}
    else
      {:task =>  "sanger-wxs-variant-calling", :inputs => options}
    end
  end

  dep :ARGO_sanger_pre
  extension "vcf"
  task :ARGO_sanger_filters => :text do
    job = step(:ARGO_sanger_pre).join
    vcf_files = job.file("publish_dir").glob("*/*/*.vcf.gz")
    vcf = vcf_files.shift
    while new = vcf_files.shift
      Open.write(self.tmp_path, HTS.job(:join_vcfs, nil, :vcf1 => vcf, :vcf2 => new).exec)
    end
    nil
  end

  dep :ARGO_sanger_filters
  extension "vcf"
  task :ARGO_sanger => :text do
    TSV.traverse step(:ARGO_sanger_filters), :into => :stream, :type => :array do |line|
      next line if line[0] =~ /^#/
      
      chr = line.split("\t").first
      next unless chr =~ /^(chr)?[0-9MTXY]+$/
      next unless line =~ /PASS/

      line
    end
  end


end
