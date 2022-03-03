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
  dep_task :ARGO_sanger, ARGO, "sanger-wgs-variant-calling", :tumour_aln_cram => :indexed_BAM, :normal_aln_cram => :indexed_BAM_normal do |sample,options,dependencies|
    options = Sample.add_sample_options sample, options

    metadata = dependencies.flatten.select{|d| d.task_name.to_s == "ARGO_metadata" }.first
    options[:tumour_aln_metadata] = metadata.file('tumor.json')
    options[:normal_aln_metadata] = metadata.file('normal.json')


    bam = dependencies.flatten.select{|d| d.task_name.to_s == "indexed_BAM" }.first
    bam_normal = dependencies.flatten.select{|d| d.task_name.to_s == "indexed_BAM_normal" }.first

    options[:tumour_aln_cram] = bam.file('index')[bam.name + '.*am']
    options[:normal_aln_cram] = bam_normal.file('index')[bam_normal.name + '.*am']

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

    #tar_gz = Rbbt.var.Sanger.tarized_genomes["hg38.tar.gz"].find
    #if ! File.exists? tar_gz
    #  Open.mkdir File.dirname(tar_gz)
    #  Misc.in_dir File.dirname(options["ref_genome_fa"]) do
    #    CMD.cmd_log("tar cvfz #{tar_gz} .")
    #  end
    #end

    options["ref_genome_fa-generateBas"] = options["ref_genome_fa"]

    {:inputs => options}
  end

end
