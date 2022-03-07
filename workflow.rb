require 'rbbt-util'
require 'rbbt/workflow'
require 'rbbt/workflow/integration/nextflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/MODULE'

Workflow.require_workflow "HTS"

module ARGO
  extend Workflow

 def self.link_FASTA(file)
    dir = File.dirname(file)

    Dir.glob(File.join(dir, '*.gz.*')).each do |source|
      target = source.sub('.gz','')
      next if File.exists?(target)
      Open.link source, target 
    end

    uncompressed = file.sub('.gz', '')

    if ! File.exists?(uncompressed)
      CMD.cmd("zcat '#{file}' > #{uncompressed}")
    end

    file
  end

  extension :json
  self.nextflow Rbbt.modules["argo-qc-tools"].cutadapt.find(:lib), :text => "output_dir/*.json"

  extension :json
  self.nextflow Rbbt.modules["argo-qc-tools"].fastqc.find(:lib), :text => "qc_metrics.json"

  extension :cram
  self.nextflow Rbbt.modules["dna-seq-processing-wfs"]["dna-seq-alignment"].find(:lib), :binary => "grch38-aligned.merged.cram"

  extension "vcf.gz"
  self.nextflow Rbbt.modules["gatk-mutect2-variant-calling"]["gatk-mutect2-variant-calling"].find(:lib), :text => "publish_dir/M2_pGenVarSnv/out/*.vcf.gz"

  input "ref_genome_fa-generateBas", :file, "Ill defined Nextflow parameter", nil, :nofile => true
  extension "vcf.gz"
  self.nextflow Rbbt.modules["sanger-wgs-variant-calling"]["sanger-wgs-variant-calling"].find(:lib), :text => "publish_dir/M2_pGenVarSnv/out/*.vcf.gz"

  input "ref_genome_fa-generateBas", :file, "Ill defined Nextflow parameter", nil, :nofile => true
  self.nextflow Rbbt.modules["sanger-wxs-variant-calling"]["sanger-wxs-variant-calling"].find(:lib), :text => "publish_dir/M2_pGenVarSnv/out/*.vcf.gz"

end

require 'ARGO/tasks/sample.rb' if defined? Sample

#require 'rbbt/knowledge_base/MODULE'
#require 'rbbt/entity/MODULE'

