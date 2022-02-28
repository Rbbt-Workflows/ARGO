module Sample
  module ARGOMetadata
    def self.analysis(analysis, study)
      JSON.parse <<-EOF
{
  "analysisId": "#{analysis}",
  "studyId": "#{study}",
  "analysisState": "Rbbt",
  "analysisType": {
    "name": "sequencing_experiment",
    "version": 10
  },
  "workflow": {
    "name": "sequencing-data-submission",
    "run_id": "nasty_montalcini",
    "version": "0.4.0",
    "short_name": "seq-submission"
  },
  "experiment": {
    "platform": "ILLUMINA",
    "platform_model": "HiSeq 2000",
    "sequencing_date": "2014-12-12",
    "experimental_strategy": "WGS",
    "sequencing_center": "EXT"
  },
  "read_group_count": 1,
  "submitter_sequencing_experiment_id": "#{analysis}"
}
      EOF
    end

    def self.donor(donor, study, gender = "Unkown")
      JSON.parse <<-EOF
{
  "donorId": "#{donor}",
  "submitterDonorId": "#{donor}",
  "studyId": "#{study}",
  "donorGender": "#{gender}",
  "info": {}
}
      EOF
    end

    def self.sample(study, donor, sample, specimen, specimen_type)
      sample = JSON.parse <<-EOF
{
  "info": {},
  "sampleId": "#{sample}",
  "specimenId": "#{specimen}",
  "submitterSampleId": "#{sample}",
  "sampleType": "DNA",
  "specimen": {
    "info": {},
    "specimenId": "#{specimen}",
    "donorId": "#{donor}",
    "submitterSpecimenId": "#{specimen}",
    "tumourNormalDesignation": "#{specimen_type}",
    "specimenTissueSource": "Solid tissue",
    "specimenType": "#{specimen_type}"
    }
}
      EOF

      sample["donor"] = self.donor(donor, study)

      sample
    end

    def self.file(study, analysis, filename, md5 = "0001")
      md5 = CMD.cmd('md5sum', filename).read.strip.split(" ").first

      JSON.parse <<-EOF
{
  "info": {},
  "objectId": "#{[study, analysis, filename] * ":"}",
  "studyId": "#{study}",
  "analysisId": "#{analysis}",
  "fileName": "#{File.basename(filename)}",
  "fileSize": #{File.size(filename)},
  "fileType": "FASTQ",
  "fileMd5sum": "#{md5}",
  "fileAccess": "controlled"
}
      EOF
    end

    def self.read_group(read_group_name, fastq1, fastq2, insert_size, library_name, read_length)
      JSON.parse(<<-EOF)
{
  "file_r1": "#{File.basename(fastq1)}",
  "file_r2": "#{File.basename(fastq2)}",
  "insert_size": #{insert_size},
  "library_name": "#{library_name}",
  "is_paired_end": true,
  "platform_unit": "74_8b",
  "read_length_r1": #{read_length},
  "read_length_r2": #{read_length},
  "sample_barcode": null,
  "submitter_read_group_id": "#{read_group_name}"
}
    EOF
    end

    def self.common_prefix(a)
      return "" unless (a.length > 0)
      result = 0
      (0 ... a.first.length).each do |k|
        all_matched = true
        character = a.first[k]
        a.each{ |str| all_matched &= (character == str[k]) }
        break unless all_matched
        result+=1
      end
      a.first.slice(0,result)
    end


    def self.read_groups(study, sample)
      sample_files = Sample.sample_files(sample)

      Misc.zip_fields(sample_files["FASTQ"] || []).collect do |f1,f2|
        name = common_prefix([File.basename(f1), File.basename(f2)])
        name = sample if name.empty?

        library_name, _sep, name = name.partition(".")
        name = library_name if name.nil? || name.empty?

        insert_size = 283
        read_length = 150

        read_group(name, f1, f2, insert_size, library_name, read_length)
      end
    end

    def self.alignment(study, sample, type)
      analysis_code = "Rbbt"
      analysis = Sample::ARGOMetadata.analysis(analysis_code, study)

      sample_files = Sample.sample_files sample

      analysis["samples"] = [ Sample::ARGOMetadata.sample(study, sample, sample, sample, type) ]

      analysis["files"] = sample_files["FASTQ"].collect do |fastqs|
        fastqs.collect do |fastq|
          Sample::ARGOMetadata.file(study, analysis_code, fastq)
        end
      end.flatten if sample_files["FASTQ"]

      analysis["read_groups"] = Sample::ARGOMetadata.read_groups(study, sample)
      analysis["read_group_count"] = analysis["read_groups"].length

      analysis.to_json
    end

  end

  input :sequencing_center, :string, "SEQUENCING_CENTER BAM field", "DefaultSequencingCenter"
  input :study, :string, "Study code", "Rbbt"
  extension :json
  task :ARGO_metadata => :text do |sequencing_center,study|
    analysis_code = "Rbbt"
    analysis = Sample::ARGOMetadata.analysis(analysis_code, study)

    study = Sample.sample_study sample
    normal_sample = Sample.matched_normal sample, study

    sample_files = Sample.sample_files sample
    normal_sample_files = Sample.sample_files normal_sample if normal_sample

    analysis["samples"] = [Sample::ARGOMetadata.sample(study, sample, sample, sample, "Tumor")]

    analysis["samples"] += [Sample::ARGOMetadata.sample(study, normal_sample, normal_sample, normal_sample, "Normal")] if normal_sample

    analysis["files"] = sample_files["FASTQ"].collect do |fastqs|
      fastqs.collect do |fastq|
        Sample::ARGOMetadata.file(study, analysis_code, fastq)
      end
    end.flatten if sample_files["FASTQ"]

    analysis["files"] += normal_sample_files["FASTQ"].collect do |fastqs|
      fastqs.collect do |fastq|
        Sample::ARGOMetadata.file(study, analysis_code, fastq)
      end
    end.flatten if normal_sample_files && normal_sample_files["FASTQ"]

    analysis["read_groups"] = Sample::ARGOMetadata.read_groups(study, sample)
    analysis["read_groups"] += Sample::ARGOMetadata.read_groups(study, normal_sample) if normal_sample
    analysis["read_group_count"] = analysis["read_groups"].length

    Open.write(file('tumor.json'), Sample::ARGOMetadata.alignment(study, sample, 'Tumor'))
    Open.write(file('normal.json'), Sample::ARGOMetadata.alignment(study, normal_sample, 'Normal')) if normal_sample

    tsv = TSV.setup({}, "type~submitter_id,uniform_id#:type=:single")
    tsv["donor"] = [sample, sample]
    tsv["specimen"] = [sample, sample]
    tsv["sample"] = [sample, sample]

    Open.write(file('tumor_extra.tsv'), tsv.to_s(:header_hash => "", :preamble => ""))

    tsv["donor"] = [normal_sample, normal_sample]
    tsv["specimen"] = [normal_sample, normal_sample]
    tsv["sample"] = [normal_sample, normal_sample]
    Open.write(file('normal_extra.tsv'), tsv.to_s(:header_hash => "", :preamble => ""))

    analysis.to_json
  end

end

require 'ARGO/tasks/sample/BAM'
require 'ARGO/tasks/sample/Mutect2'
require 'ARGO/tasks/sample/Sanger'
