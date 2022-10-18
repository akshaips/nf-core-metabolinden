#!/usr/bin/env nextflow
/*



========================================================================================
                         nf-core/metabolinden
========================================================================================
 nf-core/metabolinden Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/metabolinden
----------------------------------------------------------------------------------------
*/

log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+
def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/metabolinden --input '*.mzML' -profile docker"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}
////////////////////////////////////////////////////
/* --     randomize files     -- */
////////////////////////////////////////////////////

def group_files( input,samples_in_chunks=2,randomize=false,seed=42){
old_key=input[0]
input=input[1]

if(samples_in_chunks<=1 || samples_in_chunks>=input.size())
{
  step_wise_linking=false
  return input
}

def range_of_vec = (0..input.size()-1).toList()
if(randomize==true){
  range_of_vec.shuffle(new Random(seed))
}



if((range_of_vec.size()%samples_in_chunks)<2)
{
  while((range_of_vec.size()%samples_in_chunks)<2)
  {
    samples_in_chunks = (samples_in_chunks+1)
  }
  println 'The minimum number of samples must be one. Adjusting samples_in_chunks to '+ (samples_in_chunks)


}

chunks=range_of_vec.collate(samples_in_chunks)
def vector_sizes=[]
groups=[]
for(n in 0..chunks.size()-1)
{
  vector_sizes.add(chunks[n].size())
  aa=["chunk"+(n+1)+"_chunkend_"+old_key]*chunks[n].size()
  for(i in chunks[n])
  {
    groups.add(["chunk"+(n+1)+old_key,old_key,file(input[i])])
  }
}
return(groups)


}

////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////
Channel.fromPath(params.input,checkIfExists: true)
.map{def key = "start"
return tuple(key, it)}
 .set { mzml_input}//.map { tag, files -> tuple( groupKey(tag, files.size()), files ) }
//.transpose()

/*
* Create a channel for centroiding parameters
*/

if(params.need_centroiding==true)
{
Channel.fromPath(params.peak_picker_param,checkIfExists: true)
.set { peak_picker_param }
}
/*
* Create a channel for recalibration parameters
*/

if(params.need_recalibration==true)
{
Channel.fromPath(params.peak_recalibration_param,checkIfExists: true)
.into { peak_recalibration_param;c3 }

/*
* Create a channel for recalibration standards
*/

Channel.fromPath(params.recalibration_masses,checkIfExists: true)
.into { recalibration_masses;c4 }
}

/*
* Create a channel for feature_detection
*/

if(params.need_quantification==true)
{
Channel.fromPath(params.feature_finder_param,checkIfExists: true)
.set { feature_finder_param }
}else{
  params.need_alignment=false
  params.need_linking=false
  params.need_identification=false

}
/*
* Create a channel for feature aligment parameters
*/
if(params.need_alignment==true)
{
Channel.fromPath(params.feature_alignment_param,checkIfExists: true)
.set { feature_alignment_param }
}

/*
* Create a channel for feature linker parameters
*/

if(params.need_linking==true)
{
  Channel.fromPath(params.feature_linker_param,checkIfExists: true)
.set { feature_linker_param }

if(params.step_wise_linking==true)
if(params.use_same_setting_for_second_linking==false)
{
  Channel.fromPath(params.feature_linker_param2,checkIfExists: true)
  .set { feature_linker_param2 }
}



}

/*
* Create a channel for feature linker parameters
*/

if(params.need_identification==true)
{
  Channel.fromPath(params.identification_input,checkIfExists: true)
.into { identification_input;identification_input_qc }
}

if(params.need_qc==true)
{
  Channel.fromPath(params.qc_file,checkIfExists: true)
.set { experimental_design }
}

/*
* Create a channel for feature linker parameters
*/


// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, 'Specify correct --awsqueue and --awsregion parameters on AWSBatch!'
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, 'Outdir not on S3 - specify S3 Bucket to run on AWSBatch!'
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, 'Specify a local tracedir or run without trace! S3 cannot be used for tracefiles.'
}

// Stage config files
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)
ch_output_qc_rmd = channel.fromPath("$projectDir/assets/qc.Rmd", checkIfExists: true)


////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

// Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = workflow.runName
// TODO nf-core: Report custom parameters here
summary['Input']            = params.input
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
}

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-metabolinden-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/metabolinden Workflow Summary'
    section_href: 'https://github.com/nf-core/metabolinden'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf('.csv') > 0) filename
                      else null
        }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file 'software_versions.csv'

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    OpenMSInfo |  grep -oP -m 1 '([0-9][.][0-9][.][0-9])' &> v_openms.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}


///////////////////////////////////////////////////////////
/* --         main functions of the workflow          -- */
//////////////////////////////////////////////////////////

/*
 * Step 1.  Do centroiding if needed
 */

if(params.need_centroiding==true){
  process process_peak_picker_openms  {
    label 'openms'
    //label 'process_low'
    tag "${mzMLFile} using ${setting_file} parameter"
    publishDir "${params.outdir}/process_peak_picker_pos_openms", mode: params.publish_dir_mode, enabled: params.publishDir_intermediate

    input:
    set val(key), file(mzMLFile) from mzml_input
    each file(setting_file) from peak_picker_param

    output:
    tuple val("${key}_${setting_file.baseName}"), file("output_${key}_${setting_file.baseName}/${mzMLFile}") into recalibration_channel

    script:
    """
    mkdir "output_${key}_${setting_file.baseName}"
    PeakPickerHiRes -in $mzMLFile -out "output_${key}_${setting_file.baseName}/$mzMLFile" -ini $setting_file
    """
  }

  }else{
    recalibration_channel=mzml_input
    log.info "skipping centroiding!"
  }

  /*
   * Step 2.  Do recalibration if needed
   */



  if(params.need_recalibration==true){
    process process_recalibration_openms  {
      label 'openms'
      //label 'process_low'
      tag "processing ${mzMLFile} using ${setting_file}"
      publishDir "${params.outdir}/process_recalibration_openms", mode: params.publish_dir_mode, enabled: params.publishDir_intermediate

      input:
      set val(key), file(mzMLFile) from recalibration_channel
      each file(setting_file) from peak_recalibration_param
      each file(recal_masses) from recalibration_masses


      output:
      set val("${key}_${setting_file.baseName}"), file("output_${key}_${setting_file.baseName}/${mzMLFile}") into quant_feature_detection
      set val("${key}_${setting_file.baseName}"), file("output_${key}_${setting_file.baseName}/*.png") into recal_plots
      set val("${key}_${setting_file.baseName}"), file("output_${key}_${setting_file.baseName}/*.csv") into recal_info

      script:
      """
      mkdir "output_${key}_${setting_file.baseName}"

      InternalCalibration -in $mzMLFile -out "output_${key}_${setting_file.baseName}/$mzMLFile" \\
      -ini $setting_file -cal:lock_in $recal_masses -quality_control:models_plot "output_${key}_${setting_file.baseName}/${mzMLFile.baseName}_models_plot.png" \\
      -quality_control:residuals_plot "output_${key}_${setting_file.baseName}/${mzMLFile.baseName}_residuals_plot.png" \\
      -quality_control:residuals "output_${key}_${setting_file.baseName}/${mzMLFile.baseName}_residuals.csv" \\
      -quality_control:models "output_${key}_${setting_file.baseName}/${mzMLFile.baseName}_models.csv"
      """
    }
    }else{
      quant_feature_detection=recalibration_channel
    }


    /*
     * Step 3.  Do quantification if needed
     */

    if(params.need_quantification==true){
      process process_masstrace_detection_openms  {
        label 'openms'
        //label 'process_low'
        tag "processing ${mzMLFile} using ${setting_file}"
        publishDir "${params.outdir}/process_masstrace_detection_pos_openms", mode: params.publish_dir_mode, enabled: params.publishDir_intermediate

        input:
        set val(key), file(mzMLFile) from quant_feature_detection
        each file(setting_file) from feature_finder_param

        output:
        set val("${key}_${setting_file.baseName}"), file("output_${key}_${setting_file.baseName}/${mzMLFile}.featureXML") into quant_feature_alignment_tmp

        script:
        """
        mkdir "output_${key}_${setting_file.baseName}"
        FeatureFinderMetabo -in $mzMLFile -out "output_${key}_${setting_file.baseName}/${mzMLFile}.featureXML" -ini $setting_file
        """
      }
      }else{
        quant_feature_alignment_tmp=quant_feature_detection
      }

      /*
       * Step 4.  Do alignment if needed
       */

       quant_feature_alignment_tmp
       .groupTuple().set{alignment_input}


if(params.need_alignment==true)
{
  process process_masstrace_alignment_openms  {
    label 'openms'
    //label 'process_low'
    tag "$key"
    publishDir "${params.outdir}/process_masstrace_alignment_openms", mode: params.publish_dir_mode, enabled: params.publishDir_intermediate

    input:
    set val(key), file(mzMLFile) from alignment_input
    each file(setting_file) from feature_alignment_param

    output:
    set val("${key}_${setting_file.baseName}"), file("output_${key}_${setting_file.baseName}/*.*") into feature_linker_input_tmp

    script:
    def inputs_aggregated = mzMLFile.collect{ "$it" }.join(" ")
    def output_aggregated = mzMLFile.collect{ "\"output_${key}_${setting_file.baseName}/$it\"" }.join(" ")
    """
    mkdir "output_${key}_${setting_file.baseName}"
    MapAlignerPoseClustering -in $inputs_aggregated -out $output_aggregated -ini $setting_file
    """
  }
}else{

  feature_linker_input_tmp=alignment_input
}


/*
 * Step 5.  Do linking if needed
 */
feature_linker_input_tmp.set{feature_linker_input}


if(params.need_linking==true)
{
  if(params.step_wise_linking==false)
  {
    process process_masstrace_linker_openms  {
      label 'openms'
      //label 'process_low'
      tag "$setting_file"
      publishDir "${params.outdir}/process_masstrace_linker_openms", mode: params.publish_dir_mode, enabled: params.publishDir_intermediate

      input:
      set val(key), file(mzMLFile) from feature_linker_input
      each file(setting_file) from feature_linker_param

      output:
      set val("${key}_${setting_file.baseName}"), file("output_${key}_${setting_file.baseName}/*.*") into feature_identification_input

      script:
      def inputs_aggregated = mzMLFile.collect{ "$it" }.join(" ")
      """
      mkdir "output_${key}_${setting_file.baseName}"
      FeatureLinkerUnlabeledQT -in $inputs_aggregated -out "output_${key}_${setting_file.baseName}/linked_features.consensusXML" -ini $setting_file
      """
    }
  }else{

    feature_linker_input.set{chunker}


    process group_files_chunks {
  //label 'process_low'
      input:
      set val(key), file(mzMLFile) from chunker

      output:
      set val(a), val(key), file(otfile) into chunked_output


      script:
      def old_key=key
      def input = mzMLFile.collect()

      def samples_in_chunks=params.number_of_files
      if(samples_in_chunks<=1 || samples_in_chunks>=input.size())
      {
        step_wise_linking=false
        return input
      }

      def range_of_vec = (0..input.size()-1).toList()
      if(params.randomize_files==true){
        range_of_vec.shuffle(new Random(params.seed_for_linking))
      }



      if((range_of_vec.size()%samples_in_chunks)<2)
      {
        while((range_of_vec.size()%samples_in_chunks)<2)
        {
          samples_in_chunks = (samples_in_chunks+1)
        }
        println 'The minimum number of samples must be one. Adjusting samples_in_chunks to '+ (samples_in_chunks)


      }

      def chunks=range_of_vec.collate(samples_in_chunks)
      def vector_sizes=[]
      def groups=[]
      def ids=[]
      def outfiles=[]
      for(n in 0..chunks.size()-1)
      {
        vector_sizes.add(chunks[n].size())
        aa=["chunk"+(n+1)+"_chunkend_"+old_key]*chunks[n].size()
        for(i in chunks[n])
        {
          ids.add("chunk"+(n+1)+old_key)
          outfiles.add(input[i])
        }
      }
      a=ids
      otfile=outfiles

      """
      echo "dummy"
      """
    }



    chunked_output
         .transpose().groupTuple().set{chunked_input_for_linking}


    process process_masstrace_linker_openms_chunks  {
      label 'openms'
      //label 'process_low'
      tag "${key.unique().join("")}"
      publishDir "${params.outdir}/process_masstrace_linker_openms_chunks", mode: params.publish_dir_mode, enabled: params.publishDir_intermediate

      input:
      set val(key_group),val(key), file(mzMLFile) from chunked_input_for_linking
      each file(setting_file) from feature_linker_param

      output:
      set val("${key.unique().join("")}_${setting_file.baseName}"), file("output_${key_group}_${setting_file.baseName}/*.*"),file("setting_${key_group}_${setting_file.baseName}/*.ini") into feature_linking2_input_tmp

      script:
      def inputs_aggregated = mzMLFile.collect{ "$it" }.join(" ")
      """
      mkdir "output_${key_group}_${setting_file.baseName}"
      mkdir "setting_${key_group}_${setting_file.baseName}"
      FeatureLinkerUnlabeledQT -in $inputs_aggregated -out "output_${key_group}_${setting_file.baseName}/${key_group}_linked_features.consensusXML" -ini $setting_file
      cp ${setting_file} "setting_${key_group}_${setting_file.baseName}"

      """
    }
  feature_linking2_input_tmp.groupTuple(by:0).map{a,b,c->tuple(a,b,c[0])}.set{feature_linking2_input}

if(params.use_same_setting_for_second_linking==false)
{
  process process_masstrace_linker_openms_chunks_join  {
    label 'openms'
    //label 'process_low'
    tag "$st_file"
    publishDir "${params.outdir}/process_masstrace_linker_openms_chunks_join", mode: params.publish_dir_mode, enabled: params.publishDir_intermediate

    input:
    set val(key), file(mzMLFile),file(st_file) from feature_linking2_input
    each file(setting_file) from feature_linker_param2

    output:
    set val("${key}_${st_file.baseName}"), file("output_${key}_${st_file.baseName}/*.*") into feature_identification_input

    script:
    def inputs_aggregated = mzMLFile.collect{ "$it" }.join(" ")
    """
    mkdir "output_${key}_${st_file.baseName}"
    FeatureLinkerUnlabeledQT -in $inputs_aggregated -out "output_${key}_${st_file.baseName}/linked_features.consensusXML" -ini $setting_file -keep_subelements
    """
  }
}else{
  process process_masstrace_linker_openms_chunks_join_samesetting  {
    label 'openms'
    //label 'process_low'
    tag "${st_file.baseName}"
    publishDir "${params.outdir}/process_masstrace_linker_openms_chunks_join", mode: params.publish_dir_mode, enabled: params.publishDir_intermediate

    input:
    set val(key), file(mzMLFile),file(st_file) from feature_linking2_input

    output:
    set val("${key}_${st_file.baseName}"), file("output_${key}_${st_file.baseName}/*.*") into feature_identification_input

    script:
    def inputs_aggregated = mzMLFile.collect{ "$it" }.join(" ")
    def setting = st_file.collect()

    """
    mkdir "output_${key}_${st_file.baseName}"
    FeatureLinkerUnlabeledQT -in $inputs_aggregated -out "output_${key}_${st_file.baseName}/linked_features.consensusXML" -ini $st_file -keep_subelements
    """
  }
}


  }

}else{

  feature_identification_input=feature_linker_input
}

/*
 * Step 5.  Do identification if needed
 */
if(params.need_identification==true)
{
process convert_library_to_idXML  {
  label 'openms'
  //label 'process_low'
  tag "$libfile"
  publishDir "${params.outdir}/convert_library_to_idXML", mode: params.publish_dir_mode, enabled: params.publishDir_intermediate
  echo true
  input:
  file libfile from identification_input


  output:
  file "${libfile.baseName}/${libfile.baseName}.idXML" into libsearch_database

  """
  #!/usr/bin/env python3.9

  import pyopenms as py
  import csv
  import re
  import numpy as np
  import os
  import xml.etree.ElementTree as ET
  from pathlib import Path


  os.mkdir("${libfile.baseName}")
  inputs=list(filter((lambda x: re.search(r'tsv\$', x)), os.listdir()))
  if len(inputs)>1:
    sys.exit('Expect one tsv file')

  inputs="${libfile}"
  file_stem=Path(inputs).stem

  header=[]
  line_nr=0
  use_rt="$params.identification_use_rt"=="true"
  convert_to_seconds="$params.identification_convert_rt_to_seconds"=="true"

  fm=[]
  pm=[]
  id_trace=1
  with open(inputs) as tsvfile:
      tsvreader = csv.reader(tsvfile, delimiter="\t")
      for line in tsvreader:
          if line_nr==0:
            header=line
            all_mzs=[i for i, x in enumerate(header) if "$params.identification_mz_column" in x]
            all_rt=[i for i, x in enumerate(header) if "$params.identification_rt_column" in x]
            if(use_rt==True):
              rest=[i for i, x in enumerate(header) if "$params.identification_rt_column" not in x and "$params.identification_mz_column" not in x]
              rt_range=[0]
            else:
              rest=[i for i, x in enumerate(header) if "$params.identification_mz_column" not in x]
              rt_range=np.arange ($params.identification_min_rt, $params.identification_max_rt, $params.identification_scan_time)
              all_rt=[i for i, x in enumerate(header) if "$params.identification_rt_column" in x]
              #all_mzs=np.repeat(all_mzs,len(rt_range))
              convert_to_seconds=False
          else:
            mz=np.float(line[all_mzs[0]])
            rt_ex=np.float(line[all_rt[0]])
            for i in range(0,len(rt_range)):
              if(use_rt==False):
                rt_ex=np.float(rt_range[i])
              if(convert_to_seconds==True):
                rt_ex=rt_ex*60
              pid= py.PeptideIdentification()
              pid.setMZ(mz)
              pid.setRT(rt_ex);
              pid.setIdentifier(header[all_mzs[0]]+"_"+header[all_rt[0]]+"_"+str(mz)+"_"+str(rt_ex))
              #pid.setMetaValue("mz_adduct",header[all_mzs[i]])
              #pid.setMetaValue("cc_adduct",header[all_ccs[i]])
              #pid.setMetaValue("mz_adduct_value",line[all_mzs[i]])
              #pid.setMetaValue("cc_adduct_value",line[all_ccs[i]])
              pid.setIdentifier(str(id_trace))
              hit=py.PeptideHit()
              hit.setSequence(py.AASequence.fromString("METABOLITE"))
              phit=py.ProteinHit()
              phit.setAccession(str(id_trace))
              pe=py.PeptideEvidence()
              pe.setProteinAccession(str(id_trace));
              hit.addPeptideEvidence(pe)
              pid.insertHit(hit)

              pid.setHits([hit])
              ppp=py.ProteinIdentification()
              ppp.insertHit(phit)
              ppp.setIdentifier(str(id_trace))
              id_trace=id_trace+1
              for j in range(0,len(rest)):
                pid.setMetaValue(header[rest[j]],line[rest[j]])
              fm.append(pid)
              pm.append(ppp)
          line_nr=line_nr+1

  f=py.IdXMLFile()
  f.store("${libfile.baseName}/"+file_stem+".idXML",pm,fm)


  """
}


process process_masstrace_matchlib_openms  {
  label 'openms'
  //label 'process_low'
  tag "$input"
  publishDir "${params.outdir}/process_masstrace_matchlib_openms", mode: params.publish_dir_mode, enabled: params.publishDir_intermediate

  input:
  set val(key), file(input) from feature_identification_input
  each file(idfile) from libsearch_database


  output:
  set val("${key}_${idfile.baseName}"), file("output_${key}_${idfile.baseName}/$input") into qc_input

  """
  mkdir "output_${key}_${idfile.baseName}"
  IDMapper -in $input -id $idfile -out "output_${key}_${idfile.baseName}/$input" -mz_reference 'precursor' \\
  -mz_measure 'ppm' -rt_tolerance $params.internal_database_rt_tolerance -mz_tolerance $params.internal_database_ppm_tolerance
  """
}


}else{

qc_input=feature_identification_input
}

/*
 * Step 6.  Export the data
 */
 if(params.need_exporting==true)
 {
   process process_feature_exporter_openms  {
     label 'openms'
     //label 'process_low'
     tag "$consensusXML - $key"
     publishDir "${params.outdir}/process_feature_exporter_openms", mode: params.publish_dir_mode, enabled: params.publishDir_intermediate

     input:
     set val(key), file(consensusXML) from qc_input

     output:
     set val(key), file("output_${key}/${consensusXML.baseName}.tsv") into output_formatting

     """
     mkdir "output_${key}"
     TextExporter -in $consensusXML -out "output_${key}/${consensusXML.baseName}.tsv" -feature:add_metavalues 0 -id:add_metavalues 0
     """
   }

   process process_feature_output_formatting_openms  {
     label 'openms'
     //label 'process_low'
     tag "$consensusXML - $key"
     publishDir "${params.outdir}/process_feature_output_formatting_openms", mode: params.publish_dir_mode, enabled: params.publishDir_intermediate

     input:
     set val(key), file(consensusXML) from output_formatting

     output:
     set val(key), file("output_${key}/*.tsv") into sum_calc

     """
     #!/usr/bin/env Rscript

     dir.create("output_${key}")
     filepath<-"$consensusXML"

    output_file_feature <- "output_${key}/${consensusXML.baseName}_quantification.tsv"
    output_file_id <- "output_${key}/${consensusXML.baseName}_identification.tsv"
    # Open the connection
    con = file(filepath, "r")

    # read line by line until the feature headers
    # consensus header

    consensus_header <- NA
    id_header <- NA
    line = readLines(con, n = 1)
    con_pre <- FALSE
    f_line_nr <- 0

    while ( TRUE ) {
      line = readLines(con, n = 1)

      if ( length(line) == 0 ) {
        break
      }
      #data[max(grep(data,pattern = "#CONSENSUS",fixed = T))]
      line_split<-strsplit(line,split = "\t")[[1]]
      if(line_split[1]=="#CONSENSUS")
      {
        consensus_header <- c("Met_ID",line_split)
      }


      if(line_split[1]=="#PEPTIDE")
      {
        id_header <- c("Met_ID",line_split)

      }


      if(line_split[1]=="MAP")
      {
        map_id<-line_split[2]
        map_name <- basename(line_split[3])
        consensus_header<-gsub(pattern = paste("_",map_id,"\$",sep = ""),replacement = paste("_",map_name,sep=""),x = consensus_header)
      }



      if(line_split[1]=="CONSENSUS")
      {
        f_line_nr <- f_line_nr + 1
        if(f_line_nr==1)
        {
          write(paste(c(consensus_header),collapse  = "\t"),file = output_file_feature,append = TRUE)
          write(paste(c(id_header),collapse  = "\t"),file = output_file_id,append = TRUE)
        }

        write(paste(c(f_line_nr,line_split),collapse  = "\t"),file = output_file_feature,append = TRUE)
      }

      if(line_split[1]=="PEPTIDE")
      {
        write(paste(c(f_line_nr,line_split),collapse  = "\t"),file = output_file_id,append = TRUE)
      }

    }

     """
   }



   if(params.need_qc==true)
   {

     process process_qc_report  {
       label 'r'
       //label 'process_low'
       tag "$consensusXML - $key"
       publishDir "${params.outdir}/process_qc_report", mode: params.publish_dir_mode, enabled: params.publishDir_intermediate
       stageInMode 'copy'
       input:
       set val(key), file(consensusXML) from sum_calc
       each file(rmd_file) from ch_output_qc_rmd
       each file(exp_des) from experimental_design
       each file(id_in) from identification_input_qc
       output:
       set val(key), file("output_${key}/qc.html") into out

       """
       mkdir "output_${key}"
       #$consensusXML
       #$rmd_file
       #Rscript -e "rmarkdown::render('example.Rmd',params=list(args = myarg))"

       Rscript -e \\
"rmarkdown::render('$rmd_file',output_file='qc.html',params=list(input_abundance='${consensusXML[1]}',
input_id='${consensusXML[0]}',input_exp_design='$exp_des',input_internal_standards='$id_in',id_name='$params.identification_name_column',
id_rt='$params.identification_rt_column',id_use_rt=as.logical('$params.identification_use_rt'),
id_mz='$params.identification_mz_column',id_check_rt='$params.identification_use_rt',color_by='$params.covariate_qc',
interactive=as.logical('$params.interactive_plot')))"
cp qc.html output_${key}/qc.html

       """
     }
   }
 }


/*
 * Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file 'results_description.html'

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/metabolinden] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/metabolinden] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.max_multiqc_email_size)
    // On success try attach the multiqc report


    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: null, mqcMaxSize: 1 ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/metabolinden] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            mail_cmd.execute() << email_html
            log.info "[nf-core/metabolinden] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/metabolinden]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/metabolinden]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

workflow.onError {
    // Print unexpected parameters - easiest is to just rerun validation
    NfcoreSchema.validateParameters(params, json_schema, log)
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = 'hostname'.execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "${c_red}====================================================${c_reset}\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "${c_red}====================================================${c_reset}\n"
                }
            }
        }
    }
}
