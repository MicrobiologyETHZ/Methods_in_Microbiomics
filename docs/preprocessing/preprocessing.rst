===================
Data Preprocessing
===================

-----------------------
General Considerations
-----------------------
- Adapters
- Contaminants
- Quality Filter

----------------
Example Workflow
----------------
.. image:: ../images/Data-Preprocessing.png



-----------------
Example command
-----------------

.. code-block:: console

    bbduk.sh -Xmx1G pigz=t bgzip=f usejni=t in=<forward_fastq> in2=<reverse_fastq> \
    out=stdout.fq outm=<output_adapter_matched> outs=<output_adapter_singletons>  \
    refstats=<output_adapter_stats> statscolumns=5 overwrite=t ref=<input.adapters> \
    ktrim=r k=23 mink=11 hdist=1  2 >> <log_file> | \
    bbduk.sh -Xmx1G usejni=t pigz=t bgzip=f interleaved=true overwrite=t \
    in=stdin.fq out=stdout.fq outm={output.phix_matched} outs={output.phix_singletons} \
    ref={input.phix} k=31 hdist=1 refstats={output.phix_stats} statscolumns=5 2>> {log.log} | \
    bbduk.sh -Xmx1G pigz=t bgzip=f usejni=t overwrite=t interleaved=true \
    in=stdin.fq fastawrap=10000 out1={output.fq1_clean} out2={output.fq2_clean} \
    outm={output.qc_failed} outs={output.qc_singletons} minlength={params.minlen} \
    qtrim=rl maq={params.maq} maxns=1  stats={output.qc_stats} statscolumns=5 trimq={params.trimq}  2>> {log.log};



--------------------
Other Considerations
--------------------

.. image:: ../images/Data-Preprocessing2.png


------------------
Further Reading
------------------

