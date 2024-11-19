process CONCAT_RAW_COUNTS {
    label 'process_low'

    input:
    tuple val(meta), path(raw_counts)

    output:
    path("all_raw_counts.tsv"), emit: all_raw_counts
    path "versions.yml",    emit: versions

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import os, platform

    config_dir = "${params.base_dir}/${params.outdir}/htseq/"

    main_df = pd.DataFrame(columns=['miRNA'])

    for file in os.listdir(config_dir):
        if file.endswith(".tsv"):
            df = pd.read_csv(os.path.join(config_dir, file), sep='\\t', header=None, names=['miRNA', file[:-4]])
            main_df = main_df.merge(df, on='miRNA', how='outer')


    main_df.to_csv('all_raw_counts.tsv', sep='\\t', header=True, index=False)

    with open("versions.yml", "w") as file:
        file.write(f"${task.process}:\n")
        file.write(f"\tpython: {platform.python_version()}")
    """
}