databases = {
    "pfam": "data/pfam/36.0/pfam_a.hmm",
    "superfamily": "data/superfamily/1.75/hmmlib_1.75",
    "smart": "data/smart/9.0/smart.HMMs",
    "gene3d": "data/gene3d/4.3.0/gene3d_main.hmm"
}

rule hmmsearch:
    input:
        fasta = "sample.fasta",
        hmm = databases[wildcards.database]
    output:
        "resistify_result/{database}.txt"
    threads: 4
    resources:
        mem_mb=4000,
        partition="short"
    shell:
    "hmmsearch -E 0.00001 --domtblout {output} {input.hmm} {input.fasta} > /dev/null"

rule resistify:
    input:
        expand("resistify_result/{database}.txt", database=databases.keys())
    output:
        "resistify_result/result.txt"
    threads: 1
    resources:
        mem_mb=8000,
        partition="short"
    shell:
        "resistify.py {input} > {output}"

rule all:
    input:
        "resistify_result/result.txt"