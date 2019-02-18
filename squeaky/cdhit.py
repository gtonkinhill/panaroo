import subprocess
import tempfile
import os


def run_cdhit(input_file, output_file,
    id=0.95,
    n_cpu=1,
    s=0.0, # length difference cutoff (%), default 0.0
    aL=0.0, # alignment coverage for the longer sequence
    AL=99999999, # alignment coverage control for the longer sequence
    aS=0.0, # alignment coverage for the shorter sequence
    AS=99999999, # alignment coverage control for the shorter sequence
    accurate=True, # use the slower but more accurate options
    use_local=False, #whether to use local or global sequence alignment
    quiet=False):

    cmd = ("cd-hit" +
        " -T " + str(n_cpu) +
        " -i " + input_file +
        " -o " + output_file +
        " -c " + str(id) +
        " -s " + str(s) +
        " -aL " + str(aL) +
        " -AL " + str(AL) +
        " -aS " + str(aS) +
        " -AS " + str(AS) +
        " -M 0 -d 999"
        )

    if use_local:
        cmd += " -G 0"

    if accurate:
        cmd += " -g 1 -n 2"

    if not quiet:
        print("running cmd: " + cmd)
    else:
        cmd += " > /dev/null"

    subprocess.run(cmd, shell=True, check=True)

    return


def run_cdhit_est(input_file, output_file,
        id=0.95,
        n_cpu=1,
        s=0.0, # length difference cutoff (%), default 0.0
        aL=0.0, # alignment coverage for the longer sequence
        AL=99999999, # alignment coverage control for the longer sequence
        aS=0.0, # alignment coverage for the shorter sequence
        AS=99999999, # alignment coverage control for the shorter sequence
        accurate=True, # use the slower but more accurate options
        use_local=False, #whether to use local or global sequence alignment
        strand=1, # default do both +/+ & +/- alignments if set to 0, only +/+
        quiet=False):

        cmd = ("cd-hit-est" +
            " -T " + str(n_cpu) +
            " -i " + input_file +
            " -o " + output_file +
            " -c " + str(id) +
            " -s " + str(s) +
            " -aL " + str(aL) +
            " -AL " + str(AL) +
            " -aS " + str(aS) +
            " -AS " + str(AS) +
            " -r " + str(strand) +
            " -M 0 -d 999"
            )

        if use_local:
            cmd += " -G 0"

        if accurate:
            cmd += " -g 1 -n 2"

        if not quiet:
            print("running cmd: " + cmd)
        else:
            cmd += " > /dev/null"

        subprocess.run(cmd, shell=True, check=True)

        return

def cluster_nodes_cdhit(G, nodes, outdir, id=0.95, dna=False,
    s=0.0, # length difference cutoff (%), default 0.0
    aL=0.0, # alignment coverage for the longer sequence
    AL=99999999, # alignment coverage control for the longer sequence
    aS=0.0, # alignment coverage for the shorter sequence
    AS=99999999, # alignment coverage control for the shorter sequence
    accurate=True, # use the slower but more accurate options
    use_local=False, #whether to use local or global sequence alignment
    strand=1, # default do both +/+ & +/- alignments if set to 0, only +/+
    quiet=False):

    # create the files we will need
    centroids_to_nodes = {}
    temp_input_file = tempfile.NamedTemporaryFile(delete = False, dir=outdir)
    temp_input_file.close()
    temp_output_file = tempfile.NamedTemporaryFile(delete = False, dir=outdir)
    temp_output_file.close()

    with open(temp_input_file.name, 'w') as outfile:
        for node in nodes:
            outfile.write(">" + G.node[node]["centroid"] + "\n")
            centroids_to_nodes[G.node[node]["centroid"]] = node
            if dna:
                outfile.write(G.node[node]["dna"] + "\n")
            else:
                outfile.write(G.node[node]["protein"] + "\n")

    # run cd-hit
    if dna:
        run_cdhit_est(input_file=temp_input_file.name,
            output_file=temp_output_file.name,
            id=id, s=s, aL=aL, AL=AL, aS=AS, accurate=accurate,
            use_local=use_local, strand=strand, quiet=quiet)
    else:
        run_cdhit(input_file=temp_input_file.name,
            output_file=temp_output_file.name,
            id=id, s=s, aL=aL, AL=AL, aS=AS, accurate=accurate,
            use_local=use_local, quiet=quiet)

    # process the output
    clusters = []
    with open(temp_output_file.name + ".clstr", 'rU') as infile:
        c = []
        for line in infile:
            if line[0]==">":
                clusters.append(c)
                c=[]
            else:
                c.append(centroids_to_nodes[line.split(">")[1].split("...")[0]])
        clusters.append(c)
    clusters = clusters[1:]

    # remove temporary files
    os.remove(temp_input_file.name)
    os.remove(temp_output_file.name)
    os.remove(temp_output_file.name + ".clstr")

    return clusters
