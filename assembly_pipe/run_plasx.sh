THREADS=20
PREFIX=~/borg/assembly/96plex/test/test2
contigs=~/borg/assembly/96plex/96plex_p5_3/96plex_p5.hifiasm.p_ctg.rename.fa


# anvi-gen-contigs-database -L 0 -T $THREADS --project-name $PREFIX -f $contigs -o $PREFIX.db

# anvi-export-gene-calls --gene-caller prodigal -c $PREFIX.db -o $PREFIX-gene-calls.txt

# anvi-run-ncbi-cogs -T $THREADS --cog-version COG14 --cog-data-dir /home/shuaiw/borg/paper/anvio_db/COG_2014 -c $PREFIX.db

# anvi-run-pfams -T $THREADS --pfam-data-dir /home/shuaiw/borg/paper/anvio_db/Pfam_v32 -c $PREFIX.db

# anvi-export-functions --annotation-sources COG14_FUNCTION,Pfam -c $PREFIX.db -o $PREFIX-cogs-and-pfams.txt

/home/sdiamond/miniconda3/envs/plasx/bin/plasx search_de_novo_families \
    -g $PREFIX-gene-calls.txt \
    -o $PREFIX-de-novo-families.txt \
    --threads $THREADS \
    --overwrite

/home/sdiamond/miniconda3/envs/plasx/bin/plasx predict -a $PREFIX-cogs-and-pfams.txt $PREFIX-de-novo-families.txt \
    -g $PREFIX-gene-calls.txt \
    -o $PREFIX-plasmid_nr-scores.txt \
    --overwrite



