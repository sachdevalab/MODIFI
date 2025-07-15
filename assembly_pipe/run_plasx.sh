THREADS=20
PREFIX=~/borg/assembly/96plex/test/test
contigs=~/borg/assembly/96plex/96plex_p5_3/96plex_p5.hifiasm.p_ctg.rename.fa

# # Create an anvio contigs database from the fasta file
# # - The `-L 0` parameter ensures that contigs remain intact and aren't split
# anvi-gen-contigs-database -L 0 -T $THREADS --project-name $PREFIX -f $contigs -o $PREFIX.db

# # Export gene calls (including amino acid sequences) to text file
# anvi-export-gene-calls --gene-caller prodigal -c $PREFIX.db -o $PREFIX-gene-calls.txt

# # Annotate COGs
# anvi-run-ncbi-cogs -T $THREADS --cog-version COG14 --cog-data-dir /home/shuaiw/borg/paper/anvio_db/COG_2014 -c $PREFIX.db

# # Annotate Pfams
# anvi-run-pfams -T $THREADS --pfam-data-dir /home/shuaiw/borg/paper/anvio_db/Pfam_v32 -c $PREFIX.db

# # Export functions to text file
# anvi-export-functions --annotation-sources COG14_FUNCTION,Pfam -c $PREFIX.db -o $PREFIX-cogs-and-pfams.txt

# ~1 hr if THREADS=4. ~5 min if THREADS=128.
# - For faster processing, set THREADS to the number of CPU cores available.
# - This command requires a high amount of RAM (at least ~60Gb). If your machine has low RAM, then you need to set the `--splits` flag to a high number.
#   This will split the de novo families into chunks, reducing the max RAM usage. E.g. if you have only ~8Gb, we recommend setting `--splits` to 32 or higher.
plasx search_de_novo_families \
    -g $PREFIX-gene-calls.txt \
    -o $PREFIX-de-novo-families.txt \
    --threads $THREADS \
    --splits 32 \
    -db /home/shuaiw/borg/paper/anvio_db/plasx_db/PlasX_mmseqs_profiles \
    --tmp ~/borg/assembly/96plex/test/test/tmp \
    --overwrite
