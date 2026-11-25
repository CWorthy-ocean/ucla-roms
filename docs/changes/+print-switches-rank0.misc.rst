Only MPI rank 0 prints the list of activated CPP switches at start-up. Before, every rank printed it, and the copies could interleave with rank 0's own start-up report.
