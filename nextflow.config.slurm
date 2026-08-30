// Nextflow overlay for the SLURM cluster.
// Apply alongside the base nextflow.config with:
// nextflow -c nextflow.config.slurm run <pipeline> ...

process {
    executor = 'slurm'
    queue = 'compute'
    
    // Default resources per process
    cpus = 4
    memory = '32 GB'
    time = '24 h'
    
}

// Enable resume by default
resume = true

// Preserve cached work so resume remains possible
cleanup = false

// Maximum parallel jobs
executor {
    queueSize = 10
    submitRateLimit = '10/1min'
}
