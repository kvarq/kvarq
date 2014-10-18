
default_config = {
    'quality' : 13,
    'errors' : 2,
    'minimum overlap' : 25,
    'minimum readlength' : 25,
#    'stop median coverage' : 30,
    'threads' : 8,
    'spacing': 25,
}

def config_params(config, fastq):
    return dict(
            nthreads=config['threads'],
            maxerrors=config['errors'],
            minreadlength=config['minimum readlength'],
            minoverlap=config['minimum overlap'],
            Amin=fastq.Q2A(config['quality']),
            Azero=fastq.Azero
        )

