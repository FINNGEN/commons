# Common scripts for data munging and miscellaneous bioinformatics

This repository curates the miscellaneous scripts that do not deserve their own repository. Everyone is encouraged to place their own little helper scripts and tools here for toher's enjoyment.

## Organization
Every helper script should have its own subdirectory, and inside that directory should be a mardown file (README.md) that describes how to use the script. Here is an example layout:
```
commons
├── helper_script_1
│   ├── README.md
│   ├── scripts
│   │   └── helper.sh
│   ├── wdl
│   │   └── helper.wdl
│   ├── docker
│   │   └── Dockerfile
│   └───data
│       └── datafile
├── helper_script_2
│   ├── README.md
│   └── scripts
│       └── helper.sh
└── README.md
```

## Script promotion
In case your script evolves and becomes larger and unwieldy to maintain in the commons repository, it should be promoted into its own repository.
