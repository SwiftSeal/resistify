# Resistify

Resistify is a lightweight program designed to classify NLRs by their protein domain architecture.
I have created this program as an alternative to several similar programmes for a couple of reasons.
 
The first is to move away from using InterProScan as a dependency.
While InterProScan is a useful resource for annotating protein domains, it very feature-rich and can be challenging to set up on a new system.
It's distribution isn't well supported by conda which is an additional challenge when integrating it into automated workflows.
The actual required output of InterProScan for a program like this can be replicated by using `hmmsearch` against the core InterProScan databases.

Secondly, I've created this to be as free of dependencies as possible.
This allows Resistify to be easily distributed and quickly installed.

Thirdly, useability.
I'm aiming to design Resistify so that the process of annotation (i.e., `hmmsearch` runs) and parsing are independent.
This allows users to execute each process independently and integrate it easily into a pipeline.

I am grateful to the authors of [NLRtracker](https://github.com/slt666666/NLRtracker) for their list of NLR-related annotations which are used in this program.

## Installation

Currently, you'll need to clone this repo to get started.
A pip/conda alternative is in progress!
Will make one when it's at version 0.1.

The databases are fairly large and need to be moved to the `resistify/data/` directory first.
I've provided a script, `get_interproscan.sh`, which simplifies this process.
Run it from the base directory to get started.
I might move this functionality to a python function soon.

