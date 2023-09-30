# [`environment/`](environment/)

This directory contains the environment files for the project.

## [`requirements.txt`](./requirements.txt)

This file contains the list of packages used in the project. The file can be used to create a virtual environment with `pip install -r requirements.txt`.

Alternatively, the Docker image can be used to run the project.

## [`Dockerfile`](./Dockerfile)

This file contains the instructions to create a Docker image for the project. The image can be created with, e.g., `docker build -t text2model ./environment` from the root directory of the project.
