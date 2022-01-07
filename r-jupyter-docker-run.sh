#!/bin/bash

docker run -dit -p 18888:8888 -p 54321:54321 -v ${PWD}:/home/jovyan/work --name r-jupyter -e GRANT_SUDO=yes -e JUPYTER_ENABLE_LAB=yes --user root  jupyter/r-notebook bash
