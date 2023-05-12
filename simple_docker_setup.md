# Simple docker instructions

First install docker and docker desktop

## Create image

Start from a barebones base image for what you want to create, for example: `FROM rocker/r-ver:4.1.0` or `FROM alpine:latest`

## Build container
```
docker build -t username/container_name . -f ./Dockerfile
```

## Test container 
```
docker run -it --rm=true username/container_name /bin/bash
```
* `-it` interactive flag
* `--rm=true` after we exit, this will clean up the runnining container so Docker uses less disk space.
* `username/container_name` which container to start
* `/bin/bash` tells Docker that when the container starts, we want a command line (bash) inside to run commands

or 

```
docker run username/container <tool> -h
```
to print the help command for `<tool>`
  

## Pushing to DockerHub

#### Log in

```
docker login --username=<username>
```  
Enter password

#### List docker images

```
docker images
```

#### Tag your image

```
docker tag <tag_id> <username>/<container_name>:<tagname>
```

#### Push image to repository

```
docker push <username>/<container_name>
```

Sources:
  https://chtc.cs.wisc.edu/uw-research-computing/docker-test.html. 
  https://jsta.github.io/r-docker-tutorial/04-Dockerhub.html
