# Definitions
repository = humanpangenomics
identifier = mitohifi
version = 0.0.1
git_commit ?= $(shell git log --pretty=oneline -n 1 | cut -f1 -d " ")
name = ${repository}/${identifier}
tag = ${version}--${git_commit}

# Steps
build:
	# do the docker build
	docker build -t ${name}:${tag} .
	docker tag ${name}:${tag} ${name}:latest
	docker push ${name}:${tag}
	docker push ${name}:latest