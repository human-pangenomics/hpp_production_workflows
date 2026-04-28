# Definitions
repository = humanpangenomics
identifier = chrom_assignment
version = 1.0
git_commit ?= $(shell git log --pretty=oneline -n 1 | cut -f1 -d " ")
name = ${repository}/${identifier}
tag = ${version}--${git_commit}

# Steps
build:
	docker buildx create --name multiarch --driver docker-container --use
	
	docker buildx build --platform linux/amd64,linux/arm64 \
		-t ${name}:${tag} \
		-t ${name}:latest \
		--push .
	
	docker buildx rm multiarch

