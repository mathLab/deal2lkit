
.PHONY: base

all: base deal2lkit

base: base/Dockerfile
	docker build -t mathlab/deal2lkit-base base
	docker push mathlab/deal2lkit-base

deal2lkit: base deal2lkit/Dockerfile
	docker build -t mathlab/deal2lkit:master-debugrelease deal2lkit
	docker push mathlab/deal2lkit:master-debugrelease
	docker tag mathlab/deal2lkit:master-debugrelease mathlab/deal2lkit:latest
	docker push mathlab/deal2lkit:latest
