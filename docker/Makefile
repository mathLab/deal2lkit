
.PHONY: base

all: base deal2lkit

base: base-8.4.2/Dockerfile base-8.5.pre.4/Dockerfile
	docker build -t mathlab/deal2lkit-base:8.4.2 base-8.4.2
	docker push mathlab/deal2lkit-base:8.4.2
	docker build -t mathlab/deal2lkit-base:8.5.pre.4 base-8.5.pre.4
	docker push mathlab/deal2lkit-base:8.5.pre.4

deal2lkit: base deal2lkit/Dockerfile
	docker build -t mathlab/deal2lkit:master-debugrelease deal2lkit
	docker push mathlab/deal2lkit:master-debugrelease
	docker tag mathlab/deal2lkit:master-debugrelease mathlab/deal2lkit:latest
	docker push mathlab/deal2lkit:latest
