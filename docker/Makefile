
.PHONY: base

all: base deal2lkit

base: base-8.4.2/Dockerfile base-8.5.pre.4/Dockerfile
	docker build -t mathlab/deal2lkit-base:8.4.2 base-8.4.2
	docker push mathlab/deal2lkit-base:8.4.2
	docker build -t mathlab/deal2lkit-base:8.5.pre.4 base-8.5.pre.4
	docker push mathlab/deal2lkit-base:8.5.pre.4

deal2lkit: base deal2lkit-8.4.2/Dockerfile deal2lkit-8.5.pre.4/Dockerfile
	docker build -t mathlab/deal2lkit:8.4.2-debugrelease deal2lkit-8.4.2
	docker push mathlab/deal2lkit:8.4.2-debugrelease
	docker build -t mathlab/deal2lkit:8.5.pre.4-debugrelease deal2lkit-8.5.pre.4
	docker push mathlab/deal2lkit:8.5.pre.4-debugrelease
