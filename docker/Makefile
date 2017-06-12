.SECONDARY: 
	
Dockerfile.base.%: base-template
	sed "s/xxx/$*/" base-template > Dockerfile.base.$*

Dockerfile.deal2lkit.%: Dockerfile.base.% deal2lkit-template
	sed "s/xxx/$*/" deal2lkit-template > Dockerfile.deal2lkit.$*

locks/base.%: Dockerfile.base.%
	docker build -t mathlab/deal2lkit-base:$* -f Dockerfile.base.$* .
	touch locks/base.$*

locks/deal2lkit.%: Dockerfile.deal2lkit.%
	docker build -t mathlab/deal2lkit:$*-debugrelease -f Dockerfile.deal2lkit.$* .
	touch locks/deal2lkit.$*

push-%: locks/base.% locks/deal2lkit.%
	docker push mathlab/deal2lkit-base:$*
	docker push mathlab/deal2lkit:$*-debugrelease
