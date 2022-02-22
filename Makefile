image=ljocha/gc-derivatization
port=9000


build:
	docker build -t ${image} .


run:
	docker run -p ${port}:${port} -u $(shell id -u) -w /work -v ${PWD}:/work ${image} jupyter notebook --ip 0.0.0.0 --port ${port}

bash:
	docker run -ti -p ${port}:${port} -u $(shell id -u) -w /work -v ${PWD}:/work ${image} bash
