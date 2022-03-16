image=ljocha/gc-derivatization
port=9000

package-build:
	python3 -m build

docker-build: 
	docker build -t ${image} .

docker-run:
	docker run -p ${port}:${port} -u $(shell id -u) -w /work -v ${PWD}:/work ${image} jupyter notebook --ip 0.0.0.0 --port ${port}

docker-bash:
	docker run -ti -p ${port}:${port} -u $(shell id -u) -w /work -v ${PWD}:/work ${image} bash



