init:
	pip3 install -r requirements.txt

setup:
	setup.py

test:
	pp g3d slice "${PP_PATH}/hdf5explfiles"
	pp raw2hist 2d "${PP_PATH}/hdf5explfiles" --psv x y

.PHONY: init test