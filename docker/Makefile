VERSION = 1.2.1
TAG1 = us.gcr.io/broad-dsp-lrma/lr-pav:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-pav:latest

all: build

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)

