# syntax=docker/dockerfile:1

FROM python:3.9

WORKDIR .

COPY . .

EXPOSE 5101

ENTRYPOINT FLASK_APP=flask_app.py flask run --host=0.0.0.0
