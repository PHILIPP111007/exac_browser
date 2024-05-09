FROM python:3.12.3-slim
WORKDIR /usr/src/app

ADD ./exac_data/ ./exac_data

# do not create files .pyc
ENV PYTHONDONTWRITEBYTECODE 1
# do not put stdout and stderr in a buffer
ENV PYTHONUNBUFFERED 1

RUN pip install --upgrade pip

COPY ./app/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY ./app/ .
